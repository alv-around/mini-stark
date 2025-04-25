use crate::error::{ProverError, VerifierError};
use crate::fiatshamir::DigestReader;
use crate::merkle::{MerklePath, MerkleRoot, MerkleTree, MerkleTreeConfig, Tree};
use ark_ff::FftField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::core_api::BlockSizeUser;
use digest::{Digest, FixedOutputReset};
use itertools::Itertools;
use log::{debug, info, trace};
use nimue::plugins::ark::{FieldChallenges, FieldReader, FieldWriter};
use nimue::{Arthur, ByteChallenges, ByteWriter, DigestBridge, Merlin};
use std::iter::zip;

#[derive(Clone)]
pub struct FriProof<D: Digest, F: FftField> {
    points: Vec<Vec<[(F, F); 3]>>,
    queries: Vec<Vec<[MerklePath<D, F>; 2]>>,
    quotients: Vec<Vec<Vec<F>>>,
}

#[derive(Clone)]
pub struct FriConfig<D: Digest, F: FftField> {
    pub(crate) queries: usize,
    pub(crate) merkle_config: MerkleTreeConfig<D, F>,
    pub(crate) blowup_factor: usize,
    pub(crate) rounds: usize,
}

pub struct Fri<D: Digest, F: FftField>(FriConfig<D, F>);

impl<D, F> Fri<D, F>
where
    F: FftField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
{
    pub fn new(config: FriConfig<D, F>) -> Self {
        let FriConfig {
            queries,
            blowup_factor,
            rounds,
            ..
        } = config;
        info!(
            "*******\nFRI Prover initialized with following config:\nqueries: {} | blowup factor: {} | rounds: {}\n*******\n",
            queries, blowup_factor, rounds
        );
        Self(config)
    }

    pub fn prove(
        &self,
        transcript: &mut Merlin<DigestBridge<D>>,
        poly: DensePolynomial<F>,
    ) -> Result<(FriProof<D, F>, Vec<u8>), ProverError> {
        // INFO: here we don't need degree padding as we know poly is "full"
        let fri_rounds = self.commit_phase(transcript, poly)?;
        let (proof, _) = self.query_phase(transcript, fri_rounds)?;
        Ok((proof, transcript.transcript().to_vec()))
    }

    fn commit_phase(
        &self,
        transcript: &mut Merlin<DigestBridge<D>>,
        poly: DensePolynomial<F>,
    ) -> Result<Vec<FriRound<D, F>>, ProverError> {
        info!(
            "FRI proving: commit phase - folding poly by {} {} times",
            self.0.merkle_config.inner_children, self.0.rounds
        );
        let mut fri_rounds = Vec::with_capacity(self.0.rounds);
        let round_domain_size = (poly.degree() + 1) * self.0.blowup_factor;

        // commit to the first round without folding
        let mut previous_round = FriRound::new(
            poly.clone(),
            round_domain_size,
            self.0.merkle_config.clone(),
        );
        fri_rounds.push(previous_round.clone());

        // in the sucsesive rounds fold poly nomial according to alpha
        for _ in 1..self.0.rounds {
            trace!("previous poly round coeffs: {:?}", previous_round.poly);

            // deep coefficients
            let [z]: [F; 1] = transcript.challenge_scalars()?;
            let deep_coeffs = previous_round.get_deep_coeffs(z);
            let denominator = DensePolynomial::from_coefficients_vec(vec![-z, F::ONE]);
            let deep_poly = DensePolynomial::from_coefficients_slice(&deep_coeffs);
            trace!("deep poly: {:?}", deep_poly);
            transcript.add_scalars(&deep_coeffs)?;

            let [alpha]: [F; 1] = transcript.challenge_scalars()?;
            let folded_poly = previous_round.fold_poly(alpha);
            trace!("folded poly coeffs: {:?}", folded_poly);
            let deep_value =
                DensePolynomial::from_coefficients_slice(&[deep_poly.evaluate(&alpha)]);
            let round_poly = (folded_poly - deep_value) / denominator;
            trace!("round poly: {:?}", round_poly);

            let domain_size = previous_round.next_round_domain_size();
            previous_round =
                FriRound::<D, _>::new(round_poly, domain_size, self.0.merkle_config.clone());
            let round_commit = previous_round.commit.root();
            transcript.add_bytes(&round_commit)?;
            fri_rounds.push(previous_round.clone());
        }

        Ok(fri_rounds)
    }

    fn query_phase(
        &self,
        transcript: &mut Merlin<DigestBridge<D>>,
        fri_rounds: Vec<FriRound<D, F>>,
    ) -> Result<(FriProof<D, F>, Vec<u8>), ProverError> {
        info!("FRI prover: starting query phase");
        let mut betas = vec![0u8; 8 * self.0.queries]; // usize is 64-bits
        transcript.fill_challenge_bytes(&mut betas)?;
        let betas = betas
            .chunks_exact(8)
            .map(|a| usize::from_le_bytes(a.try_into().unwrap()))
            .collect::<Vec<usize>>();

        let mut queries = Vec::new();
        let mut points = Vec::new();
        let mut quotients = Vec::new();

        for (round_i, (previous, round)) in fri_rounds.into_iter().tuple_windows().enumerate() {
            println!("Prove Round {round_i}");
            assert_eq!(
                previous.domain.size() / self.0.merkle_config.inner_children,
                round.domain.size()
            );

            let mut round_queries = Vec::new();
            let mut round_points = Vec::new();
            let mut round_quotients = Vec::new();
            for query in &mut betas.iter() {
                let mut beta = *query;
                if beta > previous.domain.size() {
                    beta %= previous.domain.size();
                }

                let x1 = previous.domain.element(beta);
                let x2 = previous.domain.element(round.domain.size() + beta);
                let x3 = round.domain.element(beta);
                let y1 = previous.poly.evaluate(&x1);
                let y2 = previous.poly.evaluate(&x2);
                let y3 = round.poly.evaluate(&x3);
                round_points.push([(x1, y1), (x2, y2), (x3, y3)]);
                assert_eq!(x3, previous.domain.element(2 * beta));

                // quotienting
                // g(x) = ax + b
                let a = (y2 - y1) / (x2 - x1);
                let b = y1 - a * x1;
                let g = DensePolynomial::from_coefficients_vec(vec![b, a]);

                // q(x) = f(x) - g(x) / Z(x)
                let numerator = previous.poly.clone() - g;
                let vanishing_poly = Self::calculate_vanishing_poly(&[x1, x2]);
                let q = numerator / vanishing_poly;
                round_quotients.push(q.to_vec());

                // merkle commits
                let proof1 = previous.commit.generate_proof(&y1)?;
                let proof2 = previous.commit.generate_proof(&y2)?;
                round_queries.push([proof1, proof2]);
            }

            points.push(round_points);
            queries.push(round_queries);
            quotients.push(round_quotients);
            debug!("FRI prover - query phase: round {round_i} achieved");
        }

        Ok((
            FriProof {
                points,
                queries,
                quotients,
            },
            transcript.transcript().to_vec(),
        ))
    }

    pub fn verify(
        &self,
        proof: FriProof<D, F>,
        arthur: &mut Arthur<'_, DigestBridge<D>, u8>,
    ) -> Result<bool, VerifierError> {
        info!(
            "*******\n
            Starting FRI Verification with following config:
            rounds: {} | blowup_factor: {} \n
            *******\n",
            self.0.rounds, self.0.blowup_factor,
        );

        let Transcript(commits, alphas, betas, deep_queries, deep_polys) =
            self.read_proof_transcript(arthur).unwrap();
        assert_eq!(commits.len(), self.0.rounds - 1);
        assert_eq!(commits.len(), proof.points.len());

        let domain = Radix2EvaluationDomain::<F>::new(1 << self.0.rounds).unwrap();
        let mut prev_x3s = betas.iter().map(|a| domain.element(*a)).collect::<Vec<F>>();
        for (i, (round_points, round_queries)) in zip(proof.points, proof.queries).enumerate() {
            debug!("FRI Verifier: verification Round {}", i + 1);

            for (j, ([(x1, y1), (x2, y2), (x3, y3)], [path1, path2])) in
                zip(round_points, round_queries).enumerate()
            {
                assert_eq!(x1, prev_x3s[j]);
                assert_eq!(-x1, x2);
                assert_eq!(x1.pow([2]), x3);

                let quotient =
                    DensePolynomial::from_coefficients_vec(proof.quotients[i][j].clone());
                let vanishing_poly = Self::calculate_vanishing_poly(&[x1, x2, x3]);
                let total_degree = quotient.degree() + vanishing_poly.degree();
                assert!(total_degree >= 2);
                assert!(total_degree <= 1 << (self.0.rounds - i));
                let _ = quotient / vanishing_poly;
                // linearity test
                let a = (y2 - y1) / (x2 - x1);
                let b = y1 - a * x1;
                let deep_adjusted_y =
                    y3 * (x3 - deep_queries[i]) + deep_polys[i].evaluate(&alphas[i]);
                let g = DensePolynomial::from_coefficients_vec(vec![b, a]);
                assert_eq!(g.evaluate(&alphas[i]), deep_adjusted_y);

                assert!(path1.leaf_neighbours.contains(&y1));
                commits[i].check_proof::<_>(path1);
                assert!(path2.leaf_neighbours.contains(&y2));
                commits[i].check_proof::<_>(path2);
                prev_x3s[j] = x3;
            }
        }

        Ok(true)
    }

    fn read_proof_transcript(
        &self,
        arthur: &mut Arthur<'_, DigestBridge<D>, u8>,
    ) -> Result<Transcript<D, F>, VerifierError> {
        debug!("FRI Verifier: reading proof transcript");
        let mut commits = Vec::new();
        let mut alphas = Vec::new();
        let mut deep_queries = Vec::new();
        let mut deep_polys = Vec::new();
        let domain_size = 1 << self.0.rounds;

        for _ in 1..self.0.rounds {
            // DEEP coefficinets
            let [z]: [F; 1] = arthur.challenge_scalars()?;
            deep_queries.push(z);
            let b_coeffs: [F; 2] = arthur.next_scalars()?;
            let b = DensePolynomial::from_coefficients_slice(&b_coeffs);
            deep_polys.push(b);

            let [alpha]: [F; 1] = arthur.challenge_scalars()?;
            alphas.push(alpha);
            let digest = arthur.next_digest()?;
            commits.push(MerkleRoot(digest));
        }

        let mut betas = vec![0u8; 8 * self.0.queries];
        arthur.fill_challenge_bytes(&mut betas)?;
        let betas = betas
            .chunks_exact(8)
            .map(|a| usize::from_le_bytes(a.try_into().unwrap()))
            .map(|a| if a > domain_size { a % domain_size } else { a })
            .collect();

        Ok(Transcript(commits, alphas, betas, deep_queries, deep_polys))
    }

    fn calculate_vanishing_poly(roots: &[F]) -> DensePolynomial<F> {
        roots
            .iter()
            .map(|i| DensePolynomial::from_coefficients_slice(&[-*i, F::ONE]))
            .reduce(|acc, e| acc * e)
            .unwrap()
    }
}

pub(super) struct Transcript<D: Digest, F: FftField>(
    Vec<MerkleRoot<D>>,      // round commits
    Vec<F>,                  // alphas
    Vec<usize>,              // betas
    Vec<F>,                  // deep queries
    Vec<DensePolynomial<F>>, // deep polys
);

#[derive(Clone)]
pub(super) struct FriRound<D: Digest, F: FftField> {
    poly: DensePolynomial<F>,
    commit: MerkleTree<D, F>,
    domain: Radix2EvaluationDomain<F>,
    split_factor: usize,
    splited_polys: Vec<DensePolynomial<F>>,
}

impl<D, F> FriRound<D, F>
where
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    F: FftField,
{
    fn new(poly: DensePolynomial<F>, domain_size: usize, config: MerkleTreeConfig<D, F>) -> Self {
        let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
        let split_factor = config.inner_children;
        let splited_polys = Self::split_poly(&poly, split_factor);
        let commit = FriRound::<D, F>::codeword_commit(&poly, domain, config);

        Self {
            poly,
            commit,
            domain,
            split_factor,
            splited_polys,
        }
    }

    fn split_poly(poly: &DensePolynomial<F>, split_factor: usize) -> Vec<DensePolynomial<F>> {
        let mut coeffs_vectors = Vec::<Vec<F>>::new();
        for _ in 0..split_factor {
            coeffs_vectors.push(Vec::<F>::new());
        }
        for (i, element) in poly.coeffs().iter().enumerate() {
            let index = i % split_factor;
            coeffs_vectors[index].push(*element);
        }

        coeffs_vectors
            .iter()
            .map(|coeffs: &Vec<F>| DensePolynomial::from_coefficients_slice(coeffs))
            .collect::<Vec<DensePolynomial<F>>>()
    }

    fn codeword_commit(
        poly: &DensePolynomial<F>,
        domain: Radix2EvaluationDomain<F>,
        config: MerkleTreeConfig<D, F>,
    ) -> MerkleTree<D, F> {
        let leafs = poly.evaluate_over_domain_by_ref(domain);
        MerkleTree::<D, _>::new(&leafs.evals, config)
    }

    fn get_deep_coeffs(&self, z: F) -> [F; 2] {
        [
            self.splited_polys[0].evaluate(&z),
            self.splited_polys[1].evaluate(&z),
        ]
    }

    fn fold_poly(&self, alpha: F) -> DensePolynomial<F> {
        self.splited_polys
            .iter()
            .enumerate()
            .map(|(i, poly)| {
                poly.naive_mul(&DensePolynomial::from_coefficients_slice(&[
                    alpha.pow([i as u64])
                ]))
            })
            .reduce(|acc, e| acc + e)
            .unwrap()
    }

    fn next_round_domain_size(&self) -> usize {
        self.domain.size() / self.split_factor
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::fiatshamir::FriIOPattern;
    use crate::field::Goldilocks;
    use crate::merkle::MerkleTreeConfig;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use nimue::DigestBridge;
    use nimue::IOPattern;
    use sha2::{self, Sha256};
    use std::marker::PhantomData;

    #[test]
    #[ignore]
    fn test_split_fold_poly() {}

    #[test]
    fn test_fri_prover_new() {
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let queries = 3;
        let rounds = 3;
        let io: IOPattern<DigestBridge<Sha256>> =
            FriIOPattern::<_, Goldilocks>::new_fri("🍟", 3, queries);
        let mut transcript = io.to_merlin();

        let merkle_config = MerkleTreeConfig {
            leafs_per_node: 2,
            inner_children: 2,
            _digest: PhantomData::<Sha256>,
            _field: PhantomData::<Goldilocks>,
        };

        let config = FriConfig {
            rounds,
            queries,
            merkle_config,
            blowup_factor: 2,
        };

        let fri = Fri::new(config);
        assert_eq!(fri.0.rounds, 3);

        let _fri_proof = fri.prove(&mut transcript, poly);
    }

    #[test_log::test]
    fn test_fri_new() {
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let queries = 1;
        let rounds = 3;
        let io: IOPattern<DigestBridge<Sha256>> =
            FriIOPattern::<_, Goldilocks>::new_fri("🍟", rounds, 2);
        let mut transcript = io.to_merlin();

        let merkle_config = MerkleTreeConfig {
            leafs_per_node: 2,
            inner_children: 2,
            _digest: PhantomData::<Sha256>,
            _field: PhantomData::<Goldilocks>,
        };

        let config = FriConfig {
            queries,
            rounds,
            merkle_config,
            blowup_factor: 2,
        };

        let fri = Fri::<Sha256, _>::new(config);
        let (proof, transcript) = fri.prove(&mut transcript, poly).unwrap();
        let mut arthur = io.to_arthur(&transcript);
        assert!(fri.verify(proof, &mut arthur).unwrap());
    }
}
