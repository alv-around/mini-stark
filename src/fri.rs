use crate::error::{ProverError, VerifierError};
use crate::fiatshamir::DigestReader;
use crate::merkle::{MerklePath, MerkleRoot, MerkleTree, MerkleTreeConfig, Tree};
use crate::Hash;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use ark_poly::EvaluationDomain;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::core_api::BlockSizeUser;
use digest::{Digest, FixedOutputReset};
use itertools::Itertools;
use log::{debug, info, trace};
use nimue::plugins::ark::FieldChallenges;
use nimue::{Arthur, ByteChallenges, ByteWriter, DigestBridge, Merlin};
use std::iter::zip;

#[derive(Clone)]
pub struct FriProof<D: Digest, F: PrimeField> {
    points: Vec<Vec<[(F, F); 3]>>,
    queries: Vec<Vec<[MerklePath<D, F>; 2]>>,
    quotients: Vec<Vec<Vec<F>>>,
}

#[derive(Clone)]
pub struct FriConfig<D: Digest, F: PrimeField> {
    pub(crate) queries: usize,
    pub(crate) merkle_config: MerkleTreeConfig<D, F>,
    pub(crate) blowup_factor: usize,
    pub(crate) rounds: usize,
}

pub struct Fri<D: Digest, F: PrimeField>(FriConfig<D, F>);

impl<D, F> Fri<D, F>
where
    F: PrimeField,
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
    ) -> Result<(FriProof<D, F>, Hash<D>, Vec<u8>), ProverError> {
        // INFO: here we don't need degree padding as we know poly is "full"
        let (commit, fri_rounds) = self.commit_phase(transcript, poly)?;
        let (proof, _) = self.query_phase(transcript, fri_rounds)?;
        // FIXME: change interface so commit is not needed
        Ok((proof, commit[0].clone(), transcript.transcript().to_vec()))
    }

    fn commit_phase(
        &self,
        transcript: &mut Merlin<DigestBridge<D>>,
        poly: DensePolynomial<F>,
    ) -> Result<(Vec<Hash<D>>, Vec<FriRound<D, F>>), ProverError> {
        info!(
            "FRI proving: commit phase - folding poly by {} {} times",
            self.0.merkle_config.inner_children, self.0.rounds
        );
        let mut commits = Vec::with_capacity(self.0.rounds);
        let mut fri_rounds = Vec::with_capacity(self.0.rounds);

        // commit to the first round without folding
        let first_round = FriRound::new(
            poly.clone(),
            self.0.blowup_factor,
            self.0.merkle_config.clone(),
        );
        let first_commit = first_round.commit.root();
        transcript.add_bytes(&first_commit)?;
        fri_rounds.push(first_round);
        commits.push(first_commit);

        // in the sucsesive rounds fold poly nomial according to alpha
        let mut previous_poly = poly;
        for _ in 1..self.0.rounds {
            let [alpha]: [F; 1] = transcript.challenge_scalars()?;
            let folded_poly = FriRound::<D, _>::split_and_fold(
                &previous_poly,
                alpha,
                self.0.merkle_config.inner_children,
            );

            trace!("previous poly round coeffs: {:?}", previous_poly);
            trace!("foded poly coeffs: {:?}", folded_poly);
            trace!("folded poly degree:{}", folded_poly.degree());
            previous_poly = folded_poly.clone();
            let round = FriRound::<D, _>::new(
                folded_poly,
                self.0.blowup_factor,
                self.0.merkle_config.clone(),
            );
            let round_commit = round.commit.root();
            transcript.add_bytes(&round_commit)?;
            fri_rounds.push(round);
            commits.push(round_commit);
        }

        Ok((commits, fri_rounds))
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

        let Transcript(commits, alphas, betas) = self.read_proof_transcript(arthur).unwrap();
        // assert_eq!(1 << commits.len(), self.domain_size);
        assert_eq!(commits.len(), self.0.rounds);
        assert_eq!(commits.len() - 1, proof.points.len());

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
                let g = DensePolynomial::from_coefficients_vec(vec![b, a]);
                assert_eq!(g.evaluate(&alphas[i]), y3);

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
        let domain_size = 1 << self.0.rounds;

        for _ in 1..self.0.rounds {
            let digest = arthur.next_digest()?;
            commits.push(MerkleRoot(digest));

            let alpha: [F; 1] = arthur.challenge_scalars()?;
            alphas.push(alpha[0]);
        }

        let digest = arthur.next_digest()?;
        commits.push(MerkleRoot(digest));

        let mut betas = vec![0u8; 8 * self.0.queries];
        arthur.fill_challenge_bytes(&mut betas)?;
        let betas = betas
            .chunks_exact(8)
            .map(|a| usize::from_le_bytes(a.try_into().unwrap()))
            .map(|a| if a > domain_size { a % domain_size } else { a })
            .collect();

        Ok(Transcript(commits, alphas, betas))
    }

    fn calculate_vanishing_poly(roots: &[F]) -> DensePolynomial<F> {
        roots
            .iter()
            .map(|i| DensePolynomial::from_coefficients_slice(&[-*i, F::ONE]))
            .reduce(|acc, e| acc * e)
            .unwrap()
    }
}

pub(super) struct Transcript<D: Digest, F: PrimeField>(Vec<MerkleRoot<D>>, Vec<F>, Vec<usize>);

#[derive(Clone)]
pub(super) struct FriRound<D: Digest, F: PrimeField> {
    poly: DensePolynomial<F>,
    commit: MerkleTree<D, F>,
    domain: Radix2EvaluationDomain<F>,
}

impl<D, F> FriRound<D, F>
where
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    F: PrimeField,
{
    fn new(poly: DensePolynomial<F>, blowup_factor: usize, config: MerkleTreeConfig<D, F>) -> Self {
        let domain_size = (poly.degree() + 1) * blowup_factor;
        let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
        let commit = FriRound::<D, F>::codeword_commit(&poly, domain, config);

        Self {
            commit,
            poly,
            domain,
        }
    }

    fn split_and_fold(
        poly: &DensePolynomial<F>,
        alpha: F,
        split_factor: usize,
    ) -> DensePolynomial<F> {
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
            .enumerate()
            .map(|(i, poly)| {
                poly.naive_mul(&DensePolynomial::from_coefficients_slice(&[
                    alpha.pow([i as u64])
                ]))
            })
            .reduce(|acc, e| acc + e)
            .unwrap()
    }

    fn codeword_commit(
        poly: &DensePolynomial<F>,
        domain: Radix2EvaluationDomain<F>,
        config: MerkleTreeConfig<D, F>,
    ) -> MerkleTree<D, F> {
        let leafs = poly.evaluate_over_domain_by_ref(domain);
        MerkleTree::<D, _>::new(&leafs.evals, config)
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

    // use super::verifier::FriVerifier;

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

    #[test]
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

        let fri = Fri::<Sha256, _>::new(config.clone());
        let _fri_proof = fri.prove(&mut transcript, poly);

        // let verifier = FriVerifier::<Sha256, Goldilocks>::new(MerkleRoot(commit), degree, config);
        // let mut arthur = io.to_arthur(&transcript);
        // assert!(verifier.verify(proof, &mut arthur).unwrap());
    }
}
