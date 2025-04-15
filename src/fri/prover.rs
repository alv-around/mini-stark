use super::{FriConfig, FriProof};
use crate::merkle::{MerklePath, MerkleTree, MerkleTreeConfig, Tree};
use crate::Hash;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use ark_poly::EvaluationDomain;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::core_api::BlockSizeUser;
use digest::{Digest, FixedOutputReset};
use nimue::plugins::ark::FieldChallenges;
use nimue::DigestBridge;
use nimue::{ByteChallenges, ByteWriter, Merlin};

pub struct FriProver<'a, D: Digest, F: PrimeField>
where
    F: PrimeField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
{
    round_num: usize,
    transcript: &'a mut Merlin<DigestBridge<D>>,
    rounds: Vec<FriRound<D, F>>,
    queries: usize,
}

impl<'a, D, F> FriProver<'a, D, F>
where
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    F: PrimeField,
{
    pub fn new(
        transcript: &'a mut Merlin<DigestBridge<D>>,
        poly: DensePolynomial<F>,
        config: FriConfig<D, F>,
    ) -> Self {
        let FriConfig {
            queries,
            merkle_config,
            blowup_factor,
        } = config;
        let d = poly.degree();
        let domain = Radix2EvaluationDomain::<F>::new(d * blowup_factor).unwrap();
        let domain_size = domain.size as usize;
        let round_num = domain.log_size_of_group as usize;

        // degree padding
        println!("FRI: domain:{}, d: {}", domain_size, d);
        let power_offset = (domain_size / blowup_factor) - d;
        let mut poly_offset = poly;
        if power_offset > 0 {
            let x_power =
                SparsePolynomial::<F>::from_coefficients_vec(vec![(power_offset, F::ONE)]);
            let shift = DensePolynomial::from(x_power).naive_mul(&poly_offset);
            poly_offset = poly_offset + shift;
            println!("FRI: degree padding performed");
        }

        let mut rounds = Vec::<FriRound<D, F>>::with_capacity(round_num);
        let first_round = FriRound::new(poly_offset, domain_size, merkle_config);
        rounds.push(first_round);

        Self {
            queries,
            round_num,
            rounds,
            transcript,
        }
    }

    pub fn prove(mut self) -> (FriProof<D, F>, Vec<u8>) {
        self.commit_phase();
        self.query_phase().unwrap()
    }

    pub fn commit_phase(&mut self) -> Vec<Hash<D>> {
        assert_eq!(self.rounds.len(), 1);

        let mut commits = Vec::new();
        println!("number of rounds: {}", self.round_num);
        for i in 1..self.round_num {
            let previous_round = &self.rounds[i - 1];
            let commit = previous_round.commit.root();
            self.transcript.add_bytes(&commit).unwrap();
            commits.push(commit);
            let previous_round_config = previous_round.config.clone();

            let alpha: [F; 1] = self.transcript.challenge_scalars().unwrap();
            let folded_poly = FriRound::<D, _>::split_and_fold(
                &previous_round.poly.clone(),
                alpha[0],
                previous_round_config.inner_children,
            );
            let domain_size = folded_poly.degree() + 1;

            println!("previous poly round coeffs: {:?}", previous_round.poly);
            println!("foded poly coeffs: {:?}", folded_poly);
            println!("folded poly degree:{}", folded_poly.degree());
            let round = FriRound::<D, _>::new(folded_poly, domain_size, previous_round_config);
            self.rounds.push(round);
        }

        let previous_round = &self.rounds.last().unwrap();
        let commit = previous_round.commit.root();
        self.transcript.add_bytes(&commit).unwrap();
        commits.push(commit);

        commits
    }

    pub fn query_phase(&mut self) -> Result<(FriProof<D, F>, Vec<u8>), &str> {
        let mut betas = vec![0u8; 8 * self.queries]; // usize is 64-bits
        self.transcript.fill_challenge_bytes(&mut betas).unwrap();
        let betas = betas
            .chunks_exact(8)
            .map(|a| usize::from_le_bytes(a.try_into().unwrap()))
            .collect::<Vec<usize>>();

        let mut queries = Vec::new();
        let mut points = Vec::new();
        let mut quotients = Vec::new();

        let mut rounds_iter = self.rounds.iter_mut();
        let previous_round = rounds_iter.next().unwrap();
        let mut previous_poly = &previous_round.poly;
        let mut previous_commit = &previous_round.commit;
        let mut previous_domain = &previous_round.domain;
        let config = previous_round.config.clone();
        for round in rounds_iter {
            assert_eq!(
                previous_domain.size() / config.inner_children,
                round.domain.size()
            );

            let mut round_queries = Vec::new();
            let mut round_points = Vec::new();
            let mut round_quotients = Vec::new();
            for query in &mut betas.iter() {
                let mut beta = *query;
                if beta > previous_domain.size() {
                    beta %= previous_domain.size();
                }

                let x1 = previous_domain.element(beta);
                let x2 = previous_domain.element(round.domain.size() + beta);
                let x3 = round.domain.element(beta);
                let y1 = previous_poly.evaluate(&x1);
                let y2 = previous_poly.evaluate(&x2);
                let y3 = round.poly.evaluate(&x3);
                round_points.push([(x1, y1), (x2, y2), (x3, y3)]);
                assert_eq!(x3, previous_domain.element(2 * beta));

                // quotienting
                // g(x) = ax + b
                let a = (y2 - y1) / (x2 - x1);
                let b = y1 - a * x1;
                let g = DensePolynomial::from_coefficients_vec(vec![b, a]);

                // q(x) = f(x) - g(x) / Z(x)
                let numerator = previous_poly.clone() - g;
                let vanishing_poly = FriProver::<D, F>::calculate_vanishing_poly(&[x1, x2]);
                let q = numerator / vanishing_poly;
                println!("quotient: {:?}", q);
                round_quotients.push(q.to_vec());

                // merkle commits
                let proof1 = previous_commit.generate_proof(&y1).unwrap();
                let proof2 = previous_commit.generate_proof(&y2).unwrap();
                round_queries.push([proof1, proof2]);
            }

            points.push(round_points);
            queries.push(round_queries);
            quotients.push(round_quotients);
            previous_poly = &round.poly;
            previous_commit = &round.commit;
            previous_domain = &round.domain;
            println!("another round achieved");
        }

        Ok((
            FriProof {
                points,
                queries,
                quotients,
            },
            self.transcript.transcript().to_vec(),
        ))
    }

    pub fn get_initial_commit(&self) -> Hash<D> {
        self.rounds[0].commit.root()
    }

    pub fn query_first_commit(&self, query: usize) -> (F, MerklePath<D, F>) {
        let initial_commit = self.rounds[0].clone();
        let x = initial_commit.domain.element(query);
        let leaf = initial_commit.poly.evaluate(&x);
        let path = initial_commit.commit.generate_proof(&leaf).unwrap();
        (leaf, path)
    }

    fn calculate_vanishing_poly(roots: &[F]) -> DensePolynomial<F> {
        roots
            .iter()
            .map(|i| DensePolynomial::from_coefficients_slice(&[-*i, F::ONE]))
            .reduce(|acc, e| acc * e)
            .unwrap()
    }
}

#[derive(Clone)]
struct FriRound<D: Digest, F: PrimeField> {
    poly: DensePolynomial<F>,
    commit: MerkleTree<D, F>,
    domain: Radix2EvaluationDomain<F>,
    config: MerkleTreeConfig<D, F>,
}

impl<D, F> FriRound<D, F>
where
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    F: PrimeField,
{
    fn new(poly: DensePolynomial<F>, domain_size: usize, config: MerkleTreeConfig<D, F>) -> Self {
        let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
        let config_clone = config.clone();
        let commit = FriRound::<D, F>::codeword_commit(&poly, domain, config_clone);

        Self {
            commit,
            poly,
            domain,
            config,
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
    use crate::field::Goldilocks;
    use crate::fri::fiatshamir::FriIOPattern;
    use nimue::IOPattern;
    use sha2::{self, Sha256};
    use std::marker::PhantomData;

    #[test]
    fn test_fri_prover_new() {
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let queries = 3;
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
            queries,
            merkle_config,
            blowup_factor: 2,
        };

        let fri = FriProver::new(&mut transcript, poly, config);
        assert_eq!(fri.round_num, 3);

        let _proof = fri.prove();
    }
}
