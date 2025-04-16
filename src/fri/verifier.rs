use super::fiatshamir::DigestReader;
use super::{FriConfig, FriProof};
use crate::merkle::MerkleRoot;
use crate::util::{ceil_log2_k, logarithm_of_two_k};
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial};
use ark_std::test_rng;
use digest::core_api::BlockSizeUser;
use digest::{Digest, FixedOutputReset};
use log::{debug, info};
use nimue::plugins::ark::FieldChallenges;
use nimue::{Arthur, ByteChallenges, DigestBridge, IOPatternError};
use std::iter::zip;
use std::marker::PhantomData;

pub struct Transcript<D: Digest, F: PrimeField>(Vec<MerkleRoot<D>>, Vec<F>, Vec<usize>);

pub struct FriVerifier<D, F>
where
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    F: PrimeField,
{
    queries: usize,
    domain_size: usize,
    rounds: usize,
    commit: MerkleRoot<D>,
    marker: PhantomData<F>,
}

impl<D, F> FriVerifier<D, F>
where
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    F: PrimeField,
{
    pub fn new(commit: MerkleRoot<D>, degree: usize, config: FriConfig<D, F>) -> Self {
        let FriConfig {
            queries,
            merkle_config,
            blowup_factor,
            ..
        } = config;
        let domain_size =
            1 << ceil_log2_k((degree + 1) * blowup_factor, merkle_config.inner_children);
        let rounds = logarithm_of_two_k(domain_size, merkle_config.inner_children).unwrap();

        info!(
            "*******\n
            FRI Verifier initialized with following config:
            commit: {:?} | degree: {} | domain size: {} | rounds: {}\n 
            *******\n",
            commit.0, degree, domain_size, rounds
        );
        Self {
            queries,
            rounds,
            domain_size,
            commit,
            marker: PhantomData::<F>,
        }
    }

    pub fn read_proof_transcript(
        &self,
        arthur: &mut Arthur<'_, DigestBridge<D>, u8>,
    ) -> Result<Transcript<D, F>, IOPatternError> {
        debug!("FRI Verifier: reading proof transcript");
        let mut commits = Vec::new();
        let mut alphas = Vec::new();

        for _ in 1..self.rounds {
            let digest = arthur.next_digest().unwrap();
            commits.push(MerkleRoot(digest));

            let alpha: [F; 1] = arthur.challenge_scalars().unwrap();
            alphas.push(alpha[0]);
        }

        let digest = arthur.next_digest().unwrap();
        commits.push(MerkleRoot(digest));

        let mut betas = vec![0u8; 8 * self.queries];
        arthur.fill_challenge_bytes(&mut betas).unwrap();
        let betas = betas
            .chunks_exact(8)
            .map(|a| usize::from_le_bytes(a.try_into().unwrap()))
            .map(|a| {
                if a > self.domain_size {
                    a % self.domain_size
                } else {
                    a
                }
            })
            .collect();

        Ok(Transcript(commits, alphas, betas))
    }

    pub fn verify(
        &self,
        proof: FriProof<D, F>,
        arthur: &mut Arthur<'_, DigestBridge<D>, u8>,
    ) -> bool {
        let Transcript(commits, alphas, betas) = self.read_proof_transcript(arthur).unwrap();
        assert_eq!(1 << commits.len(), self.domain_size);
        assert_eq!(self.commit.0, commits[0].0);
        if commits.len() != self.rounds || commits.len() - 1 != proof.points.len() {
            return false;
        }

        let domain = Radix2EvaluationDomain::<F>::new(self.domain_size).unwrap();
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
                let vanishing_poly = self.calculate_vanishing_poly(&[x1, x2, x3]);
                let total_degree = quotient.degree() + vanishing_poly.degree();
                assert!(total_degree >= 2);
                assert!(total_degree <= 1 << (self.rounds - i));
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

        true
    }

    pub fn draw_random_scalar() -> F {
        F::rand(&mut test_rng())
    }

    fn calculate_vanishing_poly(&self, roots: &[F]) -> DensePolynomial<F> {
        roots
            .iter()
            .map(|i| DensePolynomial::from_coefficients_slice(&[-*i, F::ONE]))
            .reduce(|acc, e| acc * e)
            .unwrap()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::Goldilocks;
    use crate::merkle::MerkleTreeConfig;
    use crate::Hash;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use sha2::Sha256;

    #[test]
    fn test_verifier() {
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let blowup_factor = 2;
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let degree = poly.degree();
        let hash = *Hash::<Sha256>::from_slice(&[
            196, 120, 254, 173, 12, 137, 183, 149, 64, 99, 143, 132, 76, 136, 25, 217, 164, 40, 23,
            99, 175, 146, 114, 199, 243, 150, 135, 118, 182, 5, 35, 69,
        ]);

        let merkle_config = MerkleTreeConfig {
            leafs_per_node: 2,
            inner_children: 2,
            _digest: PhantomData::<Sha256>,
            _field: PhantomData::<Goldilocks>,
        };

        let config = FriConfig {
            queries: 1,
            merkle_config,
            blowup_factor,
        };

        let verifier = FriVerifier::<Sha256, Goldilocks>::new(MerkleRoot(hash), degree, config);
        assert_eq!(verifier.rounds, 3);
    }
}
