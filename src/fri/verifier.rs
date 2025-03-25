use super::FriProof;
use crate::merkle::{Hash, MerkleRoot};
use crate::util::{ceil_log2_k, logarithm_of_two_k};
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial};
use ark_std::test_rng;
use digest::core_api::BlockSizeUser;
use digest::{generic_array::GenericArray, Digest, FixedOutputReset, OutputSizeUser};
use nimue::plugins::ark::FieldChallenges;
use nimue::{ByteChallenges, ByteReader, DigestBridge, IOPattern, IOPatternError};
use std::iter::zip;
use std::marker::PhantomData;

pub struct FriVerifier<const TREE_WIDTH: usize, D, F>
where
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    F: PrimeField,
{
    transcript: IOPattern<DigestBridge<D>>,
    degree: usize,
    domain_size: usize,
    blowup_factor: usize,
    rounds: usize,
    commit: MerkleRoot<D>,
    marker: PhantomData<F>,
}

impl<const TREE_WIDTH: usize, D, F> FriVerifier<TREE_WIDTH, D, F>
where
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    F: PrimeField,
{
    pub fn new(
        transcript: IOPattern<DigestBridge<D>>,
        commit: MerkleRoot<D>,
        degree: usize,
        blowup_factor: usize,
    ) -> Self {
        let domain_size = 1 << ceil_log2_k::<TREE_WIDTH>((degree + 1) * blowup_factor);
        let rounds = logarithm_of_two_k::<TREE_WIDTH>(domain_size).unwrap();

        Self {
            transcript,
            degree,
            rounds,
            domain_size,
            blowup_factor,
            commit,
            marker: PhantomData::<F>,
        }
    }

    pub fn read_proof_transcript(
        &self,
        transcript: &[u8],
    ) -> Result<(Vec<MerkleRoot<D>>, Vec<F>, usize), IOPatternError> {
        let mut arthur = self.transcript.to_arthur(transcript);
        let mut commits = Vec::new();
        let mut alphas = Vec::new();

        for _ in 0..self.rounds - 1 {
            //TODO: move next lines to DigestReader trait
            let mut digest_bytes = vec![0u8; <D as OutputSizeUser>::output_size()];
            arthur.fill_next_bytes(&mut digest_bytes).unwrap();
            let digest = GenericArray::from_exact_iter(digest_bytes).unwrap();
            commits.push(MerkleRoot(digest));

            let alpha: [F; 1] = arthur.challenge_scalars().unwrap();
            alphas.push(alpha[0]);
        }

        let mut digest_bytes = vec![0u8; <D as OutputSizeUser>::output_size()];
        arthur.fill_next_bytes(&mut digest_bytes).unwrap();
        let digest = GenericArray::from_exact_iter(digest_bytes).unwrap();
        commits.push(MerkleRoot(digest));

        let beta = arthur.challenge_bytes().unwrap();
        let padded_beta = usize::from_le_bytes([beta, [0u8; 4]].concat().try_into().unwrap());

        Ok((commits, alphas, padded_beta))
    }

    pub fn verify(&self, proof: FriProof<D, F>) -> bool {
        let (commits, alphas, beta) = self.read_proof_transcript(proof.transcript).unwrap();

        if commits.len() != self.rounds || commits[0].0 != self.commit.0 {
            return false;
        } else if commits.len() - 1 != proof.points.len() {
            return false;
        }

        let domain = Radix2EvaluationDomain::<F>::new(self.domain_size).unwrap();
        let mut prev_x3 = domain.element(beta);
        for (i, ([(x1, y1), (x2, y2), (x3, y3)], [path1, path2, path3])) in
            zip(proof.points, proof.queries).enumerate()
        {
            println!("Round {}", i);
            assert_eq!(x1, prev_x3);
            assert_eq!(-x1, x2);
            assert_eq!(x1.pow([2]), x3);

            let quotient = DensePolynomial::from_coefficients_vec(proof.quotients[i].clone());
            let vanishing_poly = self.calculate_vanishing_poly(&[x1, x2, x3]);
            let total_degree = quotient.degree() + vanishing_poly.degree();
            assert!(total_degree >= 2);
            assert!(total_degree <= 1 << (self.rounds - i));
            // FIXME: find a better solution to divide by vanishing poly
            let _ = quotient / vanishing_poly;

            // linearity test
            let a = (y2 - y1) / (x2 - x1);
            let b = y1 - a * x1;
            let g = DensePolynomial::from_coefficients_vec(vec![b, a]);
            assert_eq!(g.evaluate(&alphas[i]), y3);

            commits[i].check_proof::<TREE_WIDTH, _>(&y1, path1);
            commits[i].check_proof::<TREE_WIDTH, _>(&y2, path2);
            commits[i].check_proof::<TREE_WIDTH, _>(&y3, path3);

            prev_x3 = x3;
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
    use crate::fri::fiatshamir::FriIOPattern;
    use crate::merkle::Hash;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use sha2::Sha256;

    const TWO: usize = 2;

    #[test]
    fn test_verifier() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let degree = poly.degree();
        let hash = *Hash::<Sha256>::from_slice(&[
            196, 120, 254, 173, 12, 137, 183, 149, 64, 99, 143, 132, 76, 136, 25, 217, 164, 40, 23,
            99, 175, 146, 114, 199, 243, 150, 135, 118, 182, 5, 35, 69,
        ]);
        let transcript: IOPattern<DigestBridge<Sha256>> =
            FriIOPattern::<_, Goldilocks>::new_fri("🍟", 3);
        let verifier = FriVerifier::<TWO, Sha256, Goldilocks>::new(
            transcript,
            MerkleRoot(hash),
            degree,
            blowup_factor,
        );
        assert_eq!(verifier.rounds, 3);
    }
}
