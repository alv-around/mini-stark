use super::FriProof;
use crate::merkle::MerkleRoot;
use crate::util::logarithm_of_two_k;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial};
use ark_std::test_rng;
use digest::Digest;
use rand::Rng;
use std::iter::zip;

pub struct FriVerifier<const TREE_WIDTH: usize, D: Digest, F> {
    degree: usize,
    blowup_factor: usize,
    alphas: Vec<F>,
    beta: Option<usize>,
    commits: Vec<MerkleRoot<D>>,
}

impl<const TREE_WIDTH: usize, D: Digest, F: PrimeField> FriVerifier<TREE_WIDTH, D, F> {
    pub fn new(commit: MerkleRoot<D>, degree: usize, blowup_factor: usize) -> Self {
        let domain_size = (degree + 1) * blowup_factor;
        let rounds = logarithm_of_two_k::<TREE_WIDTH>(domain_size).unwrap();
        // TODO: replace rng
        let mut alphas = Vec::with_capacity(rounds);
        for _ in 1..rounds {
            alphas.push(F::rand(&mut test_rng()));
        }
        let mut commits = Vec::with_capacity(rounds);
        commits.push(commit);
        Self {
            degree,
            blowup_factor,
            alphas,
            beta: None,
            commits,
        }
    }

    pub fn get_alpha(&self) -> &[F] {
        &self.alphas
    }

    pub fn commitment(&mut self, folding_commitments: Vec<MerkleRoot<D>>) -> Result<usize, String> {
        if folding_commitments.len() + 1 != self.commits.capacity() {
            return Err(format!(
                "wrong configuration: rounds don't match {} vs {}",
                folding_commitments.len(),
                self.commits.len()
            ));
        }
        for commit in folding_commitments {
            self.commits.push(commit);
        }

        let rnd = rand::rng().random_range(0..(self.degree + 1) * self.blowup_factor);
        self.beta = Some(rnd);
        Ok(rnd)
    }

    pub fn verify(&self, proof: FriProof<D, F>) -> bool {
        if self.beta.is_none() || self.commits.len() - 1 != proof.points.len() {
            return false;
        }

        let domain_size = (self.degree + 1) * self.blowup_factor;
        let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
        let mut prev_x3 = domain.element(self.beta.unwrap());
        for (i, ([(x1, y1), (x2, y2), (x3, y3)], [path1, path2, path3])) in
            zip(proof.points, proof.queries).enumerate()
        {
            println!("Round {}", i);
            assert_eq!(x1, prev_x3);
            assert_eq!(-x1, x2);
            assert_eq!(x1.pow([2]), x3);

            let quotient = DensePolynomial::from_coefficients_vec(proof.quotients[i].clone());
            let vanishing_poly = self.calculate_vanishing_poly(&[x1, x2, x3]);
            let poly = quotient * vanishing_poly;
            assert_eq!(poly.evaluate(&x1), F::ZERO);
            assert_eq!(poly.evaluate(&x2), F::ZERO);
            // FIXME: degree test
            // assert_eq!(poly.degree(), self.degree / (TREE_WIDTH ^ (i + 1)));

            // linearity test
            let a = (y2 - y1) / (x2 - x1);
            let b = y1 - a * x1;
            let g = DensePolynomial::from_coefficients_vec(vec![b, a]);
            assert_eq!(g.evaluate(&self.alphas[i]), y3);

            self.commits[i].check_proof::<TREE_WIDTH, _>(&y1, path1);
            self.commits[i].check_proof::<TREE_WIDTH, _>(&y2, path2);
            self.commits[i].check_proof::<TREE_WIDTH, _>(&y3, path3);

            prev_x3 = x3;
        }

        true
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
        let verifier =
            FriVerifier::<TWO, Sha256, Goldilocks>::new(MerkleRoot(hash), degree, blowup_factor);
        assert_eq!(verifier.commits.capacity(), 3);
    }
}
