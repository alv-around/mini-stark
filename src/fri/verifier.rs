use super::FriProof;
use crate::merkle::MerkleRoot;
use crate::util::logarithm_of_two_k;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::Digest;
use rand::Rng;
use std::iter::zip;
use std::marker::PhantomData;

pub struct FriVerifier<const TREE_WIDTH: usize, D: Digest, F> {
    degree: usize,
    blowup_factor: usize,
    alpha: usize,
    beta: Option<usize>,
    commit: MerkleRoot<D>,
    folded_commits: Vec<MerkleRoot<D>>,
    field: PhantomData<F>,
}

impl<const TREE_WIDTH: usize, D: Digest, F: PrimeField> FriVerifier<TREE_WIDTH, D, F> {
    pub fn new(commit: MerkleRoot<D>, degree: usize, blowup_factor: usize) -> Self {
        let domain_size = (degree + 1) * blowup_factor;
        let rounds = logarithm_of_two_k::<TREE_WIDTH>(domain_size).unwrap() - 1;
        let alpha = rand::rng().random_range(0..domain_size);
        Self {
            degree,
            blowup_factor,
            alpha,
            beta: None,
            commit,
            folded_commits: Vec::with_capacity(rounds),
            field: PhantomData::<F>::default(),
        }
    }

    pub fn get_alpha(&self) -> usize {
        self.alpha
    }

    pub fn commitment(&mut self, folding_commitments: Vec<MerkleRoot<D>>) -> Result<usize, String> {
        if folding_commitments.len() != self.folded_commits.capacity() {
            return Err(format!(
                "wrong configuration: rounds don't match {} vs {}",
                folding_commitments.len(),
                self.folded_commits.len()
            ));
        }
        self.folded_commits = folding_commitments;

        let rnd = rand::rng().random_range(0..(self.degree + 1) * self.blowup_factor);
        self.beta = Some(rnd);
        Ok(rnd)
    }

    pub fn verify(&self, proof: FriProof<D, F>) -> bool {
        if self.beta.is_none() {
            return false;
        } else if self.folded_commits.len() != proof.points.len() {
            return false;
        }

        let alpha = self.alpha;
        let beta = self.beta.unwrap();

        let mut index = 0usize;
        for ([(x1, y1), (x2, y2), (x3, y3)], [path1, path2, path3]) in
            zip(proof.points, proof.queries)
        {
            let quotient = DensePolynomial::from_coefficients_vec(proof.quotients[index].clone());
            let vanishing_poly = DensePolynomial::from_coefficients_vec(vec![-x1, F::ONE])
                * DensePolynomial::from_coefficients_vec(vec![-x2, F::ONE]);
            let poly = quotient * vanishing_poly;

            assert_eq!(poly.evaluate(&x1), F::ZERO);
            assert_eq!(poly.evaluate(&x2), F::ZERO);
            // TODO: check for degree of polynomial

            // FIXME: add linearity test
            // assess folding was done correctly
            // assert!(poly.degree() == 1)

            if index == 0 {
                self.commit.check_proof::<TREE_WIDTH, _>(&y1, path1);
                self.commit.check_proof::<TREE_WIDTH, _>(&y2, path2);
            } else {
                self.folded_commits[index - 1].check_proof::<TREE_WIDTH, _>(&y1, path1);
                self.folded_commits[index - 1].check_proof::<TREE_WIDTH, _>(&y2, path2);
            }
            self.folded_commits[index].check_proof::<TREE_WIDTH, _>(&y3, path3);

            index += 1;
        }

        true
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
        assert_eq!(verifier.folded_commits.capacity(), 2);
    }
}
