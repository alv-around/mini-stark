pub mod fiatshamir;
pub mod prover;
pub mod verifier;

use crate::merkle::{Hash, MerklePath};
use ark_ff::PrimeField;
use digest::Digest;

pub struct FriProof<D: Digest, F: PrimeField> {
    commits: Vec<Hash<D>>,
    points: Vec<[(F, F); 3]>,
    queries: Vec<[MerklePath<D, F>; 3]>,
    quotients: Vec<Vec<F>>,
}

#[cfg(test)]
mod test {
    use crate::field::Goldilocks;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use sha2::Sha256;

    use super::prover::FriProver;

    const TWO: usize = 2;

    #[test]
    fn test_fri_new() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let _degree = poly.degree();

        let mut fri_prover = FriProver::<TWO, Sha256, _>::new(poly, blowup_factor);
        let _proof = fri_prover.prove();
    }
}
