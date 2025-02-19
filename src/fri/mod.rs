pub mod prover;
pub mod verifier;

use crate::merkle::MerklePath;
use ark_ff::PrimeField;
use digest::Digest;

struct FriProof<D: Digest, F: PrimeField> {
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

    use super::prover::Fri;
    use super::verifier::FriVerifier;

    const TWO: usize = 2;

    #[test]
    fn test_fri_new() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let degree = poly.degree();
        let mut fri = Fri::<TWO, Sha256, _>::new(poly, blowup_factor);
        let mut verifier = FriVerifier::<TWO, Sha256, Goldilocks>::new(
            fri.generate_commit(),
            degree,
            blowup_factor,
        );

        let alpha = verifier.get_alpha();
        let commitment = fri.commit_phase(alpha);
        assert_eq!(commitment.len(), 2);

        let beta = verifier.commitment(commitment).unwrap();
        let proof_result = fri.query_phase(beta);
        assert!(proof_result.is_ok());
        let proof = proof_result.unwrap();
        assert!(verifier.verify(proof));
    }
}
