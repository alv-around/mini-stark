use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, Polynomial};
use sha2::Sha256;

use mini_starks::field::Goldilocks;
use mini_starks::fri::prover::Fri;
use mini_starks::fri::verifier::FriVerifier;

const TWO: usize = 2;

fn main() {
    let blowup_factor = 2usize;
    let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
    let poly = DensePolynomial::from_coefficients_vec(coeffs);
    let degree = poly.degree();
    let mut fri = Fri::<TWO, Sha256, _>::new(poly, blowup_factor);
    let mut verifier =
        FriVerifier::<TWO, Sha256, Goldilocks>::new(fri.generate_commit(), degree, blowup_factor);

    let alpha = verifier.get_alpha();
    let commitment = fri.commit_phase(alpha);
    assert_eq!(commitment.len(), 2);

    let beta = verifier.commitment(commitment).unwrap();
    let proof_result = fri.query_phase(beta);
    assert!(proof_result.is_ok());
    let proof = proof_result.unwrap();
    assert!(verifier.verify(proof));
}
