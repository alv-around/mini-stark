use crate::fri::{prover::FriProver, verifier::FriVerifier};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::Digest;

struct STARK<const N: usize, D: Digest, F: PrimeField> {
    polys: Vec<DensePolynomial<F>>,
    prover: FriProver<N, D, F>,
}

impl<const N: usize, D: Digest, F: PrimeField> STARK<N, D, F> {
    pub fn new(polys: Vec<DensePolynomial<F>>, blowup_factor: usize) -> Self {
        let batch_randomness = FriVerifier::<N, D, F>::draw_random_scalar();
        let mut batch_poly = DensePolynomial::<_>::from_coefficients_vec(vec![F::ZERO]);
        // parametric batching g = f_0 + r f_1 + r^2 f_2 + .. + r^(n-1) f_{n-1}
        for (i, f) in polys.iter().enumerate() {
            batch_poly = batch_poly
                + DensePolynomial::<_>::from_coefficients_vec(vec![
                    batch_randomness.pow([i as u64])
                ]) * f;
        }
        let max_d = polys.iter().map(|x| x.degree()).max();
        let prover = FriProver::<N, D, _>::new(batch_poly, blowup_factor);
        Self { polys, prover }
    }
}
