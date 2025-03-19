use crate::air::{Constrains, Provable, TraceTable, Verifiable};
use crate::fri::FriProof;
use crate::fri::{prover::Fri, verifier::FriVerifier};
use crate::merkle::{Hash, MerkleTree, Tree};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial, Radix2EvaluationDomain};
use digest::Digest;
use std::error::Error;

struct StarkConfig<const N: usize, F, D> {
    field: F,
    blowup_factor: usize,
    digest: D,
}

struct StarkProof<D: Digest, F: PrimeField> {
    trace_commit: Hash<D>,
    constrain_trace_commit: Hash<D>,
    validity_poly_commit: Hash<D>,
    fri_proof: FriProof<D, F>,
}

struct Stark<const N: usize, D: Digest, F: PrimeField> {
    config: StarkConfig<N, F, D>,
}

impl<const N: usize, D: Digest, F: PrimeField> Stark<N, D, F> {
    pub fn prove<T, AIR: Provable<T, F> + Verifiable<F>>(
        &self,
        claim: AIR,
        witness: T,
        r: F,
        alphas: &[F],
        beta: usize,
    ) -> Result<StarkProof<D, F>, Box<dyn Error>> {
        // // compute the lde of the trace & commit to it
        let trace = claim.trace(witness);
        let lde_domain_size = self.config.blowup_factor * trace.len();
        let mut lde_trace = TraceTable::<F>::new(lde_domain_size, trace.width());

        // TODO: add the coset trick to add zk
        let lde_domain = Radix2EvaluationDomain::new(lde_domain_size).unwrap();
        let column_polys = trace.interpolate_col_polys();
        for (i, poly) in column_polys.get_trace_polynomials().into_iter().enumerate() {
            let lde_column_poly = poly.evaluate_over_domain(lde_domain);
            lde_trace.add_col(i, lde_column_poly.evals);
        }
        let trace_commit = MerkleTree::<N, D, F>::new(lde_trace.get_data()).root();

        // // compute the constrain lde trace  and commit to it
        let constrains = claim.derive_constrains(&trace);
        let mut constrain_trace = TraceTable::<F>::new(lde_domain_size, constrains.len());
        for (i, constrain) in constrains.enumerate() {
            let constrain_lde = constrain.evaluate_over_domain(lde_domain);
            constrain_trace.add_col(i, constrain_lde.evals);
        }
        let constrain_trace_commit = MerkleTree::<N, D, F>::new(constrain_trace.get_data()).root();

        // // compute the validity polynomial and commit it
        // parametric batching g = f_0 + r f_1 + r^2 f_2 + .. + r^(n-1) f_{n-1}
        let mut constrain_poly = DensePolynomial::<_>::from_coefficients_vec(vec![F::ZERO]);
        for i in 0..trace.width() {
            let trace_poly = column_polys.get_trace_poly(i);
            constrain_poly = constrain_poly
                + DensePolynomial::<_>::from_coefficients_vec(vec![r.pow([i as u64])]) * trace_poly;
        }
        let validity_poly = constrain_poly
            .clone()
            .evaluate_over_domain(lde_domain)
            .evals;
        let validity_poly_commit = MerkleTree::<N, D, F>::new(&validity_poly).root();

        // // Make the low degree test FRI
        let mut prover = Fri::<N, D, _>::new(constrain_poly, self.config.blowup_factor);
        prover.commit_phase(alphas);
        let fri_proof = prover.query_phase(beta)?;

        Ok(StarkProof {
            trace_commit,
            constrain_trace_commit,
            validity_poly_commit,
            fri_proof,
        })
    }

    pub fn verify(constrains: Constrains<F>, proof: StarkProof<D, F>, z: F) -> bool {
        let one = constrains.get_domain().element(0);
        for i in 0..constrains.get_boundary_constrain_number() {
            let boundary_constrain = constrains.get_boundary_constrain(i);
            assert_eq!(boundary_constrain.evaluate(&one), F::ZERO);
        }
        false
    }
}
