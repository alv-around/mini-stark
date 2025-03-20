use crate::air::{Constrains, Provable, TraceTable, Verifiable};
use crate::fri::FriProof;
use crate::fri::{prover::Fri, verifier::FriVerifier};
use crate::merkle::{Hash, MerkleRoot, MerkleTree, Tree};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain};
use digest::Digest;
use std::error::Error;
use std::marker::PhantomData;

pub struct StarkProof<D: Digest, F: PrimeField> {
    pub degree: usize,
    trace_commit: Hash<D>,
    constrain_trace_commit: Hash<D>,
    validity_poly_commit: MerkleRoot<D>,
    fri_commits: Vec<MerkleRoot<D>>,
    fri_proof: FriProof<D, F>,
}

pub struct Stark<const N: usize, D: Digest, F: PrimeField> {
    blowup_factor: usize,
    field: PhantomData<F>,
    digest: PhantomData<D>,
}

impl<const N: usize, D: Digest, F: PrimeField> Stark<N, D, F> {
    pub fn new(blowup_factor: usize) -> Self {
        Self {
            blowup_factor,
            field: PhantomData,
            digest: PhantomData,
        }
    }
    pub fn prove<T, AIR: Provable<T, F> + Verifiable<F>>(
        &self,
        claim: AIR,
        witness: T,
        r: F,
        alphas: &[F],
        beta: usize,
    ) -> Result<StarkProof<D, F>, Box<dyn Error>> {
        // // compute the lde of the trace & commit to it
        let trace = claim.trace(&witness);
        let degree = trace.len();
        let lde_domain_size = self.blowup_factor * degree;
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

        // // Make the low degree test FRI
        let mut prover = Fri::<N, D, _>::new(constrain_poly, self.blowup_factor);
        let validity_poly_commit = prover.generate_commit();
        let fri_commits = prover.commit_phase(alphas);
        let fri_proof = prover.query_phase(beta)?;

        Ok(StarkProof {
            degree,
            trace_commit,
            constrain_trace_commit,
            validity_poly_commit,
            fri_commits,
            fri_proof,
        })
    }

    pub fn verify(
        &self,
        constrains: Constrains<F>,
        proof: StarkProof<D, F>,
        alphas: Vec<F>,
        beta: usize,
    ) -> bool {
        let mut commits = Vec::new();
        commits.push(proof.validity_poly_commit);
        for commit in proof.fri_commits {
            commits.push(commit);
        }
        let fri_verifier = FriVerifier::<N, D, F>::new_with_config(
            proof.degree,
            self.blowup_factor,
            commits,
            alphas,
            beta,
        );

        fri_verifier.verify(proof.fri_proof)
    }
}
