use crate::air::{Constrains, Provable, TraceTable, Verifiable};
use crate::fri::FriProof;
use crate::fri::{prover::FriProver, verifier::FriVerifier};
use crate::merkle::{MerkleRoot, MerkleTree, Tree};
use crate::Hash;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain};
use digest::core_api::BlockSizeUser;
use digest::{Digest, FixedOutputReset};
use nimue::{DigestBridge, IOPattern, Merlin};
use std::error::Error;
use std::marker::PhantomData;

pub struct StarkProof<D: Digest, F: PrimeField> {
    pub degree: usize,
    trace_commit: Hash<D>,
    constrain_trace_commit: Hash<D>,
    fri_proof: FriProof<D, F>,
}

pub struct Stark<const N: usize, D: Digest, F: PrimeField> {
    blowup_factor: usize,
    field: PhantomData<F>,
    digest: PhantomData<D>,
}

impl<const N: usize, D, F> Stark<N, D, F>
where
    F: PrimeField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
{
    pub fn new(blowup_factor: usize) -> Self {
        Self {
            blowup_factor,
            field: PhantomData,
            digest: PhantomData,
        }
    }
    pub fn prove<T, AIR: Provable<T, F> + Verifiable<F>>(
        &self,
        merlin: Merlin<DigestBridge<D>>,
        claim: AIR,
        witness: T,
        r: F,
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
        let mut prover = FriProver::<N, D, _>::new(merlin, constrain_poly, self.blowup_factor);
        let fri_proof = prover.prove();

        Ok(StarkProof {
            degree,
            trace_commit,
            constrain_trace_commit,
            fri_proof,
        })
    }

    pub fn verify(
        &self,
        transcript: IOPattern<DigestBridge<D>>,
        constrains: Constrains<F>,
        proof: StarkProof<D, F>,
    ) -> bool {
        let fri_verifier = FriVerifier::<N, D, F>::new(
            transcript,
            MerkleRoot(proof.constrain_trace_commit),
            proof.degree,
            self.blowup_factor,
        );

        fri_verifier.verify(proof.fri_proof)
    }
}
