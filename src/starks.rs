use crate::air::{Constrains, Provable, TraceTable, Verifiable};
use crate::fri::FriProof;
use crate::fri::{fiatshamir::DigestReader, prover::FriProver, verifier::FriVerifier};
use crate::merkle::{MerklePath, MerkleRoot, MerkleTree, Tree};
use crate::Hash;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain};
use digest::core_api::BlockSizeUser;
use digest::{Digest, FixedOutputReset};
use nimue::plugins::ark::FieldChallenges;
use nimue::{Arthur, ByteChallenges, ByteWriter, DigestBridge, IOPattern, Merlin};
use std::error::Error;
use std::marker::PhantomData;

pub struct StarkProof<D: Digest, F: PrimeField> {
    pub degree: usize,
    transcript: Vec<u8>,
    trace_commit: Hash<D>,
    trace_queries: Vec<MerklePath<D, F>>,
    constrain_trace_commit: Hash<D>,
    constrain_queries: Vec<MerklePath<D, F>>,
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
        mut merlin: Merlin<DigestBridge<D>>,
        air: AIR,
        witness: T,
    ) -> Result<StarkProof<D, F>, Box<dyn Error>> {
        // // compute the lde of the trace & commit to it
        let trace = air.trace(&witness);
        let degree = trace.len();
        let lde_domain_size = self.blowup_factor * degree;

        // TODO: add the coset trick to add zk
        let mut lde_trace = TraceTable::<F>::new(lde_domain_size, trace.width());
        let lde_domain = Radix2EvaluationDomain::new(lde_domain_size).unwrap();
        let column_polys = trace.interpolate_col_polys();
        for (i, poly) in column_polys.get_trace_polynomials().into_iter().enumerate() {
            let lde_column_poly = poly.evaluate_over_domain(lde_domain);
            lde_trace.add_col(i, lde_column_poly.evals);
        }
        let trace_codeword = MerkleTree::<N, D, F>::new(lde_trace.get_data());
        let trace_commit = trace_codeword.root();
        // TODO: add claim's to transcript
        // merlin.add_bytes(&claim);
        merlin.add_bytes(&trace_commit).unwrap();

        // provide trace auth paths

        // // compute the constrain lde trace  and commit to it
        let constrains = air.derive_constrains(&trace);
        let mut constrain_trace = TraceTable::<F>::new(lde_domain_size, constrains.len());
        for (i, constrain) in constrains.enumerate() {
            let constrain_lde = constrain.evaluate_over_domain(lde_domain);
            constrain_trace.add_col(i, constrain_lde.evals);
        }
        let constrain_codeword = MerkleTree::<N, D, F>::new(constrain_trace.get_data());
        let constrain_trace_commit = constrain_codeword.root();
        merlin.add_bytes(&constrain_trace_commit).unwrap();

        // // compute the validity polynomial and commit it
        // parametric batching g = f_0 + r f_1 + r^2 f_2 + .. + r^(n-1) f_{n-1}
        let r: [F; 1] = merlin.challenge_scalars().unwrap();
        let mut constrain_poly = DensePolynomial::<_>::from_coefficients_vec(vec![F::ZERO]);
        for i in 0..trace.width() {
            let trace_poly = column_polys.get_trace_poly(i);
            constrain_poly = constrain_poly
                + DensePolynomial::<_>::from_coefficients_vec(vec![r[0].pow([i as u64])])
                    * trace_poly;
        }

        let _validity_poly = constrain_poly
            .clone()
            .evaluate_over_domain(lde_domain)
            .evals;

        // // Make the low degree test FRI
        let prover = FriProver::<N, D, _>::new(&mut merlin, constrain_poly, 1);
        let fri_proof = prover.prove();

        // 3. Queries
        let mut trace_queries = Vec::new();
        let mut constrain_queries = Vec::new();
        {
            let rand_bytes: [u8; 8] = merlin.challenge_bytes().unwrap();
            let query = usize::from_le_bytes(rand_bytes) % lde_domain_size;

            // trace queries
            let leaf = lde_trace.get_value(query, 0);
            let path = trace_codeword.generate_proof(leaf).unwrap();
            trace_queries.push(path);

            // quotients queries
            let leaf = constrain_trace.get_value(query, 0);
            let path = constrain_codeword.generate_proof(leaf).unwrap();
            constrain_queries.push(path);
        }

        let transcript = merlin.transcript().to_vec();
        Ok(StarkProof {
            degree,
            transcript,
            trace_commit,
            trace_queries,
            constrain_trace_commit,
            constrain_queries,
            fri_proof,
        })
    }

    pub fn verify(
        &self,
        transcript: IOPattern<DigestBridge<D>>,
        _constrains: Constrains<F>,
        proof: StarkProof<D, F>,
    ) -> bool {
        let mut arthur: Arthur<'_, DigestBridge<D>, u8> = transcript.to_arthur(&proof.transcript);
        // 1. check symbolic link to quotients ??
        let _trace_commit = arthur.next_digest().unwrap();
        let _quotient_commit = arthur.next_digest().unwrap();

        // 2. run fri
        let fri_verifier = FriVerifier::<N, D, F>::new(
            transcript,
            MerkleRoot(proof.constrain_trace_commit.clone()),
            proof.degree,
            self.blowup_factor,
        );
        assert!(fri_verifier.verify(proof.fri_proof));

        // 3. run queries
        // TODO: number of queries dependent of target security. For the moment one query
        let trace_domain = Radix2EvaluationDomain::<F>::new(proof.degree).unwrap();
        let rand_bytes: [u8; 8] = arthur.challenge_bytes().unwrap();
        let query = usize::from_le_bytes(rand_bytes);

        let leaf = trace_domain.element(query);
        let trace_root = MerkleRoot(proof.trace_commit);
        for query in proof.trace_queries.into_iter() {
            assert!(trace_root.check_proof::<N, _>(&leaf, query));
        }

        let quotient_root = MerkleRoot(proof.constrain_trace_commit);
        for query in proof.constrain_queries.into_iter() {
            assert!(quotient_root.check_proof::<N, _>(&leaf, query));
        }

        true
    }
}
