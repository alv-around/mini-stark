use crate::air::{Constrains, Provable, TraceTable};
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
    arthur: Vec<u8>,
    trace_commit: Hash<D>,
    constrain_trace_commit: Hash<D>,
    constrain_queries: Vec<MerklePath<D, F>>,
    mixed_constrain_commit: Hash<D>,
    mixed_constrain_queries: Vec<(F, MerklePath<D, F>)>,
    fri_commit: Hash<D>,
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
    pub fn prove<T, AIR: Provable<T, F>>(
        &self,
        mut merlin: Merlin<DigestBridge<D>>,
        air: AIR,
        witness: T,
    ) -> Result<StarkProof<D, F>, Box<dyn Error>> {
        // 1. compute trace polys and commit to them
        let trace = air.trace(&witness);
        let trace_polys = trace.get_trace_polys();
        let mut trace_poly_coeffs = TraceTable::<F>::new(trace.len(), trace.width());
        for (i, poly) in trace_polys.iter().enumerate() {
            trace_poly_coeffs.add_col(i, poly.coeffs.clone());
        }
        let trace_codeword = MerkleTree::<N, D, F>::new(trace_poly_coeffs.get_data());
        let trace_commit = trace_codeword.root();
        merlin.add_bytes(&trace_commit).unwrap();

        // TODO: add the coset trick to add zk
        let lde_domain_size = self.blowup_factor * trace.len();
        let lde_domain = Radix2EvaluationDomain::new(lde_domain_size).unwrap();
        let constrains = trace.derive_constrains();
        let mut constrain_trace = TraceTable::<F>::new(lde_domain_size, constrains.len());
        for (i, poly) in constrains.get_polynomials().into_iter().enumerate() {
            let evals = poly.evaluate_over_domain(lde_domain);
            constrain_trace.add_col(i, evals.evals);
        }
        let constrain_trace_codeword = MerkleTree::<N, D, F>::new(constrain_trace.get_data());
        let constrain_trace_commit = constrain_trace_codeword.root();
        merlin.add_bytes(&constrain_trace_commit).unwrap();

        // parametric batching g = f_0 + r f_1 + r^2 f_2 + .. + r^(n-1) f_{n-1}
        let r: [F; 1] = merlin.challenge_scalars().unwrap();
        let mut mixed_constrain_poly = DensePolynomial::<_>::from_coefficients_vec(vec![F::ZERO]);
        for (i, constrain_poly) in constrains.get_polynomials().iter().enumerate() {
            mixed_constrain_poly = mixed_constrain_poly
                + DensePolynomial::<_>::from_coefficients_vec(vec![r[0].pow([i as u64])])
                    * constrain_poly;
        }
        let mixed_constrain_lde = mixed_constrain_poly
            .clone()
            .evaluate_over_domain(lde_domain);
        let mut mixed_constrain_trace = TraceTable::<F>::new(lde_domain_size, 1);
        mixed_constrain_trace.add_col(0, mixed_constrain_lde.evals);
        let mixed_constrain_codeword = MerkleTree::<N, D, F>::new(mixed_constrain_trace.get_data());
        let mixed_constrain_commit = mixed_constrain_codeword.root();

        //
        let rand_bytes: [u8; 8] = merlin.challenge_bytes().unwrap();
        let query = usize::from_le_bytes(rand_bytes) % lde_domain_size;

        // // Make the low degree test FRI
        let prover = FriProver::<N, D, _>::new(&mut merlin, mixed_constrain_poly, 2);
        let fri_commit = prover.get_initial_commit();

        // 3. Queries
        let mut constrain_queries = Vec::new();
        let mut mixed_constrain_queries = Vec::new();

        {
            // constrain queries
            for i in 0..constrains.len() {
                let leaf = constrain_trace.get_value(query, i);
                let path = constrain_trace_codeword.generate_proof(leaf).unwrap();
                constrain_queries.push(path);
            }

            // validity query
            // TODO: refactor to not use prover for these queies
            // let (leaf, path) = prover.query_first_commit(query);
            let leaf = mixed_constrain_trace.get_value(query, 0);
            let path = mixed_constrain_codeword.generate_proof(leaf).unwrap();
            mixed_constrain_queries.push((*leaf, path));
        }

        // TODO: move prove when prover initialize
        let (fri_proof, _) = prover.prove();

        let arthur = merlin.transcript().to_vec();
        Ok(StarkProof {
            arthur,
            trace_commit,
            constrain_trace_commit,
            constrain_queries,
            mixed_constrain_commit,
            mixed_constrain_queries,
            fri_commit,
            fri_proof,
        })
    }

    pub fn verify(
        &self,
        transcript: IOPattern<DigestBridge<D>>,
        constrains: Constrains<F>,
        proof: StarkProof<D, F>,
    ) -> bool {
        let StarkProof {
            arthur,
            trace_commit,
            constrain_trace_commit,
            constrain_queries,
            mixed_constrain_commit,
            mixed_constrain_queries,
            fri_commit,
            fri_proof,
        } = proof;
        let mut arthur: Arthur<'_, DigestBridge<D>, u8> = transcript.to_arthur(&arthur);
        let degree = constrains.domain.size();

        // 1. check symbolic link to quotients ??
        let _zerofier = constrains.domain.vanishing_polynomial();
        assert_eq!(arthur.next_digest().unwrap(), trace_commit);
        assert_eq!(arthur.next_digest().unwrap(), constrain_trace_commit);
        let _r: [F; 1] = arthur.challenge_scalars().unwrap();

        // 3. run queries
        // TODO: number of queries dependent of target security. For the moment one query
        let _trace_domain = Radix2EvaluationDomain::<F>::new(degree).unwrap();
        let rand_bytes: [u8; 8] = arthur.challenge_bytes().unwrap();
        let _query = usize::from_le_bytes(rand_bytes);

        let mixed_constrain_root = MerkleRoot::<D>(mixed_constrain_commit.clone());
        let (mixed_const_leaf, path) = mixed_constrain_queries[0].clone();
        assert!(mixed_constrain_root.check_proof::<N, _>(&mixed_const_leaf, path));

        // let quotient_root = MerkleRoot::<D>(proof.constrain_trace_commit);
        // for query in proof.constrain_queries.into_iter() {
        //     // assert!(quotient_root.check_proof::<N, _>(&leaf, query));
        // }

        // 2. run fri
        let fri_verifier =
            FriVerifier::<N, D, F>::new(MerkleRoot(fri_commit), degree - 1, self.blowup_factor);
        assert!(fri_verifier.verify(fri_proof, &mut arthur));

        true
    }
}
