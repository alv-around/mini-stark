use crate::air::{Constrains, Matrix, Provable, TraceTable};
use crate::fri::FriProof;
use crate::fri::{fiatshamir::DigestReader, prover::FriProver, verifier::FriVerifier};
use crate::merkle::{MerklePath, MerkleRoot, MerkleTree, Tree};
use crate::Hash;
use ark_ff::{PrimeField, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial, Radix2EvaluationDomain};
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
    constrain_queries: Vec<Vec<MerklePath<D, F>>>,
    validity_commit: Hash<D>,
    validity_queries: Vec<(F, MerklePath<D, F>)>,
    fri_commit: Hash<D>,
    fri_proof: FriProof<D, F>,
}

pub struct Stark<const N: usize, D: Digest, F: PrimeField> {
    blowup_factor: usize,
    query_num: usize,
    field: PhantomData<F>,
    digest: PhantomData<D>,
}

impl<const N: usize, D, F> Stark<N, D, F>
where
    F: PrimeField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
{
    pub fn new(blowup_factor: usize, query_num: usize) -> Self {
        Self {
            blowup_factor,
            query_num,
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
        let mut trace_poly_coeffs = Matrix::<F>::new(trace.trace.len(), trace_polys.len(), None);
        for (i, poly) in trace_polys.iter().enumerate() {
            trace_poly_coeffs.add_col(i, poly.coeffs.clone());
        }
        let trace_codeword = MerkleTree::<N, D, F>::new(trace_poly_coeffs.get_data());
        let trace_commit = trace_codeword.root();
        merlin.add_bytes(&trace_commit).unwrap();

        // TODO: add the coset trick to add zk
        let lde_domain_size = self.blowup_factor * trace_poly_coeffs.len();
        let lde_domain = Radix2EvaluationDomain::new(lde_domain_size).unwrap();
        let constrains = trace.derive_constrains();
        let mut constrain_trace = Matrix::<F>::new(lde_domain_size, constrains.len(), None);
        for (i, poly) in constrains.get_polynomials().into_iter().enumerate() {
            let evals = poly.evaluate_over_domain(lde_domain);
            constrain_trace.add_col(i, evals.evals);
        }
        println!("Proving: constrain trace {:?}", constrain_trace.get_data());
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
        let (rest, validity_poly) =
            mixed_constrain_poly.divide_by_vanishing_poly(constrains.domain);
        assert_eq!(rest, DensePolynomial::zero());
        let validity_lde = validity_poly.clone().evaluate_over_domain(lde_domain);

        let mut validity_trace = Matrix::<F>::new(lde_domain_size, 1, None);
        validity_trace.add_col(0, validity_lde.evals);
        let validity_codeword = MerkleTree::<N, D, F>::new(validity_trace.get_data());
        let validity_commit = validity_codeword.root();

        // 2. Queries
        let mut query_bytes = vec![0u8; 8 * self.query_num];
        merlin.fill_challenge_bytes(&mut query_bytes).unwrap();
        let queries = query_bytes
            .chunks_exact_mut(8)
            .map(|bytes| usize::from_le_bytes(bytes.try_into().unwrap()) % lde_domain_size)
            .collect::<Vec<_>>();

        let mut constrain_queries = Vec::new();
        let mut validity_queries = Vec::new();
        for query in queries.into_iter() {
            // constrain queries
            let mut query_constrain_queries = Vec::new();
            for i in 0..constrains.len() {
                let leaf = constrain_trace.get_value(query, i);
                let path = constrain_trace_codeword.generate_proof(leaf).unwrap();
                query_constrain_queries.push(path);
            }
            constrain_queries.push(query_constrain_queries);

            // validity query
            let leaf = validity_trace.get_value(query, 0);
            let path = validity_codeword.generate_proof(leaf).unwrap();
            validity_queries.push((*leaf, path));
        }

        // 3. Make the low degree test FRI
        let prover = FriProver::<N, D, _>::new(&mut merlin, validity_poly, 2);
        let fri_commit = prover.get_initial_commit();
        let (fri_proof, _) = prover.prove();

        let arthur = merlin.transcript().to_vec();
        Ok(StarkProof {
            arthur,
            trace_commit,
            constrain_trace_commit,
            constrain_queries,
            validity_commit,
            validity_queries,
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
            validity_commit,
            validity_queries,
            fri_commit,
            fri_proof,
        } = proof;
        let mut arthur: Arthur<'_, DigestBridge<D>, u8> = transcript.to_arthur(&arthur);
        let degree = constrains.domain.size();

        // 1. check symbolic link to quotients ??
        assert_eq!(arthur.next_digest().unwrap(), trace_commit);
        assert_eq!(arthur.next_digest().unwrap(), constrain_trace_commit);
        let [r]: [F; 1] = arthur.challenge_scalars().unwrap();

        // 2. run queries
        // TODO: number of queries dependent of target security. For the moment one query
        let lde_domain = Radix2EvaluationDomain::<F>::new(degree * self.blowup_factor).unwrap();
        let zerofier = constrains
            .domain
            .vanishing_polynomial()
            .evaluate_over_domain(lde_domain);
        println!("Zerofier evals: {:?}", zerofier.evals);
        let mut query_bytes = vec![0u8; 8 * self.query_num];
        arthur.fill_challenge_bytes(&mut query_bytes).unwrap();
        let queries = query_bytes
            .chunks_exact_mut(8)
            .map(|bytes| usize::from_le_bytes(bytes.try_into().unwrap()) % lde_domain.size())
            .collect::<Vec<_>>();

        let validity_root = MerkleRoot::<D>(validity_commit.clone());
        let quotient_root = MerkleRoot::<D>(constrain_trace_commit);
        for (i, query) in queries.into_iter().enumerate() {
            let (v_x, path) = validity_queries[i].clone();
            assert!(validity_root.check_proof::<N, _>(&v_x, path));

            let mut c_x = DensePolynomial::zero();
            for (j, path) in constrain_queries[i].iter().enumerate() {
                let constrain = constrains.get_constrain_poly(j);
                let w_i = lde_domain.element(query);
                let leaf = constrain.evaluate(&w_i);
                assert!(quotient_root.check_proof::<N, _>(&leaf, path.clone()));

                c_x = c_x
                    + DensePolynomial::from_coefficients_vec(vec![r.pow([i as u64])]) * constrain;
            }
            let (rest, _quotient) = c_x.divide_by_vanishing_poly(constrains.domain);
            assert_eq!(rest, DensePolynomial::zero());
        }

        // 3. run fri
        let fri_root = MerkleRoot::<D>(fri_commit);
        let fri_verifier = FriVerifier::<N, D, F>::new(fri_root, degree - 1, self.blowup_factor);
        assert!(fri_verifier.verify(fri_proof, &mut arthur));

        true
    }
}
