use crate::air::{Constrains, Matrix, Provable};
use crate::fiatshamir::StarkIOPattern;
use crate::fri::{fiatshamir::DigestReader, prover::FriProver, verifier::FriVerifier};
use crate::fri::{FriConfig, FriProof};
use crate::merkle::{MerklePath, MerkleRoot, MerkleTree, MerkleTreeConfig, Tree};
use crate::Hash;
use ark_ff::{PrimeField, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial, Radix2EvaluationDomain};
use digest::core_api::BlockSizeUser;
use digest::{Digest, FixedOutputReset};
use nimue::plugins::ark::FieldChallenges;
use nimue::{Arthur, ByteChallenges, ByteWriter, DigestBridge, IOPattern};
use std::error::Error;
use std::marker::PhantomData;

pub struct StarkProof<D: Digest, F: PrimeField> {
    arthur: Vec<u8>,
    trace_commit: Hash<D>,
    constrain_trace_commit: Hash<D>,
    constrain_queries: Vec<MerklePath<D, F>>,
    validity_commit: Hash<D>,
    validity_queries: Vec<MerklePath<D, F>>,
    fri_commit: Hash<D>,
    fri_proof: FriProof<D, F>,
}

pub struct Stark<D, F>(StarkConfig<D, F>)
where
    F: PrimeField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone;

impl<D, F> Stark<D, F>
where
    F: PrimeField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
{
    pub fn new(config: StarkConfig<D, F>) -> Self {
        Self(config)
    }
    pub fn prove<T, AIR: Provable<T, F>>(
        &self,
        air: AIR,
        witness: T,
    ) -> Result<StarkProof<D, F>, Box<dyn Error>> {
        let mut merlin = self.0.io.to_merlin();
        // 1. compute trace polys and commit to them
        let trace = air.trace(&witness);
        let trace_domain = trace.get_domain();
        let trace_codeword =
            MerkleTree::<D, F>::new(trace.trace.get_data(), self.0.merkle_config.clone());
        let trace_commit = trace_codeword.root();
        merlin.add_bytes(&trace_commit).unwrap();

        // TODO: add the coset trick to add zk
        let lde_domain_size = self.0.blowup_factor * trace.trace.len();
        let [random_shift]: [F; 1] = merlin.challenge_scalars().unwrap();
        let lde_domain = Radix2EvaluationDomain::new(lde_domain_size)
            .unwrap()
            .get_coset(random_shift)
            .unwrap();
        let constrains = trace.derive_constrains();
        let mut constrain_trace = Matrix::<F>::new(lde_domain_size, constrains.len(), None);
        for (i, poly) in constrains.get_polynomials().into_iter().enumerate() {
            let evals = poly.evaluate_over_domain(lde_domain);
            constrain_trace.add_col(i, evals.evals);
        }
        let constrain_trace_codeword =
            MerkleTree::<D, F>::new(constrain_trace.get_data(), self.0.merkle_config.clone());
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
        let (rest, validity_poly) = mixed_constrain_poly.divide_by_vanishing_poly(trace_domain);
        assert_eq!(rest, DensePolynomial::zero());
        let validity_lde = validity_poly.clone().evaluate_over_domain(lde_domain);

        let mut validity_trace = Matrix::<F>::new(lde_domain_size, 1, None);
        validity_trace.add_col(0, validity_lde.evals);
        let validity_codeword = MerkleTree::<D, F>::new(
            validity_trace.get_data(),
            self.0.fri_config.merkle_config.clone(),
        );
        let validity_commit = validity_codeword.root();

        // 2. Queries
        let mut query_bytes = vec![0u8; 8 * self.0.query_num];
        merlin.fill_challenge_bytes(&mut query_bytes).unwrap();
        let queries = query_bytes
            .chunks_exact_mut(8)
            .map(|bytes| usize::from_le_bytes(bytes.try_into().unwrap()) % lde_domain_size)
            .collect::<Vec<_>>();

        let mut constrain_queries = Vec::new();
        let mut validity_queries = Vec::new();
        for query in queries.into_iter() {
            // constrain queries
            let leaf = constrain_trace.get_value(query, 0);
            let path = constrain_trace_codeword.generate_proof(leaf).unwrap();
            constrain_queries.push(path);

            // validity query
            let leaf = validity_trace.get_value(query, 0);
            let path = validity_codeword.generate_proof(leaf).unwrap();
            validity_queries.push(path);
        }

        // 3. Make the low degree test FRI
        let prover = FriProver::<D, _>::new(&mut merlin, validity_poly, self.0.fri_config.clone());
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

    pub fn verify(&self, constrains: Constrains<F>, proof: StarkProof<D, F>) -> bool {
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
        let mut arthur: Arthur<'_, DigestBridge<D>, u8> = self.0.io.to_arthur(&arthur);
        assert_eq!(arthur.next_digest().unwrap(), trace_commit);

        let [shift]: [F; 1] = arthur.challenge_scalars().unwrap();
        let domain = Radix2EvaluationDomain::<F>::new(self.0.degree + 1).unwrap();
        assert_eq!(arthur.next_digest().unwrap(), constrain_trace_commit);

        let [r]: [F; 1] = arthur.challenge_scalars().unwrap();

        // 2. run queries
        // TODO: number of queries dependent of target security.
        let lde_domain = Radix2EvaluationDomain::<F>::new(domain.size() * self.0.blowup_factor)
            .unwrap()
            .get_coset(shift)
            .unwrap();
        let zerofier = domain
            .vanishing_polynomial()
            .evaluate_over_domain(lde_domain);
        println!("Zerofier evals: {:?}", zerofier.evals);
        let mut query_bytes = vec![0u8; 8 * self.0.query_num];
        arthur.fill_challenge_bytes(&mut query_bytes).unwrap();
        let queries = query_bytes
            .chunks_exact_mut(8)
            .map(|bytes| usize::from_le_bytes(bytes.try_into().unwrap()) % lde_domain.size())
            .collect::<Vec<_>>();

        let validity_root = MerkleRoot::<D>(validity_commit.clone());
        let quotient_root = MerkleRoot::<D>(constrain_trace_commit);
        for (i, query) in queries.into_iter().enumerate() {
            let path = constrain_queries[i].clone();
            let mut c_x = DensePolynomial::zero();
            let mut leafs = Vec::new();
            let w_i = lde_domain.element(query);
            for (i, constrain) in constrains.get_polynomials().iter().enumerate() {
                let leaf = constrain.evaluate(&w_i);
                leafs.push(leaf);

                c_x = c_x
                    + DensePolynomial::from_coefficients_vec(vec![r.pow([i as u64])]) * constrain;
            }

            assert_eq!(leafs, path.leaf_neighbours);
            assert!(quotient_root.check_proof(path));
            let (rest, quotient) = c_x.divide_by_vanishing_poly(domain);
            assert_eq!(rest, DensePolynomial::zero());

            let path = validity_queries[i].clone();
            let evaluation = quotient.evaluate(&w_i);
            assert!(path.leaf_neighbours.contains(&evaluation));
            assert!(validity_root.check_proof(path));
        }

        // 3. run fri
        let fri_root = MerkleRoot::<D>(fri_commit);
        let fri_verifier =
            FriVerifier::<D, F>::new(fri_root, self.0.degree, self.0.fri_config.clone());
        assert!(fri_verifier.verify(fri_proof, &mut arthur));

        true
    }
}

pub struct StarkConfig<D, F>
where
    F: PrimeField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
{
    blowup_factor: usize,
    query_num: usize,
    degree: usize,
    fri_config: FriConfig<D, F>,
    merkle_config: MerkleTreeConfig<D, F>,
    io: IOPattern<DigestBridge<D>>,
}

impl<D, F> StarkConfig<D, F>
where
    F: PrimeField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
{
    pub fn new(
        blowup_factor: usize,
        query_num: usize,
        degree: usize,
        trace_columns: usize,
    ) -> Self {
        Self {
            blowup_factor,
            query_num,
            degree,
            fri_config: FriConfig {
                blowup_factor,
                merkle_config: MerkleTreeConfig {
                    leafs_per_node: 2,
                    inner_children: 2,
                    _digest: PhantomData::<D>,
                    _field: PhantomData::<F>,
                },
            },
            merkle_config: MerkleTreeConfig {
                leafs_per_node: trace_columns,
                inner_children: 2,
                _digest: PhantomData::<D>,
                _field: PhantomData::<F>,
            },
            io: StarkIOPattern::<D, F>::new_stark(degree, query_num, "🐺"),
        }
    }
}
