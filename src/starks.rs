use crate::air::{Constrains, Matrix, Provable};
use crate::error::{ProverError, VerifierError};
use crate::fiatshamir::{DigestIOWritter, StarkIOPattern};
use crate::fiatshamir::{DigestReader, FriIOPattern};
use crate::fri::{Fri, FriConfig, FriProof};
use crate::merkle::{MerkleTree, MerkleTreeConfig, Tree};
use crate::util::ceil_log2_k;
use crate::Hash;
use ark_ff::{FftField, PrimeField, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial, Radix2EvaluationDomain};
use digest::core_api::BlockSizeUser;
use digest::{Digest, FixedOutputReset};
use log::{debug, error, info};
use nimue::plugins::ark::{FieldChallenges, FieldIOPattern};
use nimue::{Arthur, ByteWriter, DigestBridge, IOPattern};
use std::iter::zip;
use std::marker::PhantomData;

pub struct StarkProof<D: Digest, F: FftField> {
    arthur: Vec<u8>,
    trace_commit: Hash<D>,
    constrain_trace_commit: Hash<D>,
    constrain_queries: Vec<Vec<F>>,
    validity_queries: Vec<F>,
    fri_proof: FriProof<D, F>,
}

pub struct Stark<D, F>(StarkConfig<D, F>)
where
    F: FftField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone;

impl<D, F> Stark<D, F>
where
    F: FftField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
{
    pub fn new(config: StarkConfig<D, F>) -> Self {
        let StarkConfig {
            steps,
            security_bits,
            blowup_factor,
            rounds,
            ..
        } = config;
        info!(
            "\n
            --------- \n 
            New STARK with following configuration:\n
            trace length: {} | security bits: {} | blowup factor: {} | rounds: {}\n
            --------- \n",
            steps, security_bits, blowup_factor, rounds
        );
        Self(config)
    }

    pub fn prove<T, AIR: Provable<T, F>>(
        &self,
        air: AIR,
        witness: T,
    ) -> Result<StarkProof<D, F>, ProverError> {
        let mut merlin = self.0.io.to_merlin();
        info!("Proving: start proving...");

        // 1.1 compute trace and commit to trace
        let trace = air.trace(&witness);
        let trace_domain = trace.get_domain();
        let trace_codeword =
            MerkleTree::<D, F>::new(trace.trace.get_data(), self.0.merkle_config.clone());
        let trace_commit = trace_codeword.root();
        merlin.add_bytes(&trace_commit)?;
        debug!(
            "Proving: 1.1 original trace ({:?}) committed",
            trace_commit.to_ascii_lowercase()
        );

        // 1.2 calculate low degree extension of the original trace
        let lde_domain_size = self.0.blowup_factor * trace_domain.size();
        let [random_shift]: [F; 1] = merlin.challenge_scalars()?;
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
        merlin.add_bytes(&constrain_trace_commit)?;
        debug!(
            "domain size: {} | lde domain size: {} | blowup factor: {}",
            trace_domain.size(),
            lde_domain_size,
            self.0.blowup_factor
        );
        debug!(
            "Proving: 1.2 constrain trace ({:?}) committed",
            constrain_trace_commit.to_ascii_lowercase()
        );

        // 1.3 mix constrains polynomial into the validity polynomial
        let r: [F; 1] = merlin.challenge_scalars()?;
        debug!("random variable for mixing r: {:?}", r);
        let mut mixed_constrain_poly = DensePolynomial::<_>::from_coefficients_vec(vec![F::ZERO]);
        // parametric batching g = f_0 + r f_1 + r^2 f_2 + .. + r^(n-1) f_{n-1}
        for (i, constrain_poly) in constrains.get_polynomials().iter().enumerate() {
            mixed_constrain_poly = mixed_constrain_poly
                + DensePolynomial::<_>::from_coefficients_vec(vec![r[0].pow([i as u64])])
                    * constrain_poly;
        }
        let (rest, validity_poly) = mixed_constrain_poly.divide_by_vanishing_poly(trace_domain);
        assert_eq!(rest, DensePolynomial::zero());
        debug!("Proving: 1.3 validity poly from mixing constrains",);

        // 2. DEEP-ALI: receive random queries from verifier and prove the queries correspond to the
        //    previous work.
        let mut queries = vec![F::ZERO; self.0.constrain_queries];
        merlin.fill_challenge_scalars(&mut queries)?;
        debug!(
            "Proving: 2. draw {} trace queries",
            self.0.constrain_queries
        );

        // 2.1 Link the committed constrain trace to the validity trace
        let mut constrain_queries = Vec::new();
        let mut validity_queries = Vec::new();
        for query in queries.into_iter() {
            // constrain queries
            let constrain_query: Vec<F> = constrains
                .get_polynomials()
                .iter()
                .map(|constrain| constrain.evaluate(&query))
                .collect();
            constrain_queries.push(constrain_query);

            // validity query
            let validity_query = validity_poly.evaluate(&query);
            validity_queries.push(validity_query);
        }
        debug!("Proving: 2.1 queries to committed traces provided");

        // 3. DEEP-IOPP: Make the low degree test FRI on the validity polynomial
        let fri = Fri::<D, _>::new(self.0.fri_config.clone());
        let (fri_proof, _) = fri.prove(&mut merlin, validity_poly)?;
        debug!("Proving: 3 FRI test proved");

        info!("Proving: Finished successfully!");
        let arthur = merlin.transcript().to_vec();
        Ok(StarkProof {
            arthur,
            trace_commit,
            constrain_trace_commit,
            constrain_queries,
            validity_queries,
            fri_proof,
        })
    }

    pub fn verify(
        &self,
        constrains: Constrains<F>,
        proof: StarkProof<D, F>,
    ) -> Result<bool, VerifierError> {
        info!("Verification: start verification...");
        let StarkProof {
            arthur,
            trace_commit,
            constrain_trace_commit,
            constrain_queries,
            validity_queries,
            fri_proof,
        } = proof;
        // 1. assert proof commits match transcript and calculate coset
        let mut arthur: Arthur<'_, DigestBridge<D>, u8> = self.0.io.to_arthur(&arthur);
        assert_eq!(arthur.next_digest()?, trace_commit);

        let [_shift]: [F; 1] = arthur.challenge_scalars()?;
        let domain = Radix2EvaluationDomain::<F>::new(self.0.degree + 1).unwrap();
        assert_eq!(arthur.next_digest()?, constrain_trace_commit);

        let [r]: [F; 1] = arthur.challenge_scalars()?;
        debug!("Verification: 1. proof commits match transcript");

        // 2. build validity polynomial assert it matches the the query values provided in
        //    the proof
        let mut queries = vec![F::ZERO; self.0.constrain_queries];
        arthur.fill_challenge_scalars(&mut queries)?;
        debug!("Verification: 2.1 queries from transcript retrieved");

        // 2.2 assert that the queries provided from constrain trace match commit
        // build with this queries the validity polynomial
        for (query, (constrain_query, validity_query)) in
            zip(queries, zip(constrain_queries, validity_queries))
        {
            let mut c_x = DensePolynomial::zero();
            for (i, (constrain, constrain_eval)) in
                zip(constrains.get_polynomials(), constrain_query).enumerate()
            {
                assert_eq!(constrain.evaluate(&query), constrain_eval);
                c_x = c_x
                    + DensePolynomial::from_coefficients_vec(vec![r.pow([i as u64])]) * constrain;
            }

            let (rest, quotient) = c_x.divide_by_vanishing_poly(domain);
            assert_eq!(rest, DensePolynomial::zero());

            let evaluation = quotient.evaluate(&query);
            assert_eq!(evaluation, validity_query);
        }
        debug!("Verification: 2.2 linking between validity and constrain polynomials successfull");

        // 3. run fri
        let fri_verifier = Fri::<D, F>::new(self.0.fri_config.clone());
        assert!(fri_verifier.verify(fri_proof, &mut arthur).unwrap());
        debug!("Verification: 3. FRI verification passed");

        info!("Verification: proof verification completed");
        Ok(true)
    }
}

pub struct StarkConfig<D, F>
where
    F: FftField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    IOPattern<DigestBridge<D>>:
        StarkIOPattern<D, F> + FriIOPattern<D, F> + FieldIOPattern<F> + DigestIOWritter<D>,
{
    #[allow(dead_code)]
    security_bits: usize,
    steps: usize,
    blowup_factor: usize,
    rounds: usize,
    constrain_queries: usize,
    degree: usize,
    fri_config: FriConfig<D, F>,
    merkle_config: MerkleTreeConfig<D, F>,
    io: IOPattern<DigestBridge<D>>,
}

impl<D, F> StarkConfig<D, F>
where
    F: FftField + PrimeField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    IOPattern<DigestBridge<D>>:
        StarkIOPattern<D, F> + FriIOPattern<D, F> + FieldIOPattern<F> + DigestIOWritter<D>,
{
    pub fn new(
        security_bits: usize,
        blowup_factor: usize,
        steps: usize,
        trace_columns: usize,
    ) -> Self {
        let (constrain_queries, fri_queries) =
            Self::num_queries_from_config(security_bits, blowup_factor, steps);
        let degree = steps - 1;
        let rounds = ceil_log2_k((steps * blowup_factor) + 1, 2);

        Self {
            security_bits,
            steps,
            blowup_factor,
            rounds,
            constrain_queries,
            degree,
            fri_config: FriConfig {
                queries: fri_queries,
                blowup_factor,
                rounds,
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
            io: <IOPattern<DigestBridge<D>> as StarkIOPattern<D, F>>::new_stark(
                rounds,
                constrain_queries,
                fri_queries,
                "🐺",
            ),
        }
    }

    fn num_queries_from_config(
        security_bits: usize,
        blowup_factor: usize,
        steps: usize,
    ) -> (usize, usize) {
        if security_bits < 20 {
            error!("STARK Config: security bits has to be at least 20");
            panic!("");
        }
        let log_steps = ceil_log2_k(steps, 2);
        let linking_queries = security_bits.div_ceil(F::MODULUS_BIT_SIZE as usize - log_steps);

        let rounds = ceil_log2_k(steps * blowup_factor, 2);
        let rho = 1f64 / blowup_factor as f64;
        let denominator = (2f64 / (1f64 + rho)).log2();
        let total_fri_queries = (security_bits as f64) / denominator;
        let round_fri_queries = (total_fri_queries / (rounds as f64)).ceil() as usize;

        (linking_queries, round_fri_queries)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::Goldilocks;
    use sha2::Sha256;

    #[test]
    #[should_panic]
    fn test_stark_config_with_low_security_bits() {
        // test query_number is at least 1
        StarkConfig::<Sha256, Goldilocks>::num_queries_from_config(1, 4, 128);
    }

    #[test]
    fn test_stark_config_query_numbers() {
        let blowup_factor = 4;
        let steps = 129;
        assert_eq!(Goldilocks::MODULUS_BIT_SIZE, 64);

        // test query_number is at least 1
        let (constrain_queries, fri_queries) =
            StarkConfig::<Sha256, Goldilocks>::num_queries_from_config(20, blowup_factor, steps);
        assert_eq!(constrain_queries, 1);
        assert_eq!(fri_queries, 3);

        let (constrain_queries, fri_queries) =
            StarkConfig::<Sha256, Goldilocks>::num_queries_from_config(20, 2, 9);
        assert_eq!(constrain_queries, 1);
        assert_eq!(fri_queries, 10);

        // test that queries grow linearly with security bits
        let (constrain_queries, fri_queries) =
            StarkConfig::<Sha256, Goldilocks>::num_queries_from_config(128, blowup_factor, steps);
        assert_eq!(constrain_queries, 3);
        assert_eq!(fri_queries, 19);
        let (constrain_queries, fri_queries) =
            StarkConfig::<Sha256, Goldilocks>::num_queries_from_config(256, blowup_factor, 513);
        assert_eq!(constrain_queries, 5);
        assert_eq!(fri_queries, 32);
    }
}
