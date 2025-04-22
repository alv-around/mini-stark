use crate::Hash;
use ark_ff::Field;
use digest::core_api::BlockSizeUser;
use digest::{generic_array::GenericArray, Digest, FixedOutputReset, OutputSizeUser};
use nimue::{
    plugins::ark::FieldIOPattern, Arthur, ByteIOPattern, DigestBridge, DuplexHash, IOPattern,
    ProofResult,
};
use std::mem::size_of;

#[allow(dead_code)]
pub trait UsizeIOWritter: ByteIOPattern + Sized {
    fn add_usize(self, count: usize, label: &str) -> Self {
        self.add_bytes(count * 8, label)
    }
}

pub trait DigestIOWritter<D: Digest> {
    fn add_digest(self, count: usize, label: &str) -> Self;
}

impl<D> DigestIOWritter<D> for IOPattern<DigestBridge<D>>
where
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    IOPattern<DigestBridge<D>>: ByteIOPattern,
{
    fn add_digest(self, count: usize, label: &str) -> Self {
        self.add_bytes(count * <D as OutputSizeUser>::output_size(), label)
    }
}

pub trait StarkIOPattern<D: Digest, F: Field> {
    fn new_stark(
        domain_size_log: usize,
        constrain_queries: usize,
        fri_queries: usize,
        domsep: &str,
    ) -> Self;
}

impl<D, F> StarkIOPattern<D, F> for IOPattern<DigestBridge<D>>
where
    F: Field,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    Self: FieldIOPattern<F> + DigestIOWritter<D> + FriIOPattern<D, F>,
{
    fn new_stark(
        rounds: usize,
        constrain_queries: usize,
        fri_queries: usize,
        domsep: &str,
    ) -> Self {
        IOPattern::new(domsep)
            .add_digest(1, "commit to original trace")
            .challenge_scalars(1, "ZK: pick random shift of domain")
            .add_digest(1, "commit to quotients")
            .challenge_scalars(1, "batching: retrieve random scalar r")
            .challenge_bytes(8 * constrain_queries, "retrive random queries")
            .add_fri(rounds, fri_queries)
    }
}

#[allow(dead_code)]
pub trait UsizeReader {
    fn next_usize(&mut self, length: usize) -> ProofResult<Vec<usize>>;
}

impl<H: DuplexHash> UsizeReader for Arthur<'_, H, u8> {
    fn next_usize(&mut self, length: usize) -> ProofResult<Vec<usize>> {
        let usize_size = size_of::<usize>();
        let mut byte_challenges = vec![0u8; length * usize_size];
        self.fill_next_units(&mut byte_challenges)?;
        let challenges: Vec<usize> = byte_challenges
            .chunks_exact(usize_size)
            .map(|x| usize::from_le_bytes(x.try_into().unwrap()))
            .collect();
        Ok(challenges)
    }
}

pub trait FriIOPattern<D: Digest, F: Field> {
    fn new_fri(domsep: &str, rounds: usize, queries: usize) -> Self;
    fn add_fri(self, rounds: usize, queries: usize) -> Self;
}

impl<D, F> FriIOPattern<D, F> for IOPattern<DigestBridge<D>>
where
    F: Field,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    IOPattern<DigestBridge<D>>: FieldIOPattern<F> + DigestIOWritter<D>,
{
    fn new_fri(domsep: &str, rounds: usize, queries: usize) -> Self {
        IOPattern::new(domsep).add_fri(rounds, queries)
    }

    fn add_fri(self, rounds: usize, queries: usize) -> Self {
        let mut this = self;
        for _ in 0..rounds - 1 {
            this = this
                .challenge_scalars(1, "random scalar challenge: polynomial folding")
                .add_digest(1, "add merkle commit: commit to fri round");
        }

        this = this
            // .add_digest(1, "add merkle commit: commit to last fri round")
            .challenge_bytes(
                8 * queries,
                "query phase: choose a random element in the domain",
            );

        this
    }
}

pub trait DigestReader<D: Digest> {
    fn next_digest(&mut self) -> ProofResult<Hash<D>>;
}

impl<D> DigestReader<D> for Arthur<'_, DigestBridge<D>, u8>
where
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
{
    fn next_digest(&mut self) -> ProofResult<Hash<D>> {
        let mut digest_bytes = vec![0u8; <D as OutputSizeUser>::output_size()];
        self.fill_next_units(&mut digest_bytes)?;
        Ok(GenericArray::from_exact_iter(digest_bytes).unwrap())
    }
}
