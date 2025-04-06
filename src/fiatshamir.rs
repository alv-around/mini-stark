use crate::fri::fiatshamir::FriIOPattern;
use ark_ff::Field;
use digest::core_api::BlockSizeUser;
use digest::{Digest, FixedOutputReset, OutputSizeUser};
use nimue::{
    plugins::ark::FieldIOPattern, Arthur, ByteIOPattern, DigestBridge, DuplexHash, IOPattern,
    ProofResult,
};
use std::mem::size_of;

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
    fn new_stark(domain_size_log: usize, num_queries: usize, domsep: &str) -> Self;
}

impl<D, F> StarkIOPattern<D, F> for IOPattern<DigestBridge<D>>
where
    F: Field,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    IOPattern<DigestBridge<D>>: FieldIOPattern<F> + DigestIOWritter<D> + FriIOPattern<D, F>,
{
    fn new_stark(domain_size_log: usize, num_queries: usize, domsep: &str) -> Self {
        IOPattern::<DigestBridge<D>>::new(domsep)
            .add_digest(2, "commit to trace & quotients")
            .challenge_scalars(1, "batching: retrieve random scalar r")
            .challenge_bytes(8 * num_queries, "retrive random queries")
            .add_fri(domain_size_log)
    }
}

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
