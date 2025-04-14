use crate::fiatshamir::DigestIOWritter;
use crate::Hash;
use ark_ff::Field;
use digest::core_api::BlockSizeUser;
use digest::{generic_array::GenericArray, Digest, FixedOutputReset, OutputSizeUser};
pub use nimue::Arthur;
use nimue::ProofResult;
use nimue::{
    plugins::ark::{FieldIOPattern, IOPattern},
    ByteIOPattern, DigestBridge,
};

pub trait FriIOPattern<D: Digest, F: Field> {
    fn new_fri(domsep: &str, round_numbers: usize, queries: usize) -> Self;
    fn add_fri(self, round_numbers: usize, queries: usize) -> Self;
}

impl<D, F> FriIOPattern<D, F> for IOPattern<DigestBridge<D>>
where
    F: Field,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    IOPattern<DigestBridge<D>>: FieldIOPattern<F> + DigestIOWritter<D>,
{
    fn new_fri(domsep: &str, round_numbers: usize, queries: usize) -> Self {
        IOPattern::new(domsep).add_fri(round_numbers, queries)
    }

    fn add_fri(mut self, round_numbers: usize, queries: usize) -> Self {
        for _ in 0..round_numbers - 1 {
            self = self
                .add_digest(1, "add merkle commit: commit to fri round")
                .challenge_scalars(1, "random scalar challenge: polynomial folding");
        }

        self = self
            .add_digest(1, "add merkle commit: commit to last fri round")
            .challenge_bytes(
                8 * queries,
                "query phase: choose a random element in the domain",
            );

        self
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
