use ark_ff::Field;
use digest::core_api::BlockSizeUser;
use digest::{Digest, FixedOutputReset};
use nimue::{
    plugins::ark::{FieldIOPattern, IOPattern, ProofResult},
    ByteIOPattern, DigestBridge,
};

pub trait FriIOPattern<D: Digest, F: Field> {
    fn new_fri(domsep: &str, round_numbers: usize) -> Self;
    fn add_fri(self, round_numbers: usize) -> Self;
}

pub trait DigestWriter<D: Digest> {
    fn add_digest(&mut self, digest: D) -> ProofResult<()>;
}

pub trait DigestReader<D: Digest> {
    fn read_digest(&mut self) -> ProofResult<D>;
}

impl<D, F> FriIOPattern<D, F> for IOPattern<DigestBridge<D>>
where
    F: Field,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    IOPattern<DigestBridge<D>>: FieldIOPattern<F>,
{
    fn new_fri(domsep: &str, round_numbers: usize) -> Self {
        IOPattern::new(domsep).add_fri(round_numbers)
    }

    fn add_fri(mut self, round_numbers: usize) -> Self {
        for _ in 0..round_numbers - 1 {
            self = self
                .add_bytes(
                    <D as digest::OutputSizeUser>::output_size(),
                    "add merkle commit: commit to fri round",
                )
                .challenge_scalars(1, "random scalar challenge: polynomial folding");
        }

        self = self
            .add_bytes(
                <D as digest::OutputSizeUser>::output_size(),
                "add merkle commit: commit to last fri round",
            )
            .challenge_bytes(4, "random scallar challenge: query_phase");

        self
    }
}
