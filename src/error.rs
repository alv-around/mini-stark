use error_set::error_set;
use nimue::{IOPatternError, ProofError};

error_set! {
    ProverError = {
        TranscriptError1(IOPatternError),
        TranscriptError2(ProofError),
    } || MerkleProofError;
    VerifierError = {
        TranscriptError1(ProofError),
        TranscriptError2(IOPatternError),
    };
    MerkleProofError = {
        #[display("Error generating Merkle proof: {msg}")]
        LeafNotFound{
            msg: &'static str
        },
        #[display("Error generating Merkle proof: {msg}")]
        OutOfRangeError{
            msg: &'static str},
    };
}
