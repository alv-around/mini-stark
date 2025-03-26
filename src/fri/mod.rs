pub mod fiatshamir;
pub mod prover;
pub mod verifier;

use crate::merkle::MerklePath;
use ark_ff::PrimeField;
use digest::Digest;

pub struct FriProof<'a, D: Digest, F: PrimeField> {
    transcript: &'a [u8],
    points: Vec<[(F, F); 3]>,
    queries: Vec<[MerklePath<D, F>; 3]>,
    quotients: Vec<Vec<F>>,
}

#[cfg(test)]
mod test {
    use crate::field::Goldilocks;
    use crate::fri::fiatshamir::FriIOPattern;
    use crate::merkle::MerkleRoot;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use nimue::{DigestBridge, IOPattern};
    use sha2::Sha256;

    use super::prover::FriProver;
    use super::verifier::FriVerifier;

    const TWO: usize = 2;

    #[test]
    fn test_fri_new() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let degree = poly.degree();
        let transcript: IOPattern<DigestBridge<Sha256>> =
            FriIOPattern::<_, Goldilocks>::new_fri("🍟", 3);
        let mut fri_prover =
            FriProver::<TWO, Sha256, _>::new(transcript.to_merlin(), poly, blowup_factor);

        let commit = fri_prover.get_initial_commit();

        let proof = fri_prover.prove();
        let verifier = FriVerifier::<TWO, Sha256, Goldilocks>::new(
            transcript,
            MerkleRoot(commit),
            degree,
            blowup_factor,
        );
        assert!(verifier.verify(proof));
    }
}
