pub mod fiatshamir;
pub mod prover;
pub mod verifier;

use crate::merkle::{Hash, MerklePath};
use ark_ff::PrimeField;
use digest::Digest;

pub struct FriProof<'a, D: Digest, F: PrimeField> {
    transcript: &'a [u8],
    commits: Vec<Hash<D>>,
    points: Vec<[(F, F); 3]>,
    queries: Vec<[MerklePath<D, F>; 3]>,
    quotients: Vec<Vec<F>>,
}

#[cfg(test)]
mod test {
    use crate::field::Goldilocks;
    use crate::fri::fiatshamir::FriIOPattern;
    use crate::merkle::{Hash, MerkleRoot};
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
        let hash = *Hash::<Sha256>::from_slice(&[
            196, 120, 254, 173, 12, 137, 183, 149, 64, 99, 143, 132, 76, 136, 25, 217, 164, 40, 23,
            99, 175, 146, 114, 199, 243, 150, 135, 118, 182, 5, 35, 69,
        ]);

        let mut fri_prover =
            FriProver::<TWO, Sha256, _>::new(transcript.to_merlin(), poly, blowup_factor);
        let proof = fri_prover.prove();
        let verifier = FriVerifier::<TWO, Sha256, Goldilocks>::new(
            transcript,
            MerkleRoot(hash),
            degree,
            blowup_factor,
        );
        assert!(verifier.verify(proof));
    }
}
