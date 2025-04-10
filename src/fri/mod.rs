pub mod fiatshamir;
pub mod prover;
pub mod verifier;

use crate::merkle::{MerklePath, MerkleTreeConfig};
use ark_ff::PrimeField;
use digest::Digest;

#[derive(Clone)]
pub struct FriConfig<D: Digest, F: PrimeField> {
    merkle_config: MerkleTreeConfig<D, F>,
    blowup_factor: usize,
}

#[derive(Clone)]
pub struct FriProof<D: Digest, F: PrimeField> {
    // transcript: Vec<u8>,
    points: Vec<[(F, F); 3]>,
    queries: Vec<[MerklePath<D, F>; 2]>,
    quotients: Vec<Vec<F>>,
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::Goldilocks;
    use crate::fri::fiatshamir::FriIOPattern;
    use crate::merkle::{MerkleRoot, MerkleTreeConfig};
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use nimue::{DigestBridge, IOPattern};
    use sha2::Sha256;
    use std::marker::PhantomData;

    use super::prover::FriProver;
    use super::verifier::FriVerifier;

    #[test]
    fn test_fri_new() {
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let degree = poly.degree();
        let io: IOPattern<DigestBridge<Sha256>> = FriIOPattern::<_, Goldilocks>::new_fri("🍟", 3);
        let mut transcript = io.to_merlin();

        let merkle_config = MerkleTreeConfig {
            leafs_per_node: 2,
            inner_children: 2,
            _digest: PhantomData::<Sha256>,
            _field: PhantomData::<Goldilocks>,
        };

        let config = FriConfig {
            merkle_config,
            blowup_factor: 2,
        };

        let fri_prover = FriProver::<Sha256, _>::new(&mut transcript, poly, config.clone());

        let commit = fri_prover.get_initial_commit();

        let (proof, transcript) = fri_prover.prove();
        let verifier = FriVerifier::<Sha256, Goldilocks>::new(MerkleRoot(commit), degree, config);

        let mut arthur = io.to_arthur(&transcript);
        assert!(verifier.verify(proof, &mut arthur));
    }
}
