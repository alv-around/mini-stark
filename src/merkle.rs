use ark_ff::PrimeField;
use digest::{generic_array::GenericArray, Digest, OutputSizeUser};
use std::ops::Range;

pub trait Tree<const N: usize> {
    type Input;
    type Inner;

    fn new(inputs: &[Self::Input]) -> Self;
    fn root(&self) -> Self::Inner;
    // fn get_parent(&self, node: &Self::Inner) -> &Self::Inner;
    // fn get_children(&self, node: &Self::Inner) -> &[Self::Inner; N];
    fn calculate_from_external_nodes(children: &[Self::Input]) -> Self::Inner;
    fn calculate_from_inner_nodes(children: &[Self::Inner]) -> Self::Inner;
}

pub type Hash<D> = GenericArray<u8, <D as OutputSizeUser>::OutputSize>;

pub struct MerkleTree<const N: usize, D: Digest, F: PrimeField> {
    external_leafs: Vec<F>,
    inner_nodes: Vec<Hash<D>>,
    levels: usize,
}

impl<const N: usize, D: Digest, F: PrimeField> Tree<N> for MerkleTree<N, D, F> {
    type Input = F;
    type Inner = Hash<D>;

    fn new(inputs: &[F]) -> Self {
        let leaf_num = inputs.len();
        let log2_inputs = leaf_num.trailing_zeros() as usize;
        // TODO: refactor
        let levels = match N {
            2 => log2_inputs,
            4 => log2_inputs >> 1,
            8 => log2_inputs >> 2,
            16 => log2_inputs >> 3,
            _ => panic!("MerkleTree not implemented for width: {}", N),
        };
        // assert tree is full
        assert_eq!(
            levels * N,
            leaf_num,
            "Tree is not full! input length must be a power of {N}"
        );

        let numerator = 1 - (N as i64).pow(levels as u32);
        let denominator = 1 - N as i64;
        let node_num = (numerator / denominator) as usize;
        let mut inner_nodes = Vec::with_capacity(node_num);

        for external in inputs.chunks(N) {
            let parent = MerkleTree::<N, D, _>::calculate_from_external_nodes(external);
            inner_nodes.push(parent);
        }

        let mut distance = node_num;
        for _ in leaf_num..node_num {
            let children = &inner_nodes[distance..distance + N];
            let parent = MerkleTree::<N, D, F>::calculate_from_inner_nodes(children);
            inner_nodes.push(parent);

            distance -= N - 1;
        }

        Self {
            external_leafs: inputs.to_vec(),
            inner_nodes,
            levels,
        }
    }

    fn root(&self) -> Hash<D> {
        let root = self.inner_nodes.last().unwrap().clone();
        root
    }

    fn calculate_from_external_nodes(children: &[F]) -> Hash<D> {
        let mut hasher = D::new();
        for child in children {
            hasher.update(child.to_string());
        }
        hasher.finalize()
    }

    fn calculate_from_inner_nodes(children: &[Hash<D>]) -> Hash<D> {
        let mut hasher = D::new();
        for child in children {
            hasher.update(child)
        }
        hasher.finalize()
    }
}

impl<const N: usize, F: PrimeField, D: Digest> MerkleTree<N, D, F> {
    fn get_neighbor_idx(&self, index: usize) -> Range<usize> {
        let remainder = index % N;
        let start_idx = index - remainder;
        let end_idx = index + N;
        start_idx..end_idx
    }

    fn get_parent_idx(&self, index: usize) -> usize {
        let leafs_num = self.external_leafs.len();
        if index < leafs_num {
            return index / N;
        } else {
            let n = self.inner_nodes.len();
            let inner_index = index - leafs_num;
            return (n - 1) - ((n - inner_index - 2) / N);
        }
    }

    fn get_leaf_index(&self, node: &F) -> Option<usize> {
        for (i, value) in self.external_leafs.iter().enumerate() {
            if *node == *value {
                return Some(i);
            }
        }
        None
    }

    fn get_leaf_neighbours(&self, index: usize) -> Vec<F> {
        let neighbors_idx = self.get_neighbor_idx(index);
        self.external_leafs[neighbors_idx].to_vec()
    }

    fn get_inner_neighbours(&self, index: usize) -> Vec<Hash<D>> {
        let neighbors_idx = self.get_neighbor_idx(index);
        self.inner_nodes[neighbors_idx].to_vec()
    }

    // TODO: add const LogN to code
    fn calculate_path(&self, index: usize) -> Vec<Vec<Hash<D>>> {
        let mut path = Vec::new();
        let mut current_idx = index;
        for _ in 1..self.levels {
            let neighbor = self.get_inner_neighbours(current_idx);
            path.push(neighbor);

            let parent = self.get_parent_idx(current_idx);
            current_idx = parent;
        }

        path
    }

    pub fn generate_proof(&self, leaf: &F) -> Result<MerklePath<D, F>, &str> {
        let leaf_index = self.get_leaf_index(&leaf);
        if leaf_index.is_none() {
            return Err("leaf is not included in the tree");
        }
        let leaf_index = leaf_index.unwrap();

        let leaf_neighbours = self.get_leaf_neighbours(leaf_index);
        let leaf_parent = self.get_parent_idx(leaf_index);
        let path = self.calculate_path(leaf_parent);
        Ok(MerklePath {
            leaf_neighbours,
            path,
        })
    }
}

pub struct MerklePath<D: Digest, F: PrimeField> {
    leaf_neighbours: Vec<F>,
    path: Vec<Vec<Hash<D>>>,
}

pub struct MerkleRoot<D: Digest>(Hash<D>);

impl<D: Digest> MerkleRoot<D> {
    pub fn check_proof<const N: usize, F: PrimeField>(
        &self,
        leaf: &F,
        proof: MerklePath<D, F>,
    ) -> bool {
        if !proof.leaf_neighbours.contains(leaf) {
            return false;
        };

        let mut previous =
            MerkleTree::<N, D, F>::calculate_from_external_nodes(&proof.leaf_neighbours);

        for level in proof.path {
            if !level.contains(&previous) {
                return false;
            }

            previous = MerkleTree::<N, D, F>::calculate_from_inner_nodes(&level);
        }

        if previous == self.0 {
            return true;
        }
        false
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::Goldilocks;
    use sha2::Sha256;

    use std::panic;

    const TWO: usize = 2;
    const FOUR: usize = 4;
    const EIGHT: usize = 8;
    const SIXTEEN: usize = 16;
    const WIDTHS: [usize; 4] = [TWO, FOUR, EIGHT, SIXTEEN];

    fn make_tree<const N: usize>() -> MerkleTree<N, Sha256, Goldilocks> {
        let leafs: Vec<Goldilocks> = (0..16).map(|i| Goldilocks::from(i)).collect();
        MerkleTree::<N, Sha256, _>::new(&leafs)
    }

    #[test]
    fn test_panic_for_not_full_trees() {
        let leafs = vec![
            Goldilocks::from(0),
            Goldilocks::from(1),
            Goldilocks::from(2),
        ];

        // TODO: write macro to avoid code repetition
        let result = panic::catch_unwind(|| {
            MerkleTree::<TWO, Sha256, _>::new(&leafs);
        });
        assert!(result.is_err(), "Tree of width: {} to panic", TWO);

        let result = panic::catch_unwind(|| {
            MerkleTree::<FOUR, Sha256, _>::new(&leafs);
        });
        assert!(result.is_err(), "Tree of width: {} to panic", FOUR);

        let result = panic::catch_unwind(|| {
            MerkleTree::<EIGHT, Sha256, _>::new(&leafs);
        });
        assert!(result.is_err(), "Tree of width: {} to panic", EIGHT);

        let result = panic::catch_unwind(|| {
            MerkleTree::<SIXTEEN, Sha256, _>::new(&leafs);
        });
        assert!(result.is_err(), "Tree of width: {} to panic", SIXTEEN);
    }

    // #[test]
    // fn test_valid_merkle_tree() {
    //     let tree = make_tree();
    //     let bytes: [u8; 32] = [
    //         59, 130, 140, 79, 75, 72, 197, 212, 203, 85, 98, 164, 116, 236, 158, 47, 216, 213, 84,
    //         111, 174, 64, 233, 7, 50, 239, 99, 88, 146, 228, 39, 32,
    //     ];
    //
    //     let hash_value: Hash<Sha256> = Hash::<Sha256>::clone_from_slice(&bytes);
    //
    //     assert_eq!(tree.leaf_number, 8);
    //     assert_eq!(tree.nodes.len(), 15);
    //     assert_eq!(tree.get_root().0, hash_value);
    //
    //     let hash_idx2 = MerkleTree::<Sha256, _>::hash_leaf(Goldilocks::from(2));
    //     let hash_idx5 = MerkleTree::<Sha256, _>::hash_leaf(Goldilocks::from(5));
    //     assert_eq!(tree.get_leaf_index(&hash_idx2), Some(2));
    //     assert_eq!(tree.get_leaf_index(&hash_idx5), Some(5));
    // }
    //
    // #[test]
    // fn test_neighbor_index() {
    //     let tree = make_tree();
    //
    //     assert_eq!(tree.get_neighbor_idx(4).0, 5);
    //     assert_eq!(tree.get_neighbor_idx(7).0, 6);
    // }
    //
    // #[test]
    // fn test_merkle_tree_parent_index() {
    //     let tree = make_tree();
    //
    //     // assert parents
    //     assert_eq!(tree.get_parent_idx(0), 8);
    //     assert_eq!(tree.get_parent_idx(3), 9);
    //     assert_eq!(tree.get_parent_idx(10), 13);
    //     assert_eq!(tree.get_parent_idx(12), 14);
    // }
    //
    // #[test]
    // fn test_merkle_path() {
    //     let tree = make_tree();
    //
    //     let (path1, _) = tree.calculate_path(1);
    //     let path_raw = vec![
    //         *Hash::<Sha256>::from_slice(&[
    //             95, 236, 235, 102, 255, 200, 111, 56, 217, 82, 120, 108, 109, 105, 108, 121, 194,
    //             219, 194, 57, 221, 78, 145, 180, 103, 41, 215, 58, 39, 251, 87, 233,
    //         ]),
    //         *Hash::<Sha256>::from_slice(&[
    //             169, 245, 179, 171, 97, 226, 131, 87, 207, 205, 20, 226, 180, 35, 151, 248, 150,
    //             174, 234, 141, 105, 152, 209, 158, 109, 168, 85, 132, 225, 80, 210, 180,
    //         ]),
    //         *Hash::<Sha256>::from_slice(&[
    //             3, 2, 201, 111, 69, 171, 190, 173, 178, 56, 120, 51, 26, 155, 164, 6, 7, 139, 208,
    //             189, 93, 194, 2, 193, 2, 175, 123, 153, 134, 36, 159, 1,
    //         ]),
    //     ];
    //     assert_eq!(path1, path_raw);
    //
    //     let (path2, _) = tree.calculate_path(4);
    //     let path_raw = vec![
    //         *Hash::<Sha256>::from_slice(&[
    //             239, 45, 18, 125, 227, 123, 148, 43, 170, 208, 97, 69, 229, 75, 12, 97, 154, 31,
    //             34, 50, 123, 46, 187, 207, 190, 199, 143, 85, 100, 175, 227, 157,
    //         ]),
    //         *Hash::<Sha256>::from_slice(&[
    //             19, 72, 67, 175, 127, 200, 242, 153, 80, 177, 225, 223, 183, 196, 151, 82, 224,
    //             247, 183, 17, 180, 88, 238, 154, 227, 197, 202, 34, 1, 102, 214, 136,
    //         ]),
    //         *Hash::<Sha256>::from_slice(&[
    //             196, 120, 254, 173, 12, 137, 183, 149, 64, 99, 143, 132, 76, 136, 25, 217, 164, 40,
    //             23, 99, 175, 146, 114, 199, 243, 150, 135, 118, 182, 5, 35, 69,
    //         ]),
    //     ];
    //     assert_eq!(path2, path_raw);
    // }
    //
    // #[test]
    // fn test_generate_check_proof() {
    //     let tree = make_tree();
    //     let proof = tree.generate_proof(Goldilocks::from(4)).unwrap();
    //     let path_raw = vec![
    //         *Hash::<Sha256>::from_slice(&[
    //             239, 45, 18, 125, 227, 123, 148, 43, 170, 208, 97, 69, 229, 75, 12, 97, 154, 31,
    //             34, 50, 123, 46, 187, 207, 190, 199, 143, 85, 100, 175, 227, 157,
    //         ]),
    //         *Hash::<Sha256>::from_slice(&[
    //             19, 72, 67, 175, 127, 200, 242, 153, 80, 177, 225, 223, 183, 196, 151, 82, 224,
    //             247, 183, 17, 180, 88, 238, 154, 227, 197, 202, 34, 1, 102, 214, 136,
    //         ]),
    //         *Hash::<Sha256>::from_slice(&[
    //             196, 120, 254, 173, 12, 137, 183, 149, 64, 99, 143, 132, 76, 136, 25, 217, 164, 40,
    //             23, 99, 175, 146, 114, 199, 243, 150, 135, 118, 182, 5, 35, 69,
    //         ]),
    //     ];
    //     assert_eq!(proof.path, path_raw);
    //
    //     let root = tree.get_root();
    //     let check = root.check_proof(proof);
    // assert!(check);
    // }
}
