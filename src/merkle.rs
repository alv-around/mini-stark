use crate::util::logarithm_of_two_k;
use ark_ff::PrimeField;
use digest::{generic_array::GenericArray, Digest, OutputSizeUser};
use std::ops::Range;

pub trait Tree<const N: usize> {
    type Input;
    type Inner;

    fn new(inputs: &[Self::Input]) -> Self;
    fn root(&self) -> Self::Inner;
    fn get_node_number(&self) -> usize;
    // fn get_parent(&self, node: &Self::Inner) -> &Self::Inner;
    // fn get_children(&self, node: &Self::Inner) -> &[Self::Inner; N];
    fn calculate_from_leafs(children: &[Self::Input]) -> Self::Inner;
    fn calculate_from_nodes(children: &[Self::Inner]) -> Self::Inner;
}

pub type Hash<D> = GenericArray<u8, <D as OutputSizeUser>::OutputSize>;

pub struct MerkleTree<const N: usize, D: Digest, F: PrimeField> {
    leafs: Vec<F>,
    nodes: Vec<Hash<D>>,
    levels: usize,
}

impl<const N: usize, D: Digest, F: PrimeField> Tree<N> for MerkleTree<N, D, F> {
    type Input = F;
    type Inner = Hash<D>;

    fn new(inputs: &[F]) -> Self {
        let leaf_num = inputs.len();
        let levels = match logarithm_of_two_k::<N>(leaf_num) {
            Ok(log) => log,
            Err(error_str) => panic!("{}", error_str),
        };

        // assert tree is full
        assert_eq!(
            N.pow(levels as u32),
            leaf_num,
            "Tree is not full! input length must be a power of {N}"
        );

        // number of nodes
        let numerator = 1 - (N as i64).pow(levels as u32);
        let denominator = 1 - N as i64;
        let node_num = (numerator / denominator) as usize;
        let mut nodes = Vec::with_capacity(node_num);
        println!("Number of nodes in the tree: {node_num}");

        let mut distance = leaf_num;
        for external in inputs.chunks(N) {
            let parent = MerkleTree::<N, D, _>::calculate_from_leafs(external);
            nodes.push(parent);
            distance -= N - 1;
        }

        let mut current_nodes = leaf_num / N;
        let nodes_left = node_num - current_nodes;
        for _ in 0..nodes_left {
            let idx = current_nodes - distance;
            let children = &nodes[idx..idx + N];
            let parent = MerkleTree::<N, D, F>::calculate_from_nodes(children);
            nodes.push(parent);
            distance -= N - 1;
            current_nodes += 1;
        }

        Self {
            leafs: inputs.to_vec(),
            nodes,
            levels,
        }
    }

    fn root(&self) -> Hash<D> {
        let root = self.nodes.last().unwrap().clone();
        root
    }

    fn get_node_number(&self) -> usize {
        self.leafs.len() + self.nodes.len()
    }

    fn calculate_from_leafs(children: &[F]) -> Hash<D> {
        let mut hasher = D::new();
        for child in children {
            hasher.update(child.to_string());
        }
        hasher.finalize()
    }

    fn calculate_from_nodes(children: &[Hash<D>]) -> Hash<D> {
        let mut hasher = D::new();
        for child in children {
            hasher.update(child)
        }
        hasher.finalize()
    }
}

impl<const N: usize, F: PrimeField, D: Digest> MerkleTree<N, D, F> {
    // TODO: better error handling
    fn get_neighbor_idx(&self, index: usize) -> Range<usize> {
        if index >= self.get_node_number() {
            panic!("index outside of tree length");
        }

        let remainder = index % N;
        let start_idx = index - remainder;
        let end_idx = start_idx + N;
        start_idx..end_idx
    }

    // TODO: better error handling
    fn get_parent_idx(&self, index: usize) -> usize {
        let leafs_num = self.leafs.len();
        if index < leafs_num {
            index / N
        } else if index < self.get_node_number() {
            let n = self.nodes.len();
            let inner_index = index - leafs_num;
            let parent = (n - 1) - ((n - inner_index - 2) / N);
            return parent;
        } else {
            panic!("index outside of tree length");
        }
    }

    fn get_leaf_index(&self, node: &F) -> Option<usize> {
        for (i, value) in self.leafs.iter().enumerate() {
            if *node == *value {
                return Some(i);
            }
        }
        None
    }

    fn get_leaf_neighbours(&self, index: usize) -> Vec<F> {
        let neighbors_idx = self.get_neighbor_idx(index);
        self.leafs[neighbors_idx].to_vec()
    }

    fn get_inner_neighbours(&self, index: usize) -> Vec<Hash<D>> {
        let neighbors_idx = self.get_neighbor_idx(index);
        self.nodes[neighbors_idx].to_vec()
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

        let mut previous = MerkleTree::<N, D, F>::calculate_from_leafs(&proof.leaf_neighbours);

        for level in proof.path {
            if !level.contains(&previous) {
                return false;
            }

            previous = MerkleTree::<N, D, F>::calculate_from_nodes(&level);
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

    fn make_tree<const N: usize>() -> MerkleTree<N, Sha256, Goldilocks> {
        let leafs: Vec<Goldilocks> = (0..16).map(Goldilocks::from).collect();
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

    #[test]
    fn test_node_calculation() {
        let tree = make_tree::<TWO>();
        assert_eq!(tree.get_node_number(), 31);
        assert_eq!(tree.leafs.len(), 16);
        assert_eq!(tree.nodes.len(), 15);

        let tree = make_tree::<FOUR>();
        assert_eq!(tree.get_node_number(), 21);
        assert_eq!(tree.leafs.len(), 16);
        assert_eq!(tree.nodes.len(), 5);

        let tree = make_tree::<SIXTEEN>();
        assert_eq!(tree.get_node_number(), 17);
        assert_eq!(tree.leafs.len(), 16);
        assert_eq!(tree.nodes.len(), 1);
    }

    #[test]
    fn test_neighbor_index() {
        let tree = make_tree::<TWO>();
        assert_eq!(tree.get_neighbor_idx(4), 4..6);
        assert_eq!(tree.get_neighbor_idx(7), 6..8);

        let tree = make_tree::<FOUR>();
        assert_eq!(tree.get_neighbor_idx(4), 4..8);
        assert_eq!(tree.get_neighbor_idx(7), 4..8);

        let tree = make_tree::<SIXTEEN>();
        assert_eq!(tree.get_neighbor_idx(4), 0..16);
        assert_eq!(tree.get_neighbor_idx(7), 0..16);
    }

    #[test]
    fn test_merkle_tree_parent_index() {
        let tree = make_tree::<TWO>();
        assert_eq!(tree.get_parent_idx(9), 4);
        assert_eq!(tree.get_parent_idx(17), 8);
        assert_eq!(tree.get_parent_idx(26), 13);
        assert_eq!(tree.get_parent_idx(28), 14);

        let tree = make_tree::<FOUR>();
        assert_eq!(tree.get_parent_idx(2), 0);
        assert_eq!(tree.get_parent_idx(9), 2);
        assert_eq!(tree.get_parent_idx(16), 4);
        assert_eq!(tree.get_parent_idx(19), 4);

        let tree = make_tree::<SIXTEEN>();
        assert_eq!(tree.get_parent_idx(0), 0);
        assert_eq!(tree.get_parent_idx(9), 0);

        // test that calling and index out of tree length panics
        let result = panic::catch_unwind(|| {
            tree.get_parent_idx(tree.get_node_number());
        });
        assert!(result.is_err());
    }

    // TODO: modify python script to reproduce hashing outcomes
    // #[test]
    // fn test_hashing_replicability()
    // let bytes: [u8; 32] = [
    //     59, 130, 140, 79, 75, 72, 197, 212, 203, 85, 98, 164, 116, 236, 158, 47, 216, 213, 84,
    //     111, 174, 64, 233, 7, 50, 239, 99, 88, 146, 228, 39, 32,
    // ];
    //
    // let hash_value: Hash<Sha256> = Hash::<Sha256>::clone_from_slice(&bytes);
    // assert_eq!(tree.root(), hash_value);
    //
    // let hash_idx2 = MerkleTree::<TWO, Sha256, _>::hash_leaf(Goldilocks::from(2));
    // let hash_idx5 = MerkleTree::<TWO, Sha256, _>::hash_leaf(Goldilocks::from(5));
    // assert_eq!(tree.get_leaf_index(&hash_idx2), Some(2));
    // assert_eq!(tree.get_leaf_index(&hash_idx5), Some(5));
    // }
    // #[test]
    // fn test_merkle_path() {
    //     let tree = make_tree::<TWO>();
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
    //     let tree = make_tree::<TWO>();
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
