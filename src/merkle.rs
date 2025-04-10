use crate::{util::logarithm_of_two_k, Hash};
use ark_ff::PrimeField;
use digest::Digest;
use std::marker::PhantomData;
use std::ops::Range;

pub trait Tree {
    type Input;
    type Inner;
    type Config;

    fn new(inputs: &[Self::Input], config: Self::Config) -> Self;
    fn root(&self) -> Self::Inner;
    fn get_node_number(&self) -> usize;
    fn calculate_from_leafs(children: &[Self::Input]) -> Self::Inner;
    fn calculate_from_nodes(children: &[Self::Inner]) -> Self::Inner;
}

#[derive(Clone)]
pub struct MerkleTreeConfig<D: Digest, F: PrimeField> {
    leafs_per_node: usize,
    inner_children: usize,
    _digest: PhantomData<D>,
    _field: PhantomData<F>,
}

#[derive(Clone)]
pub struct MerkleTree<D: Digest, F: PrimeField> {
    leafs: Vec<F>,
    nodes: Vec<Hash<D>>,
    config: MerkleTreeConfig<D, F>,
    levels: usize,
}

impl<D: Digest, F: PrimeField> Tree for MerkleTree<D, F> {
    type Input = F;
    type Inner = Hash<D>;
    type Config = MerkleTreeConfig<D, F>;

    fn new(inputs: &[F], config: MerkleTreeConfig<D, F>) -> Self {
        let MerkleTreeConfig {
            leafs_per_node,
            inner_children,
            _digest,
            _field,
        } = config;

        let leaf_num = inputs.len();
        let node_num = leaf_num / leafs_per_node;
        let levels = match logarithm_of_two_k(node_num, inner_children) {
            Ok(log) => log + 1,
            Err(error_str) => panic!("{}", error_str),
        };

        assert_eq!(
            inner_children.pow((levels - 1) as u32) * leafs_per_node,
            leaf_num,
            "Tree is not full! input length must be a power of {inner_children}"
        );

        // number of nodes
        let numerator = 1 - (inner_children as i64).pow(levels as u32);
        let denominator = 1 - inner_children as i64;
        let node_num = (numerator / denominator) as usize;
        let mut nodes = Vec::with_capacity(node_num);
        println!("Number of nodes in the tree: {node_num}");

        let mut distance = leaf_num;
        for external in inputs.chunks(leafs_per_node) {
            let parent = Self::calculate_from_leafs(external);
            nodes.push(parent);
            distance -= leafs_per_node - 1;
        }

        let mut current_nodes = leaf_num / leafs_per_node;
        println!("current_nodes: {}", current_nodes);
        let nodes_left = node_num - current_nodes;
        for _ in 0..nodes_left {
            let idx = current_nodes - distance;
            let children = &nodes[idx..idx + inner_children];
            let parent = Self::calculate_from_nodes(children);
            nodes.push(parent);
            distance -= inner_children - 1;
            current_nodes += 1;
        }

        Self {
            leafs: inputs.to_vec(),
            nodes,
            config,
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
        for child in children.iter() {
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

impl<F: PrimeField, D: Digest> MerkleTree<D, F> {
    // TODO: better error handling
    fn get_neighbor_idx(&self, index: usize) -> Range<usize> {
        if index >= self.get_node_number() {
            panic!("index outside of tree length");
        }

        let remainder = index % self.config.inner_children;
        let start_idx = index - remainder;
        let end_idx = start_idx + self.config.inner_children;
        start_idx..end_idx
    }

    // TODO: better error handling
    fn get_parent_idx(&self, index: usize) -> usize {
        let root_idx = self.get_node_number() - 1;
        match index.cmp(&root_idx) {
            std::cmp::Ordering::Greater => panic!("index outside of tree length"),
            std::cmp::Ordering::Equal => panic!("index is root node"),
            std::cmp::Ordering::Less => {
                if index < self.leafs.len() {
                    self.leafs.len() + index / self.config.leafs_per_node
                } else {
                    index + (self.get_node_number() - index + 1) / self.config.inner_children
                }
            }
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
        let shifted_index = index - self.leafs.len();
        let neighbors_idx = self.get_neighbor_idx(shifted_index);
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
            println!("current idx: {}, parent: {}", current_idx, parent);
            current_idx = parent;
        }

        path
    }

    pub fn generate_proof(&self, leaf: &F) -> Result<MerklePath<D, F>, &str> {
        let leaf_index = self.get_leaf_index(leaf);
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

#[derive(Clone, Debug)]
pub struct MerklePath<D: Digest, F: PrimeField> {
    leaf_neighbours: Vec<F>,
    path: Vec<Vec<Hash<D>>>,
}

pub struct MerkleRoot<D: Digest>(pub Hash<D>);

impl<D: Digest> MerkleRoot<D> {
    pub fn check_proof<F: PrimeField>(&self, leaf: &F, proof: MerklePath<D, F>) -> bool {
        if !proof.leaf_neighbours.contains(leaf) {
            return false;
        };

        let mut previous = MerkleTree::<D, F>::calculate_from_leafs(&proof.leaf_neighbours);

        for (i, level) in proof.path.iter().enumerate() {
            if !level.contains(&previous) {
                return false;
            }

            previous = MerkleTree::<D, F>::calculate_from_nodes(level);
            println!("Merkle Verification: Round {i} done");
        }

        if previous == self.0 {
            return true;
        }
        assert_eq!(previous, self.0);
        false
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::Goldilocks;
    use sha2::Sha256;

    use std::panic;

    const TWO: MerkleTreeConfig<Sha256, Goldilocks> = MerkleTreeConfig {
        leafs_per_node: 2,
        inner_children: 2,
        _digest: PhantomData::<Sha256>,
        _field: PhantomData::<Goldilocks>,
    };

    const TWO_FOUR: MerkleTreeConfig<Sha256, Goldilocks> = MerkleTreeConfig {
        leafs_per_node: 4,
        inner_children: 2,
        _digest: PhantomData::<Sha256>,
        _field: PhantomData::<Goldilocks>,
    };

    const FOUR: MerkleTreeConfig<Sha256, Goldilocks> = MerkleTreeConfig {
        leafs_per_node: 4,
        inner_children: 4,
        _digest: PhantomData::<Sha256>,
        _field: PhantomData::<Goldilocks>,
    };

    const SIXTEEN: MerkleTreeConfig<Sha256, Goldilocks> = MerkleTreeConfig {
        leafs_per_node: 16,
        inner_children: 16,
        _digest: PhantomData::<Sha256>,
        _field: PhantomData::<Goldilocks>,
    };

    fn make_tree(config: MerkleTreeConfig<Sha256, Goldilocks>) -> MerkleTree<Sha256, Goldilocks> {
        let leafs: Vec<Goldilocks> = (0..16).map(Goldilocks::from).collect();
        MerkleTree::<Sha256, _>::new(&leafs, config)
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
            MerkleTree::<Sha256, _>::new(&leafs, TWO);
        });
        assert!(result.is_err(), "Tree of width: 2 to panic");
    }

    #[test]
    fn test_node_calculation() {
        let tree = make_tree(TWO);
        assert_eq!(tree.get_node_number(), 31);
        assert_eq!(tree.leafs.len(), 16);
        assert_eq!(tree.nodes.len(), 15);

        let tree = make_tree(TWO_FOUR);
        assert_eq!(tree.get_node_number(), 23);
        assert_eq!(tree.leafs.len(), 16);
        assert_eq!(tree.nodes.len(), 7);

        let tree = make_tree(FOUR);
        assert_eq!(tree.get_node_number(), 21);
        assert_eq!(tree.leafs.len(), 16);
        assert_eq!(tree.nodes.len(), 5);

        let tree = make_tree(SIXTEEN);
        assert_eq!(tree.get_node_number(), 17);
        assert_eq!(tree.leafs.len(), 16);
        assert_eq!(tree.nodes.len(), 1);
    }

    #[test]
    fn test_neighbor_index() {
        let tree = make_tree(TWO);
        assert_eq!(tree.get_neighbor_idx(4), 4..6);
        assert_eq!(tree.get_neighbor_idx(7), 6..8);

        let tree = make_tree(TWO_FOUR);
        assert_eq!(tree.get_neighbor_idx(4), 4..6);
        assert_eq!(tree.get_neighbor_idx(7), 6..8);

        let tree = make_tree(FOUR);
        assert_eq!(tree.get_neighbor_idx(4), 4..8);
        assert_eq!(tree.get_neighbor_idx(7), 4..8);

        let tree = make_tree(SIXTEEN);
        assert_eq!(tree.get_neighbor_idx(4), 0..16);
        assert_eq!(tree.get_neighbor_idx(7), 0..16);
    }

    // TODO:use macro to test all configs
    #[test]
    fn test_merkle_tree_parent_index() {
        let tree = make_tree(TWO);
        // first level test
        assert_eq!(tree.get_parent_idx(1), 16);
        assert_eq!(tree.get_parent_idx(4), 18);
        assert_eq!(tree.get_parent_idx(9), 20);
        assert_eq!(tree.get_parent_idx(13), 22);
        // second level ..
        assert_eq!(tree.get_parent_idx(16), 24);
        assert_eq!(tree.get_parent_idx(18), 25);
        assert_eq!(tree.get_parent_idx(20), 26);
        assert_eq!(tree.get_parent_idx(22), 27);
        // thrid level ..
        assert_eq!(tree.get_parent_idx(24), 28);
        assert_eq!(tree.get_parent_idx(25), 28);
        assert_eq!(tree.get_parent_idx(26), 29);
        assert_eq!(tree.get_parent_idx(27), 29);
        // fourth level ..
        assert_eq!(tree.get_parent_idx(28), 30);
        assert_eq!(tree.get_parent_idx(29), 30);

        let tree = make_tree(TWO_FOUR);
        // first level test
        assert_eq!(tree.get_parent_idx(1), 16);
        assert_eq!(tree.get_parent_idx(4), 17);
        assert_eq!(tree.get_parent_idx(9), 18);
        assert_eq!(tree.get_parent_idx(13), 19);
        // second level ..
        assert_eq!(tree.get_parent_idx(16), 20);
        assert_eq!(tree.get_parent_idx(17), 20);
        assert_eq!(tree.get_parent_idx(18), 21);
        assert_eq!(tree.get_parent_idx(19), 21);
        // third level ..
        assert_eq!(tree.get_parent_idx(20), 22);
        assert_eq!(tree.get_parent_idx(21), 22);

        // test that calling and index out of tree length panics
        let result = panic::catch_unwind(|| {
            tree.get_parent_idx(tree.get_node_number());
        });
        assert!(result.is_err());
    }

    // TODO: modify python script to reproduce hashing outcomes
    #[test]
    fn test_check_proof() {
        let tree = make_tree(TWO);
        let root = tree.root();

        let f_one = Goldilocks::from(7);
        let proof = tree.generate_proof(&f_one).unwrap();
        println!("Proof: {:?}", proof);
        assert!(MerkleRoot::<Sha256>(root).check_proof::<Goldilocks>(&f_one, proof));
    }
}
