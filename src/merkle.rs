use crate::{error::MerkleProofError, util::logarithm_of_two_k, Hash};
use ark_ff::FftField;
use digest::Digest;
use log::trace;
use std::marker::PhantomData;

/// Trait defining the common interface for tree structures
pub trait Tree {
    /// Type of the leaf nodes (input elements)
    type Input;
    /// Type of the inner nodes (hashes)
    type Inner;
    /// Configuration type for the tree
    type Config;

    /// Creates a new tree from input elements and configuration
    fn new(inputs: &[Self::Input], config: Self::Config) -> Self;

    /// Returns the root hash of the tree
    fn root(&self) -> Self::Inner;

    /// Returns the total number of nodes in the tree (leaves + inner nodes)
    fn get_node_number(&self) -> usize;

    /// Computes a parent node from leaf children
    fn calculate_from_leafs(children: &[Self::Input]) -> Self::Inner;

    /// Computes a parent node from inner node children
    fn calculate_from_nodes(children: &[Self::Inner]) -> Self::Inner;
}

/// Configuration for a Merkle Tree specifying its structure
#[derive(Clone)]
pub struct MerkleTreeConfig<D: Digest, F: FftField> {
    /// Number of leaves grouped under each lowest-level inner node
    pub leafs_per_node: usize,
    /// Number of children for each inner node (except lowest level)
    pub inner_children: usize,
    /// Phantom data for the digest type
    pub _digest: PhantomData<D>,
    /// Phantom data for the field type
    pub _field: PhantomData<F>,
}

/// Merkle Tree implementation supporting configurable branching factors
/// TREE STRUCTURE EXAMPLE (leafs_per_node=4, inner_children=2, 16 leaves)
/// Level 3:                  [Node23/Root]
/// (Inner)                  /       \
//; Level 2:          [Node21]       [Node22]
///                 /    \          /    \
/// Level 1:    [Node17] [Node18] [Node19] [node20]
/// (Inner)      / \      / \      / \       / \
/// Level 0:  0 1 2 3  4 5 6 7  8 9 10 11 12 13 14 15
/// (Input)
#[derive(Clone)]
pub struct MerkleTree<D: Digest, F: FftField> {
    /// Vector containing all leaf values
    leafs: Vec<F>,
    /// Vector containing all inner nodes (hashes)
    nodes: Vec<Hash<D>>,
    /// Tree configuration
    config: MerkleTreeConfig<D, F>,
    /// Number of levels in the tree (including leaf level)
    levels: usize,
}

impl<D: Digest, F: FftField> Tree for MerkleTree<D, F> {
    type Input = F;
    type Inner = Hash<D>;
    type Config = MerkleTreeConfig<D, F>;

    /// Constructs a new Merkle Tree from input elements and configuration
    ///
    /// # Arguments
    /// * `inputs` - Slice of leaf values
    /// * `config` - Tree configuration specifying branching factors
    ///
    /// # Panics
    /// - If input length is not divisible by leafs_per_node
    /// - If the resulting tree cannot be made full with given branching factors
    fn new(inputs: &[F], config: MerkleTreeConfig<D, F>) -> Self {
        let MerkleTreeConfig {
            leafs_per_node,
            inner_children,
            _digest,
            _field,
        } = config;

        let leaf_num = inputs.len();
        let node_num = leaf_num / leafs_per_node;

        // Calculate number of levels in the tree
        let levels = match logarithm_of_two_k(node_num, inner_children) {
            Ok(log) => log + 1, // +1 to include the leaf level
            Err(error_str) => panic!("{}", error_str),
        };

        // Validate input length is compatible with configuration
        assert_eq!(leaf_num % leafs_per_node, 0);
        assert_eq!(
            inner_children.pow((levels - 1) as u32),
            leaf_num / leafs_per_node,
            "Tree is not full! input length must be a power of {inner_children}"
        );

        // TREE STRUCTURE EXAMPLE (leafs_per_node=2, inner_children=2, 16 leaves)
        // Level 3:                  [Root]
        //                          /       \
        // Level 2:          [Node28]       [Node29]
        //                 /    \          /    \
        // Level 1:    [Node24] [Node25] [Node26] [Node27]
        //             / \      / \      / \      / \
        // Level 0:  0 1 2 3  4 5 6 7  8 9 10 11 12 13 14 15

        // Calculate total number of nodes in the tree using geometric series formula
        let numerator = 1 - (inner_children as i64).pow(levels as u32);
        let denominator = 1 - inner_children as i64;
        let node_num = (numerator / denominator) as usize;
        let mut nodes = Vec::with_capacity(node_num);
        trace!("Number of inner nodes in the tree: {node_num}");

        // First pass: hash groups of leaves to create first level of inner nodes
        let mut distance = leaf_num;
        for external in inputs.chunks(leafs_per_node) {
            let parent = Self::calculate_from_leafs(external);
            nodes.push(parent);
            distance -= leafs_per_node - 1;
        }

        // Second pass: build upper levels of the tree by hashing inner nodes
        let mut current_nodes = leaf_num / leafs_per_node;
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

    /// Returns the root hash of the tree
    fn root(&self) -> Hash<D> {
        let root = self.nodes.last().unwrap().clone();
        root
    }

    /// Returns total number of nodes in the tree (leaves + inner nodes)
    fn get_node_number(&self) -> usize {
        self.leafs.len() + self.nodes.len()
    }

    /// Computes hash of a group of leaf values
    fn calculate_from_leafs(children: &[F]) -> Hash<D> {
        let mut hasher = D::new();
        for child in children.iter() {
            hasher.update(child.to_string());
        }
        hasher.finalize()
    }

    /// Computes hash of a group of inner node hashes
    fn calculate_from_nodes(children: &[Hash<D>]) -> Hash<D> {
        let mut hasher = D::new();
        for child in children {
            hasher.update(child)
        }
        hasher.finalize()
    }
}

impl<F: FftField, D: Digest> MerkleTree<D, F> {
    /// Gets the index of a node's parent in the nodes vector
    ///
    /// # Arguments
    /// * `index` - Index of the child node
    ///
    /// # Returns
    /// Result containing parent index or error if index is invalid
    fn get_parent_idx(&self, index: usize) -> Result<usize, MerkleProofError> {
        let root_idx = self.get_node_number() - 1;
        match index.cmp(&root_idx) {
            std::cmp::Ordering::Greater => Err(MerkleProofError::OutOfRangeError {
                msg: "index outside of tree length",
            }),
            std::cmp::Ordering::Equal => Err(MerkleProofError::OutOfRangeError {
                msg: "index is root node",
            }),
            std::cmp::Ordering::Less => {
                if index < self.leafs.len() {
                    // For leaf nodes, parent is in first level of inner nodes
                    Ok(self.leafs.len() + index / self.config.leafs_per_node)
                } else {
                    // For inner nodes, calculate parent position
                    Ok(index + (self.get_node_number() - index + 1) / self.config.inner_children)
                }
            }
        }
    }

    /// Finds the index of a leaf value in the tree
    ///
    /// # Arguments
    /// * `node` - Leaf value to find
    ///
    /// # Returns
    /// Result containing leaf index or error if not found
    fn get_leaf_index(&self, node: &F) -> Result<usize, MerkleProofError> {
        for (i, value) in self.leafs.iter().enumerate() {
            if *node == *value {
                return Ok(i);
            }
        }
        Err(MerkleProofError::LeafNotFound {
            msg: "leaf is not included in the tree",
        })
    }

    /// Gets the neighboring leaf values for a given leaf index
    ///
    /// These are all leaves that share the same parent in the lowest level
    fn get_leaf_neighbours(&self, index: usize) -> Vec<F> {
        let num_neighbors = self.config.leafs_per_node;
        let remainder = index % num_neighbors;
        let start_idx = index - remainder;
        let end_idx = start_idx + num_neighbors;
        self.leafs[start_idx..end_idx].to_vec()
    }

    /// Gets the neighboring inner nodes for a given inner node index
    ///
    /// These are all nodes that share the same parent at the given level
    fn get_inner_neighbours(&self, index: usize) -> Vec<Hash<D>> {
        let shifted_index = index - self.leafs.len();
        let num_neighbors = self.config.inner_children;
        let remainder = shifted_index % num_neighbors;
        let start_idx = shifted_index - remainder;
        let end_idx = start_idx + num_neighbors;
        self.nodes[start_idx..end_idx].to_vec()
    }

    /// Calculates the authentication path for a given leaf index
    ///
    /// The path consists of all sibling groups needed to reconstruct the root hash
    fn calculate_path(&self, index: usize) -> Result<Vec<Vec<Hash<D>>>, MerkleProofError> {
        let mut path = Vec::new();
        let mut current_idx = index;
        for _ in 1..self.levels {
            let neighbor = self.get_inner_neighbours(current_idx);
            path.push(neighbor);

            let parent = self.get_parent_idx(current_idx)?;
            current_idx = parent;
        }

        Ok(path)
    }

    /// Generates a Merkle proof for a given leaf value
    ///
    /// The proof contains:
    /// - The neighboring leaves at the lowest level
    /// - The sibling groups at each level up to the root
    pub fn generate_proof(&self, leaf: &F) -> Result<MerklePath<D, F>, MerkleProofError> {
        let leaf_index = self.get_leaf_index(leaf)?;

        let leaf_neighbours = self.get_leaf_neighbours(leaf_index);
        let leaf_parent = self.get_parent_idx(leaf_index)?;
        trace!(
            "generating merkle proof for leaf: {} leaf parent: {} leaf_neighbours: {:?}",
            leaf_index,
            leaf_parent,
            leaf_neighbours
        );
        let path = self.calculate_path(leaf_parent)?;
        Ok(MerklePath {
            leaf_neighbours,
            path,
        })
    }
}

/// Structure representing a Merkle proof
#[derive(Clone, Debug)]
pub struct MerklePath<D: Digest, F: FftField> {
    /// Leaf values that share the same parent as the target leaf
    pub(super) leaf_neighbours: Vec<F>,
    /// Sibling groups at each level needed to reconstruct the root
    path: Vec<Vec<Hash<D>>>,
}

/// Structure representing a Merkle root hash
#[derive(Debug)]
pub struct MerkleRoot<D: Digest>(pub Hash<D>);

impl<D: Digest> MerkleRoot<D> {
    /// Verifies a Merkle proof against this root
    ///
    /// # Arguments
    /// * `proof` - The Merkle proof to verify
    ///
    /// # Returns
    /// true if the proof is valid, false otherwise
    pub fn check_proof<F: FftField>(&self, proof: MerklePath<D, F>) -> bool {
        // Start by hashing the leaf group
        let mut previous = MerkleTree::<D, F>::calculate_from_leafs(&proof.leaf_neighbours);
        trace!(
            "Checking merkle proof for root ({:?}): {:?} | {:?}",
            self.0,
            proof.leaf_neighbours,
            proof.path
        );

        // For each level in the path, verify the hash chain
        for (i, level) in proof.path.iter().enumerate() {
            if !level.contains(&previous) {
                trace!("Merkle proof failed at level {}", i);
                return false;
            }

            previous = MerkleTree::<D, F>::calculate_from_nodes(level);
        }

        // Final check against the root hash
        if previous == self.0 {
            return true;
        }
        trace!("Merkle proof root does not match MerkleRoot");
        false
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::GoldilocksFp;
    use sha2::Sha256;
    use std::panic;

    // Test configurations
    const TWO: MerkleTreeConfig<Sha256, GoldilocksFp> = MerkleTreeConfig {
        leafs_per_node: 2,
        inner_children: 2,
        _digest: PhantomData::<Sha256>,
        _field: PhantomData::<GoldilocksFp>,
    };

    const TWO_FOUR: MerkleTreeConfig<Sha256, GoldilocksFp> = MerkleTreeConfig {
        leafs_per_node: 4,
        inner_children: 2,
        _digest: PhantomData::<Sha256>,
        _field: PhantomData::<GoldilocksFp>,
    };

    const FOUR: MerkleTreeConfig<Sha256, GoldilocksFp> = MerkleTreeConfig {
        leafs_per_node: 4,
        inner_children: 4,
        _digest: PhantomData::<Sha256>,
        _field: PhantomData::<GoldilocksFp>,
    };

    const SIXTEEN: MerkleTreeConfig<Sha256, GoldilocksFp> = MerkleTreeConfig {
        leafs_per_node: 16,
        inner_children: 16,
        _digest: PhantomData::<Sha256>,
        _field: PhantomData::<GoldilocksFp>,
    };

    fn make_tree(
        config: MerkleTreeConfig<Sha256, GoldilocksFp>,
    ) -> MerkleTree<Sha256, GoldilocksFp> {
        let leafs: Vec<GoldilocksFp> = (0..16).map(GoldilocksFp::from).collect();
        MerkleTree::<Sha256, _>::new(&leafs, config)
    }

    #[test]
    fn test_panic_for_not_full_trees() {
        let leafs = vec![
            GoldilocksFp::from(0),
            GoldilocksFp::from(1),
            GoldilocksFp::from(2),
        ];

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
    fn test_merkle_tree_parent_index() {
        let tree = make_tree(TWO);
        // first level test
        assert_eq!(tree.get_parent_idx(1).unwrap(), 16);
        assert_eq!(tree.get_parent_idx(4).unwrap(), 18);
        assert_eq!(tree.get_parent_idx(9).unwrap(), 20);
        assert_eq!(tree.get_parent_idx(13).unwrap(), 22);
        // second level ..
        assert_eq!(tree.get_parent_idx(16).unwrap(), 24);
        assert_eq!(tree.get_parent_idx(18).unwrap(), 25);
        assert_eq!(tree.get_parent_idx(20).unwrap(), 26);
        assert_eq!(tree.get_parent_idx(22).unwrap(), 27);
        // thrid level ..
        assert_eq!(tree.get_parent_idx(24).unwrap(), 28);
        assert_eq!(tree.get_parent_idx(25).unwrap(), 28);
        assert_eq!(tree.get_parent_idx(26).unwrap(), 29);
        assert_eq!(tree.get_parent_idx(27).unwrap(), 29);
        // fourth level ..
        assert_eq!(tree.get_parent_idx(28).unwrap(), 30);
        assert_eq!(tree.get_parent_idx(29).unwrap(), 30);

        let tree = make_tree(TWO_FOUR);
        // first level test
        assert_eq!(tree.get_parent_idx(1).unwrap(), 16);
        assert_eq!(tree.get_parent_idx(4).unwrap(), 17);
        assert_eq!(tree.get_parent_idx(9).unwrap(), 18);
        assert_eq!(tree.get_parent_idx(13).unwrap(), 19);
        // second level ..
        assert_eq!(tree.get_parent_idx(16).unwrap(), 20);
        assert_eq!(tree.get_parent_idx(17).unwrap(), 20);
        assert_eq!(tree.get_parent_idx(18).unwrap(), 21);
        assert_eq!(tree.get_parent_idx(19).unwrap(), 21);
        // third level ..
        assert_eq!(tree.get_parent_idx(20).unwrap(), 22);
        assert_eq!(tree.get_parent_idx(21).unwrap(), 22);

        // test that calling and index out of tree length panics
        let result = tree.get_parent_idx(tree.get_node_number());
        assert!(result.is_err());
    }

    #[test]
    fn test_check_proof() {
        let tree = make_tree(TWO);
        let root = tree.root();

        let f_one = GoldilocksFp::from(7);
        let proof = tree.generate_proof(&f_one).unwrap();
        assert!(proof.leaf_neighbours.contains(&f_one));
        assert_eq!(proof.path.len(), 3);
        assert!(MerkleRoot::<Sha256>(root).check_proof::<GoldilocksFp>(proof));

        let tree = make_tree(TWO_FOUR);
        let root = tree.root();

        let proof = tree.generate_proof(&f_one).unwrap();
        assert!(proof.leaf_neighbours.contains(&f_one));
        assert_eq!(proof.path.len(), 2);
        assert!(MerkleRoot::<Sha256>(root).check_proof::<GoldilocksFp>(proof));
    }
}
