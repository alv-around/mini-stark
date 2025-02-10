use ark_ff::Field;

#[derive(Debug)]
struct MerkleRoot<F>(F);

// MerkleTree where all nodes are stored in memory
struct MerkleTree<F> {
    nodes: Vec<F>,
    leaf_number: usize,
}

fn is_power_of_2(number: usize) -> bool {
    number & (number - 1) == 0
}

impl<F: Field> MerkleTree<F> {
    pub fn get_leafs(&self) -> &[F] {
        &self.nodes[0..self.leaf_number]
    }

    pub fn get_root(&self) -> MerkleRoot<F> {
        let capacity = self.nodes.capacity();
        MerkleRoot(self.nodes[capacity - 1])
    }

    fn calculate_path(&self, index: usize) -> Vec<F> {
        // if number is a power of 2, trailing_zeros is same as log2
        let levels = self.leaf_number.trailing_zeros() as usize;
        let mut path = Vec::with_capacity(levels);
        let mut current_idx = index;
        for _ in 0..levels {
            let neighbor: usize = self.get_neighbor_idx(current_idx);
            let parent = self.get_parent_idx(current_idx);
            path.push(neighbor);
            current_idx = parent;
        }

        path.into_iter().map(|i| self.nodes[i]).collect()
    }

    fn get_neighbor_idx(&self, index: usize) -> usize {
        let is_uneven = index % 2;
        index + 1 - (2 * is_uneven)
    }

    fn get_parent_idx(&self, index: usize) -> usize {
        let mut remaining = self.nodes.len() - index;
        if remaining % 2 == 1 {
            remaining += 1;
        }
        (remaining >> 1) + index
    }

    fn get_leaf_index(&self, node: &F) -> Option<usize> {
        for (i, value) in self.get_leafs().iter().enumerate() {
            if node == value {
                return Some(i);
            }
        }

        None
    }

    pub fn generate_tree(leafs: Vec<F>) -> Self {
        let leaf_number = leafs.len();
        assert!(
            is_power_of_2(leaf_number),
            "Number of tree leafs must be a power of 2"
        );

        let capacity = (leaf_number * 2) - 1;
        let mut nodes = Vec::with_capacity(capacity);
        for leaf in leafs {
            nodes.push(leaf);
        }

        let mut node_dist = leaf_number;
        for parent_idx in leaf_number..capacity {
            let lhs_idx = parent_idx - node_dist;
            let rhs_idx = lhs_idx + 1;
            nodes.push(nodes[lhs_idx] + nodes[rhs_idx]);
            node_dist -= 1;
        }

        Self { nodes, leaf_number }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_ff::fields::{MontBackend, MontConfig};
    use ark_ff::Fp;

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct Fp17Config;
    pub type Fp17 = Fp<MontBackend<Fp17Config, 1>, 1>;

    fn make_tree() -> MerkleTree<Fp17> {
        let leafs: Vec<Fp17> = (0..8).map(|i| Fp17::from(i)).collect();
        MerkleTree::generate_tree(leafs)
    }

    #[test]
    #[should_panic]
    fn test_arbitrary_number_of_leafs() {
        let leafs = vec![Fp17::from(0), Fp17::from(1), Fp17::from(2)];
        MerkleTree::generate_tree(leafs);
    }

    #[test]
    fn test_valid_merkle_tree() {
        let tree = make_tree();

        assert_eq!(tree.leaf_number, 8);
        assert_eq!(tree.nodes.len(), 15);
        assert_eq!(tree.get_root().0, Fp17::from(28));

        assert_eq!(tree.get_leaf_index(&tree.get_leafs()[2]), Some(2));
        assert_eq!(tree.get_leaf_index(&tree.get_leafs()[5]), Some(5));
    }

    #[test]
    fn test_neighbor_index() {
        let tree = make_tree();

        assert_eq!(tree.get_neighbor_idx(4), 5);
        assert_eq!(tree.get_neighbor_idx(7), 6);
    }

    #[test]
    fn test_merkle_tree_parent_index() {
        let tree = make_tree();

        // assert parents
        assert_eq!(tree.get_parent_idx(0), 8);
        assert_eq!(tree.get_parent_idx(3), 9);
        assert_eq!(tree.get_parent_idx(10), 13);
        assert_eq!(tree.get_parent_idx(12), 14);
    }

    #[test]
    fn test_merkle_path() {
        let tree = make_tree();

        let path1 = tree.calculate_path(1);
        assert_eq!(path1, vec![Fp17::from(0), Fp17::from(5), Fp17::from(22)]);

        let path2 = tree.calculate_path(5);
        assert_eq!(path2, vec![Fp17::from(4), Fp::from(13), Fp::from(6)]);
    }
}
