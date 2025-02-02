use ark_ff::Field;
use digest::{generic_array::GenericArray, Digest, OutputSizeUser};
use std::marker::PhantomData;

type Hash<D> = GenericArray<u8, <D as OutputSizeUser>::OutputSize>;

// MerkleTree where all nodes are stored in memory
struct MerkleTree<D: Digest, F> {
    nodes: Vec<Hash<D>>,
    leaf_number: usize,
    input: PhantomData<F>,
}

fn is_power_of_2(number: usize) -> bool {
    number & (number - 1) == 0
}

impl<F: Field, D: Digest> MerkleTree<D, F> {
    pub fn get_leafs(&self) -> &[Hash<D>] {
        &self.nodes[0..self.leaf_number]
    }

    pub fn get_root(&self) -> &Hash<D> {
        let capacity = self.nodes.capacity();
        &self.nodes[capacity - 1]
    }

    fn calculate_path(&self, index: usize) -> Vec<Hash<D>> {
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

        path.into_iter().map(|i| self.nodes[i].clone()).collect()
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

    fn get_leaf_index(&self, node: F) -> Option<usize> {
        let digest = MerkleTree::<D, F>::hash_leaf(node);
        for (i, value) in self.get_leafs().iter().enumerate() {
            if digest == *value {
                return Some(i);
            }
        }

        None
    }

    fn hash_leaf(value: F) -> Hash<D> {
        D::digest(value.to_string())
    }

    fn hash_children(children: &[Hash<D>]) -> Hash<D> {
        let mut hasher = D::new();
        for c in children {
            hasher.update(c);
        }
        hasher.finalize()
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
            let hashed_leaf = MerkleTree::<D, F>::hash_leaf(leaf);
            nodes.push(hashed_leaf);
        }

        let mut node_dist = leaf_number;
        for parent_idx in leaf_number..capacity {
            let lhs_idx = parent_idx - node_dist;
            let rhs_idx = lhs_idx + 1;
            nodes.push(MerkleTree::<D, F>::hash_children(&[
                nodes[lhs_idx].clone(),
                nodes[rhs_idx].clone(),
            ]));
            node_dist -= 1;
        }

        Self {
            nodes,
            leaf_number,
            input: PhantomData,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_ff::fields::{MontBackend, MontConfig};
    use ark_ff::Fp;
    use sha2::Sha256;

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct Fp17Config;
    pub type Fp17 = Fp<MontBackend<Fp17Config, 1>, 1>;

    fn make_tree() -> MerkleTree<Sha256, Fp17> {
        let leafs: Vec<Fp17> = (0..8).map(|i| Fp17::from(i)).collect();
        MerkleTree::generate_tree(leafs)
    }

    #[test]
    #[should_panic]
    fn test_arbitrary_number_of_leafs() {
        let leafs = vec![Fp17::from(0), Fp17::from(1), Fp17::from(2)];
        MerkleTree::<Sha256, _>::generate_tree(leafs);
    }

    #[test]
    fn test_valid_merkle_tree() {
        let tree = make_tree();
        let bytes: [u8; 32] = [
            59, 130, 140, 79, 75, 72, 197, 212, 203, 85, 98, 164, 116, 236, 158, 47, 216, 213, 84,
            111, 174, 64, 233, 7, 50, 239, 99, 88, 146, 228, 39, 32,
        ];

        let hash_value: Hash<Sha256> = Hash::<Sha256>::clone_from_slice(&bytes);

        assert_eq!(tree.leaf_number, 8);
        assert_eq!(tree.nodes.len(), 15);
        assert_eq!(*tree.get_root(), hash_value);

        assert_eq!(tree.get_leaf_index(Fp17::from(2)), Some(2));
        assert_eq!(tree.get_leaf_index(Fp17::from(5)), Some(5));
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
        let path_raw = vec![
            *Hash::<Sha256>::from_slice(&[
                95, 236, 235, 102, 255, 200, 111, 56, 217, 82, 120, 108, 109, 105, 108, 121, 194,
                219, 194, 57, 221, 78, 145, 180, 103, 41, 215, 58, 39, 251, 87, 233,
            ]),
            *Hash::<Sha256>::from_slice(&[
                169, 245, 179, 171, 97, 226, 131, 87, 207, 205, 20, 226, 180, 35, 151, 248, 150,
                174, 234, 141, 105, 152, 209, 158, 109, 168, 85, 132, 225, 80, 210, 180,
            ]),
            *Hash::<Sha256>::from_slice(&[
                3, 2, 201, 111, 69, 171, 190, 173, 178, 56, 120, 51, 26, 155, 164, 6, 7, 139, 208,
                189, 93, 194, 2, 193, 2, 175, 123, 153, 134, 36, 159, 1,
            ]),
        ];
        assert_eq!(path1, path_raw);

        let path2 = tree.calculate_path(4);
        let path_raw = vec![
            *Hash::<Sha256>::from_slice(&[
                239, 45, 18, 125, 227, 123, 148, 43, 170, 208, 97, 69, 229, 75, 12, 97, 154, 31,
                34, 50, 123, 46, 187, 207, 190, 199, 143, 85, 100, 175, 227, 157,
            ]),
            *Hash::<Sha256>::from_slice(&[
                19, 72, 67, 175, 127, 200, 242, 153, 80, 177, 225, 223, 183, 196, 151, 82, 224,
                247, 183, 17, 180, 88, 238, 154, 227, 197, 202, 34, 1, 102, 214, 136,
            ]),
            *Hash::<Sha256>::from_slice(&[
                196, 120, 254, 173, 12, 137, 183, 149, 64, 99, 143, 132, 76, 136, 25, 217, 164, 40,
                23, 99, 175, 146, 114, 199, 243, 150, 135, 118, 182, 5, 35, 69,
            ]),
        ];
        assert_eq!(path2, path_raw);
    }
}
