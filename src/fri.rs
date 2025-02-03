use super::merkle::MerkleTree;
use super::util::{get_nth_root, is_power_of_two};
use ark_ff::PrimeField;
// use ark_poly::domain::{EvaluationDomain, GeneralEvaluationDomain};
// use ark_poly::evaluations::univariate::Evaluations;
use ark_poly::univariate::DensePolynomial;
use ark_poly::Polynomial;
use digest::Digest;

struct FRI<D: Digest, F: PrimeField> {
    nth_root: F,
    merkle_tree: MerkleTree<D, F>,
    blowup_factor: usize,
    poly: DensePolynomial<F>,
    leafs: Vec<F>,
}

impl<D: Digest, F: PrimeField> FRI<D, F> {
    pub fn new(poly: DensePolynomial<F>, blowup_factor: usize) -> Self {
        if !is_power_of_two(blowup_factor) {
            panic!("blowup factor must be a power of 2");
        }

        // TODO: better handle of blowup_factor and k. They must be a power of 2
        let k = poly.degree();
        let domain_size = k * blowup_factor;
        let nth_root = get_nth_root::<F>(domain_size as u64);

        let mut domain = vec![F::ONE];
        for i in 1..domain_size {
            let result = nth_root * domain[i - 1];
            domain.push(result);
        }

        // TODO: find a better way to calculate leafs??
        let leafs: Vec<F> = domain.iter().map(|i| poly.evaluate(i)).collect();
        // let domain = EvaluationDomain::<F>::new_coset(domain_size, nth_root).unwrap();
        // let leafs: Evaluations<F, GeneralEvaluationDomain<F>> =
        //     poly.evaluate_over_domain_by_ref(domain);

        let merkle_tree = MerkleTree::<D, _>::generate_tree(&leafs);

        Self {
            poly,
            blowup_factor,
            merkle_tree,
            nth_root,
            leafs,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::Goldilocks;
    use ark_ff::Field;
    use ark_poly::DenseUVPolynomial;
    use sha2::Sha256;

    #[test]
    fn test_fri_new() {
        let blowup_factor = 2usize;
        let coeffs = vec![Goldilocks::from(1), Goldilocks::from(5)];
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let fri = FRI::<Sha256, _>::new(poly.clone(), blowup_factor);
        let leafs = fri.leafs;

        assert_eq!(leafs[0], poly.evaluate(&Goldilocks::ONE));
        assert_eq!(leafs[1], poly.evaluate(&fri.nth_root));
    }
}
