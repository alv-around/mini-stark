use super::merkle::{MerklePath, MerkleTree};
use super::util::is_power_of_two;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::Digest;

struct Fri<D: Digest, F: PrimeField> {
    blowup_factor: usize,
    rounds: usize,
    round_state: FriRound<D, F>,
}

struct FriRound<D: Digest, F: PrimeField> {
    round: usize,
    domain: Radix2EvaluationDomain<F>,
    commit: MerkleTree<D, F>,
    poly: DensePolynomial<F>,
    alpha: Option<usize>,
    beta: Option<usize>,
    poly_star: DensePolynomial<F>,
    commit_star: MerkleTree<D, F>,
}

impl<D: Digest, F: PrimeField> Fri<D, F> {
    pub fn new(coeffs: Vec<F>, blowup_factor: usize) -> Self {
        let d = coeffs.len() - 1;
        let domain_size = d * blowup_factor;
        if !is_power_of_two(domain_size) {
            panic!("blowup factor and degree of polynomial must be a power of 2");
        }

        let rounds: usize = domain_size.trailing_zeros().try_into().unwrap();
        Self {
            blowup_factor,
            rounds,
            round_state: FriRound::<D, F>::new(coeffs, domain_size),
        }
    }
}

impl<D: Digest, F: PrimeField> FriRound<D, F> {
    fn new(coeffs: Vec<F>, domain_size: usize) -> Self {
        let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
        let poly = DensePolynomial::<F>::from_coefficients_vec(coeffs);

        let leafs = poly.evaluate_over_domain_by_ref(domain);
        let commit = MerkleTree::<D, _>::generate_tree(&leafs.evals);
        Self {
            poly,
            commit,
            domain,
            round: 0,
            alpha: None,
            beta: None,
            commit_star: MerkleTree::<D, F>::generate_tree(&vec![F::ONE]),
            poly_star: DensePolynomial::<F>::from_coefficients_vec(vec![]),
        }
    }

    pub fn verify(&self, proof1: MerklePath<D, F>, proof2: MerklePath<D, F>) -> bool {
        if self.beta.is_none() || self.alpha.is_none() {
            return false;
        }
        let beta = self.beta.unwrap();
        let alpha = self.alpha.unwrap();

        //
        let n_over_2 = self.commit_star.get_leaf_number();
        let x_point1 = self.domain.element(beta);
        let x_point2 = self.domain.element(n_over_2 + beta);
        assert_eq!(proof1.leaf, self.poly.evaluate(&x_point1));
        assert_eq!(proof2.leaf, self.poly.evaluate(&x_point2));

        // Linearity check
        let y1 = proof1.leaf;
        let y2 = proof2.leaf;
        let slope = (y2 - y1) / (x_point2 - x_point1);
        let expected_y = slope * (self.domain.element(alpha) - x_point1) + y1;
        let actual_y = self.poly_star.evaluate(&self.domain.element(2 * beta));
        assert_eq!(expected_y, actual_y, "Linearity check failed");

        // check that elements where rightly committed
        assert!(self.commit.check_proof(proof1));
        assert!(self.commit_star.check_proof(proof2));

        true
    }

    pub fn split_and_fold_phase(&self, alpha: F) -> (DensePolynomial<F>, MerkleTree<D, F>) {
        let (even_poly, odd_poly) = self.split_poly();
        let folded_poly_coeffs = self.fold_poly(alpha, even_poly, odd_poly);
        let merkle_tree = MerkleTree::<D, F>::generate_tree(&folded_poly_coeffs.to_vec());

        (folded_poly_coeffs, merkle_tree)
    }

    fn split_poly(&self) -> (DensePolynomial<F>, DensePolynomial<F>) {
        let (mut even, mut odd) = (Vec::<F>::new(), Vec::<F>::new());
        for (i, element) in self.poly.coeffs().to_vec().into_iter().enumerate() {
            if i % 2 == 0 {
                even.push(element);
            } else {
                odd.push(element);
            }
        }

        let (even_poly, odd_poly) = (
            DensePolynomial::from_coefficients_vec(even),
            DensePolynomial::from_coefficients_vec(odd),
        );
        (even_poly, odd_poly)
    }

    fn fold_poly(
        &self,
        alpha: F,
        even_poly: DensePolynomial<F>,
        odd_poly: DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        let folded_domain_length = self.domain.size() >> 1;
        let folded_domain = (0..folded_domain_length).map(|i| {
            let exponent: u64 = (i * 2).try_into().unwrap();
            self.domain.group_gen().pow(&[exponent])
        });

        let folded_coeffs = folded_domain
            .map(|point| even_poly.evaluate(&point) + alpha * odd_poly.evaluate(&point))
            .collect::<Vec<F>>();

        DensePolynomial::<F>::from_coefficients_vec(folded_coeffs)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::Goldilocks;
    use sha2::Sha256;

    #[test]
    fn test_fri_new() {
        let blowup_factor = 2usize;
        let coeffs = vec![Goldilocks::from(1), Goldilocks::from(5)];
        let fri = Fri::<Sha256, _>::new(coeffs, blowup_factor);

        assert_eq!(fri.rounds, 1usize);
        assert_eq!(fri.blowup_factor, 2usize);
    }
}
