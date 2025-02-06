use super::merkle::{MerklePath, MerkleTree};
use super::util::is_power_of_two;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::Digest;

struct Fri<D: Digest, F: PrimeField> {
    rounds: usize,
    round_state: FriRound<D, F>,
}

struct FriRound<D: Digest, F: PrimeField> {
    round: usize,
    blowup_factor: usize,
    domain: Radix2EvaluationDomain<F>,
    commit: MerkleTree<D, F>,
    poly: DensePolynomial<F>,
    alpha: Option<usize>,
    beta: Option<usize>,
    poly_star: DensePolynomial<F>,
    commit_star: MerkleTree<D, F>,
    domain_star: Radix2EvaluationDomain<F>,
}

impl<D: Digest, F: PrimeField> Fri<D, F> {
    pub fn new(poly: DensePolynomial<F>, blowup_factor: usize) -> Self {
        let d = poly.degree();
        let domain_size = d * blowup_factor;
        if !is_power_of_two(domain_size) {
            panic!("blowup factor and degree of polynomial must be a power of 2");
        }

        let rounds: usize = domain_size.trailing_zeros().try_into().unwrap();
        Self {
            rounds,
            round_state: FriRound::<D, F>::new(poly, blowup_factor),
        }
    }
}

impl<D: Digest, F: PrimeField> FriRound<D, F> {
    fn new(poly: DensePolynomial<F>, blowup_factor: usize) -> Self {
        let (commit, domain) = FriRound::<D, F>::codeword_commit(&poly, blowup_factor);
        Self {
            blowup_factor,
            round: 0,
            alpha: None,
            beta: None,
            poly,
            commit,
            domain,
            commit_star: MerkleTree::<D, F>::generate_tree(&vec![F::ONE]),
            poly_star: DensePolynomial::<F>::from_coefficients_vec(vec![]),
            domain_star: Radix2EvaluationDomain::<F>::new(1).unwrap(),
        }
    }

    pub fn verify(&self, proofs: [MerklePath<D, F>; 3]) -> bool {
        if self.beta.is_none() || self.alpha.is_none() {
            return false;
        }
        let beta = self.beta.unwrap();

        // xs
        let n_over_2 = self.commit_star.get_leaf_number();
        let x_point1 = self.domain.element(beta);
        let x_point2 = self.domain.element(n_over_2 + beta);
        let x_point3 = self.domain.element(self.alpha.unwrap());

        // FIXME: solve linearity check
        // Linearity check
        let y1 = proofs[0].leaf;
        let y2 = proofs[1].leaf;
        let slope = (y2 - y1) / (x_point2 - x_point1);
        let expected_y = slope * (x_point3 - x_point1) + y1;
        let actual_y = proofs[2].leaf;
        // assert_eq!(expected_y, actual_y, "Linearity check failed");

        // check that elements where rightly committed
        let [proof1, proof2, proof3] = proofs;
        assert!(self.commit.check_proof(proof1));
        assert!(self.commit.check_proof(proof2));
        assert!(self.commit_star.check_proof(proof3));

        true
    }

    pub fn prove(&mut self, beta: usize) -> Result<[MerklePath<D, F>; 3], &str> {
        if self.alpha.is_none() || self.beta.is_some() {
            return Err("wrong time");
        }
        self.beta = Some(beta);

        let x_point1 = self.domain.element(beta);
        let y_point1 = self.poly.evaluate(&x_point1);
        let x_point2 = self.domain.element(self.domain_star.size() + beta);
        let y_point2 = self.poly.evaluate(&x_point2);
        let x_point3 = self.domain.element(2 * beta);
        let y_point3 = self.poly_star.evaluate(&x_point3);

        let proof1 = self.commit.generate_proof(y_point1).unwrap();
        let proof2 = self.commit.generate_proof(y_point2).unwrap();
        let proof3 = self.commit_star.generate_proof(y_point3).unwrap();

        Ok([proof1, proof2, proof3])
    }

    fn split_and_fold(&mut self, alpha: usize) {
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

        self.poly_star = even_poly
            + odd_poly.naive_mul(&DensePolynomial::from_coefficients_slice(&[F::from(
                alpha as u64,
            )]));
        self.alpha = Some(alpha);
        let (commit_star, domain_star) =
            FriRound::<D, F>::codeword_commit(&self.poly_star, self.blowup_factor);
        self.commit_star = commit_star;
        self.domain_star = domain_star;
    }

    fn codeword_commit(
        poly: &DensePolynomial<F>,
        blowup_factor: usize,
    ) -> (MerkleTree<D, F>, Radix2EvaluationDomain<F>) {
        let domain_size = blowup_factor * poly.degree() + 1;
        let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
        let leafs = poly.evaluate_over_domain_by_ref(domain);
        (MerkleTree::<D, _>::generate_tree(&leafs.evals), domain)
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
        let fri = Fri::<Sha256, _>::new(
            DensePolynomial::from_coefficients_vec(coeffs),
            blowup_factor,
        );

        assert_eq!(fri.rounds, 1usize);
    }

    #[test]
    fn test_split_and_fold() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_slice(&coeffs);
        let mut fri = FriRound::<Sha256, _>::new(poly, blowup_factor);

        let alpha = 1usize;
        let coeffs_result = vec![Goldilocks::from(1), Goldilocks::from(5)];
        fri.split_and_fold(alpha);
        assert_eq!(fri.poly_star.coeffs(), coeffs_result);
    }

    #[test]
    fn test_domain_folding() {
        let domain_size = 32usize;
        let folded_size = domain_size >> 1;
        assert_eq!(folded_size, 16);

        let domain = Radix2EvaluationDomain::<Goldilocks>::new(domain_size).unwrap();
        let folded_domain = Radix2EvaluationDomain::<Goldilocks>::new(folded_size).unwrap();
        assert_eq!(domain.element(2), folded_domain.element(1));
    }

    #[test]
    fn test_fri_round() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_slice(&coeffs);
        let mut fri = FriRound::<Sha256, _>::new(poly, blowup_factor);

        let alpha = 1usize;
        let beta = 1usize;
        fri.split_and_fold(alpha);
        let proof = fri.prove(beta).unwrap();
        assert!(fri.verify(proof))
    }
}
