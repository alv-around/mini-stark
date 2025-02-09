use super::merkle::{MerklePath, MerkleRoot, MerkleTree};
use super::util::is_power_of_two;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::Digest;

struct Fri<D: Digest, F: PrimeField> {
    rounds: usize,
    blowup_factor: usize,
    poly: DensePolynomial<F>,
    domain: Radix2EvaluationDomain<F>,
    commit: MerkleTree<D, F>,
    round_state: Vec<FriRound<D, F>>,
    alpha: Option<usize>,
    beta: Option<usize>,
}

struct FriProof<D: Digest, F: PrimeField> {
    queries: Vec<MerklePath<D, F>>,
}

impl<D: Digest, F: PrimeField> Fri<D, F> {
    pub fn new(poly: DensePolynomial<F>, blowup_factor: usize) -> Self {
        let d = poly.degree();
        let domain_size = (d + 1) * blowup_factor;
        if !is_power_of_two(domain_size) {
            panic!("blowup factor and degree of polynomial must be a power of 2");
        }

        let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
        let commit = FriRound::<D, F>::codeword_commit(&poly, domain.clone());

        let rounds: usize = (domain_size.trailing_zeros() - 1).try_into().unwrap();
        Self {
            rounds,
            blowup_factor,
            poly,
            round_state: Vec::new(),
            commit,
            domain,
            alpha: None,
            beta: None,
        }
    }

    pub fn commit_phase(&mut self, alpha: usize) -> Vec<MerkleRoot<D>> {
        self.alpha = Some(alpha);
        let mut oracles = Vec::new();
        let mut previous_poly = self.poly.clone();
        let mut previous_alpha = self.alpha.unwrap();
        for i in 0..self.rounds {
            let round = FriRound::<D, _>::new(previous_poly, previous_alpha, self.blowup_factor, i);
            previous_poly = round.poly.clone();
            previous_alpha = round.alpha;
            oracles.push(round.commit.get_root());
            self.round_state.push(round);
        }

        oracles
    }

    pub fn query_phase(&mut self, beta: usize) -> Result<FriProof<D, F>, &str> {
        if self.alpha.is_none() || self.beta.is_some() {
            return Err("wrong time");
        }
        self.beta = Some(beta);
        let alpha = self.alpha.unwrap();
        let mut queries = Vec::new();

        let mut previous_poly = &self.poly;
        let mut previous_commit = &self.commit;
        let mut previous_domain = &self.domain;
        for round in self.round_state.iter() {
            let x1 = previous_domain.element(beta);
            let x2 = previous_domain.element(round.domain.size() + beta);
            let x3 = self.domain.element(beta * beta);
            assert_eq!(x1, x3);

            let y1 = previous_poly.evaluate(&x1);
            let y2 = previous_poly.evaluate(&x2);
            let y3 = round.poly.evaluate(&x3);

            let z = (y1 + y2) / F::from(2)
                + F::from(alpha as u64) * (y1 - y2) / F::from(2 * beta as u64);

            assert_eq!(z, y3);

            // no need to generate proof for x3 as it will be done in the next round as x1 == x3
            let proof1 = previous_commit.generate_proof(y1).unwrap();
            let proof2 = previous_commit.generate_proof(y2).unwrap();
            queries.push(proof1);
            queries.push(proof2);

            previous_poly = &round.poly;
            previous_commit = &round.commit;
            previous_domain = &round.domain;
        }

        Ok(FriProof { queries })
    }

    pub fn verify(&self, proofs: FriProof<D, F>) -> bool {
        if self.beta.is_none() || self.alpha.is_none() {
            return false;
        }
        let beta = self.beta.unwrap();

        // // xs
        // let x1 = self.domain.element(beta);
        // let x2 = self.domain.element(self.domain_round.size() + beta);
        // let x3 = self.domain.element(2 * beta);
        //
        // // check consistency_proof
        // let [consistency1, consistency2, consistency3] = consistency_proofs;
        // let vanishing1 = FriRound::<D, _>::build_vanishing_poly(x1);
        // consistency1.check(vanishing1);
        // let vanishing2 = FriRound::<D, _>::build_vanishing_poly(x2);
        // consistency2.check(vanishing2);
        // let vanishing3 = FriRound::<D, _>::build_vanishing_poly(x3);
        // // consistency3.check(vanishing3);
        //
        // // FIXME: solve linearity check
        // // Linearity check
        // let y1 = proofs[0].leaf;
        // let y2 = proofs[1].leaf;
        // let slope = (y2 - y1) / (x2 - x_point1);
        // let expected_y = slope * (x3 - x_point1) + y1;
        // let actual_y = proofs[2].leaf;
        // // assert_eq!(expected_y, actual_y, "Linearity check failed");
        //
        // // check that elements where rightly committed
        // let [proof1, proof2, proof3] = proofs;
        // assert!(self.commit.check_proof(proof1));
        // assert!(self.commit.check_proof(proof2));
        // assert!(self.commit_star.check_proof(proof3));
        //
        false
    }
}

struct FriRound<D: Digest, F: PrimeField> {
    round: usize,
    alpha: usize,
    poly: DensePolynomial<F>,
    commit: MerkleTree<D, F>,
    domain: Radix2EvaluationDomain<F>,
}

impl<D: Digest, F: PrimeField> FriRound<D, F> {
    fn new(
        previous_poly: DensePolynomial<F>,
        alpha: usize,
        blowup_factor: usize,
        round: usize,
    ) -> Self {
        let poly = FriRound::<D, F>::split_and_fold(previous_poly, alpha);
        let domain_size = blowup_factor * (poly.degree() + 1);
        let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
        let commit = FriRound::<D, F>::codeword_commit(&poly, domain);

        Self {
            round,
            alpha,
            commit,
            poly,
            domain,
        }
    }

    fn split_and_fold(poly: DensePolynomial<F>, alpha: usize) -> DensePolynomial<F> {
        let (mut even, mut odd) = (Vec::<F>::new(), Vec::<F>::new());
        for (i, element) in poly.coeffs().iter().enumerate() {
            if i % 2 == 0 {
                even.push(*element);
            } else {
                odd.push(*element);
            }
        }

        let (even_poly, odd_poly) = (
            DensePolynomial::from_coefficients_vec(even),
            DensePolynomial::from_coefficients_vec(odd),
        );

        even_poly
            + odd_poly.naive_mul(&DensePolynomial::from_coefficients_slice(&[F::from(
                alpha as u64,
            )]))
    }

    fn codeword_commit(
        poly: &DensePolynomial<F>,
        domain: Radix2EvaluationDomain<F>,
    ) -> MerkleTree<D, F> {
        let leafs = poly.evaluate_over_domain_by_ref(domain);
        MerkleTree::<D, _>::generate_tree(&leafs.evals)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::Goldilocks;
    use ark_ff::Field;
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
    fn test_commit_phase() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_slice(&coeffs);
        let mut fri = Fri::<Sha256, _>::new(poly, blowup_factor);

        let alpha = 1usize;
        let commitment = fri.commit_phase(alpha);
        assert_eq!(fri.rounds, 2);
        assert_eq!(commitment.len(), 2);
    }

    #[test]
    fn test_query_phase() {
        let domain_size = 32usize;
        let folded_size = domain_size >> 1;
        assert_eq!(folded_size, 16);

        let domain = Radix2EvaluationDomain::<Goldilocks>::new(domain_size).unwrap();
        let folded_domain = Radix2EvaluationDomain::<Goldilocks>::new(folded_size).unwrap();
        assert_eq!(domain.element(2), folded_domain.element(1));
    }

    #[test]
    #[ignore]
    fn test_fri_round() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_slice(&coeffs);
        // let mut fri = FriRound::<Sha256, _>::new(poly, blowup_factor);

        let alpha = 1usize;
        let beta = 1usize;
        // fri.split_and_fold(alpha);
        // let (merkle_proofs, consistency_proofs) = fri.prove(beta).unwrap();
        // assert!(fri.verify(merkle_proofs, consistency_proofs));
    }

    #[test]
    #[ignore]
    fn test_quotienting() {
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_slice(&coeffs);
        let domain = Radix2EvaluationDomain::<Goldilocks>::new(coeffs.len()).unwrap();

        let x = Goldilocks::from(57);
        let y = poly.evaluate(&x);
        let vanishing_poly = DensePolynomial::from_coefficients_slice(&[-x, Goldilocks::ONE]);
        let f_minus_y = poly.clone() - DensePolynomial::from_coefficients_slice(&[y]);
        let w = f_minus_y.clone() / vanishing_poly.clone();

        assert_eq!(w * vanishing_poly, f_minus_y);
    }
}
