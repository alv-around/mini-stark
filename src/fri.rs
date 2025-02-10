use super::merkle::{MerklePath, MerkleRoot, MerkleTree};
use super::util::is_power_of_two;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::Digest;
use std::iter::zip;

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
    points: Vec<[(F, F); 3]>,
    queries: Vec<[MerklePath<D, F>; 3]>,
    quotients: Vec<DensePolynomial<F>>,
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
        self.beta = Some(beta % self.domain.size());
        let alpha = self.alpha.unwrap();
        let mut queries = Vec::new();
        let mut points = Vec::new();
        let mut quotients = Vec::new();

        let mut previous_poly = &self.poly;
        let mut previous_commit = &self.commit;
        let mut previous_domain = &self.domain;
        let mut previous_beta = self.beta.unwrap();
        for round in self.round_state.iter() {
            assert_eq!(previous_domain.size() >> 1, round.domain.size());

            let x1 = previous_domain.element(previous_beta); // f_{u-1}(z)
            let x2 = previous_domain.element(round.domain.size() + previous_beta); // f_{u-i}(-z)
            let x3 = round.domain.element(previous_beta); // f_u(z^2)
            assert_eq!(x2, x1.neg());
            assert_eq!(x3, previous_domain.element(2 * previous_beta));

            let y1 = previous_poly.evaluate(&x1);
            let y2 = previous_poly.evaluate(&x2);
            let y3 = round.poly.evaluate(&x3);

            points.push([(x1, y1), (x2, y2), (x3, y3)]);

            // quotienting
            // TODO: find a less manual way in ark-poly
            let a = (y2 - y1) / (x2 - x1);
            let b = y1 - a * x1;
            let g = DensePolynomial::from_coefficients_vec(vec![a, b]);
            let numerator = previous_poly.clone() - g;
            let vanishing_poly = DensePolynomial::from_coefficients_vec(vec![-x1, F::ONE])
                * DensePolynomial::from_coefficients_vec(vec![-x2, F::ONE]);
            let q = numerator / vanishing_poly;
            println!("quotient: {:?}", q);
            quotients.push(q);

            // merkle commits
            let proof1 = previous_commit.generate_proof(y1).unwrap();
            let proof2 = previous_commit.generate_proof(y2).unwrap();
            let proof3 = round.commit.generate_proof(y3).unwrap();
            queries.push([proof1, proof2, proof3]);

            previous_poly = &round.poly;
            previous_commit = &round.commit;
            previous_domain = &round.domain;
            previous_beta %= previous_domain.size();
            println!("another round achieved");
        }

        Ok(FriProof {
            points,
            queries,
            quotients,
        })
    }

    pub fn verify(&self, commitments: Vec<MerkleRoot<D>>, proof: FriProof<D, F>) -> bool {
        if self.beta.is_none() || self.alpha.is_none() {
            return false;
        }
        assert_eq!(commitments.len(), proof.points.len());

        let mut beta = self.beta.unwrap();
        let alpha = F::from(self.alpha.unwrap() as u64);

        let mut index = 0usize;
        for ([(x1, y1), (x2, y2), (x3, y3)], [path1, path2, path3]) in
            zip(proof.points, proof.queries)
        {
            let quotient = proof.quotients[index].clone();
            let vanishing_poly = DensePolynomial::from_coefficients_vec(vec![-x1, F::ONE])
                * DensePolynomial::from_coefficients_vec(vec![-x2, F::ONE]);
            let poly = quotient * vanishing_poly;
            // TODO: check for degree of polynomial
            assert_eq!(poly.evaluate(&x1), F::ZERO);
            assert_eq!(poly.evaluate(&x2), F::ZERO);

            // linearity_check
            let two_inv = F::from(2).inverse().unwrap();
            let x1_inv = x1.inverse().unwrap();
            let even_part = (y1 + y2) * two_inv;
            let odd_part = (y1 - y2) * two_inv * x1_inv;
            assert_eq!(y3, even_part + alpha * odd_part);

            // opening checking
            commitments[index].check_proof(path1);
            commitments[index].check_proof(path2);
            commitments[index].check_proof(path3);

            index += 1;
        }

        true
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
    use ark_ff::{AdditiveGroup, UniformRand};
    use ark_std::test_rng;
    use sha2::Sha256;

    #[test]
    fn test_fri_new() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let mut fri = Fri::<Sha256, _>::new(poly, blowup_factor);
        assert_eq!(fri.rounds, 2);

        let alpha = 3usize;
        let commitment = fri.commit_phase(alpha);
        assert_eq!(commitment.len(), 2);

        let beta = 2usize;
        let proof_result = fri.query_phase(beta);
        assert!(proof_result.is_ok());
        let proof = proof_result.unwrap();
        assert!(fri.verify(commitment, proof));
    }
}
