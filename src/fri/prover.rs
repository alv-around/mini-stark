use super::FriProof;
use crate::merkle::{MerkleRoot, MerkleTree, Tree};
use crate::util::is_power_of_two;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::Digest;

pub struct Fri<const TREE_WIDH: usize, D: Digest, F: PrimeField> {
    rounds: usize,
    blowup_factor: usize,
    poly: DensePolynomial<F>,
    domain: Radix2EvaluationDomain<F>,
    commit: MerkleTree<TREE_WIDH, D, F>,
    round_state: Vec<FriRound<TREE_WIDH, D, F>>,
    alpha: Option<usize>,
    beta: Option<usize>,
}

impl<const W: usize, D: Digest, F: PrimeField> Fri<W, D, F> {
    pub fn new(poly: DensePolynomial<F>, blowup_factor: usize) -> Self {
        let d = poly.degree();
        let domain_size = (d + 1) * blowup_factor;
        if !is_power_of_two(domain_size) {
            panic!("blowup factor and degree of polynomial must be a power of 2");
        }

        let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
        let commit = FriRound::<W, D, F>::codeword_commit(&poly, domain.clone());

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

    pub fn generate_commit(&self) -> MerkleRoot<D> {
        MerkleRoot(self.commit.root())
    }

    pub fn commit_phase(&mut self, alpha: usize) -> Vec<MerkleRoot<D>> {
        self.alpha = Some(alpha);
        let mut oracles = Vec::new();
        let mut previous_poly = self.poly.clone();
        let mut previous_alpha = self.alpha.unwrap();
        for i in 0..self.rounds {
            let round =
                FriRound::<W, D, _>::new(previous_poly, previous_alpha, self.blowup_factor, i);
            previous_poly = round.poly.clone();
            previous_alpha = round.alpha;
            oracles.push(MerkleRoot(round.commit.root()));
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
        let mut points = Vec::new();
        let mut quotients = Vec::new();

        let mut previous_poly = &self.poly;
        let mut previous_commit = &self.commit;
        let mut previous_domain = &self.domain;
        let mut previous_beta = beta;
        for round in self.round_state.iter() {
            assert_eq!(previous_domain.size() >> 1, round.domain.size());

            let x1 = previous_domain.element(previous_beta);
            let x2 = previous_domain.element(round.domain.size() + previous_beta);
            let x3 = round.domain.element(previous_beta);

            let y1 = previous_poly.evaluate(&x1);
            let y2 = previous_poly.evaluate(&x2);
            let y3 = round.poly.evaluate(&x3);

            // FIXME: solve linearity test
            // let line = round.interpolate(&[y1, y2, y3]);
            // assert_eq!(line.len(), 2);

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
            quotients.push(q.to_vec());

            // merkle commits
            let proof1 = previous_commit.generate_proof(&y1).unwrap();
            let proof2 = previous_commit.generate_proof(&y2).unwrap();
            let proof3 = round.commit.generate_proof(&y3).unwrap();
            queries.push([proof1, proof2, proof3]);

            previous_poly = &round.poly;
            previous_commit = &round.commit;
            previous_domain = &round.domain;
            previous_beta *= previous_beta;
            println!("another round achieved");
        }

        Ok(FriProof {
            points,
            queries,
            quotients,
        })
    }
}

struct FriRound<const TREE_WIDH: usize, D: Digest, F: PrimeField> {
    round: usize,
    alpha: usize,
    poly: DensePolynomial<F>,
    commit: MerkleTree<TREE_WIDH, D, F>,
    domain: Radix2EvaluationDomain<F>,
}

impl<const W: usize, D: Digest, F: PrimeField> FriRound<W, D, F> {
    fn new(
        previous_poly: DensePolynomial<F>,
        alpha: usize,
        blowup_factor: usize,
        round: usize,
    ) -> Self {
        let poly = FriRound::<W, D, F>::split_and_fold(previous_poly, alpha);
        let domain_size = blowup_factor * (poly.degree() + 1);
        let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
        let commit = FriRound::<W, D, F>::codeword_commit(&poly, domain);

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

    fn interpolate(&self, evals: &[F]) -> Vec<F> {
        self.domain.ifft(evals)
    }

    fn codeword_commit(
        poly: &DensePolynomial<F>,
        domain: Radix2EvaluationDomain<F>,
    ) -> MerkleTree<W, D, F> {
        let leafs = poly.evaluate_over_domain_by_ref(domain);
        MerkleTree::<W, D, _>::new(&leafs.evals)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::Goldilocks;
    use ark_ff::{AdditiveGroup, UniformRand};
    use ark_std::test_rng;
    use sha2::Sha256;

    const TWO: usize = 2;

    #[test]
    fn test_fri_new() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let mut fri = Fri::<TWO, Sha256, _>::new(poly, blowup_factor);
        assert_eq!(fri.rounds, 2);

        let alpha = 3usize;
        let commitment = fri.commit_phase(alpha);
        assert_eq!(commitment.len(), 2);

        let beta = 2usize;
        let proof_result = fri.query_phase(beta);
        assert!(proof_result.is_ok());
    }

    #[test]
    fn test_linearity_test() {
        let coeffs = vec![
            Goldilocks::rand(&mut test_rng()),
            Goldilocks::rand(&mut test_rng()),
            Goldilocks::ZERO,
            Goldilocks::ZERO,
        ];
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let points = (1..7)
            .map(Goldilocks::from)
            .map(|i| poly.evaluate(&i))
            .collect::<Vec<Goldilocks>>();

        let round = FriRound::<TWO, Sha256, _>::new(poly, 1, 2, 0);
        let coeffs = round.interpolate(&points);
        assert!(coeffs.len() == 2);
    }
}
