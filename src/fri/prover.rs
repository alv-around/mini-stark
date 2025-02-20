use super::FriProof;
use crate::merkle::{MerkleRoot, MerkleTree, Tree};
use crate::util::logarithm_of_two_k;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::EvaluationDomain;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::Digest;

pub struct Fri<const TREE_WIDH: usize, D: Digest, F: PrimeField> {
    rounds: usize,
    blowup_factor: usize,
    beta: Option<usize>,
    commits: Vec<FriRound<TREE_WIDH, D, F>>,
}

impl<const W: usize, D: Digest, F: PrimeField> Fri<W, D, F> {
    pub fn new(poly: DensePolynomial<F>, blowup_factor: usize) -> Self {
        let d = poly.degree();
        let domain_size = (d + 1) * blowup_factor;
        let rounds = logarithm_of_two_k::<W>(domain_size);
        if rounds.is_err() {
            panic!("blowup factor and degree of polynomial must be a power of 2");
        }
        let rounds = rounds.unwrap();
        let mut commits = Vec::<FriRound<W, D, F>>::with_capacity(rounds);
        let first_round = FriRound::new(poly, domain_size);
        commits.push(first_round);

        Self {
            rounds,
            blowup_factor,
            beta: None,
            commits,
        }
    }

    pub fn generate_commit(&self) -> MerkleRoot<D> {
        MerkleRoot(self.commits[0].commit.root())
    }

    fn domain_size(&self, poly: &DensePolynomial<F>) -> usize {
        self.blowup_factor * (poly.degree() + 1)
    }

    pub fn commit_phase(&mut self, alphas: &[F]) -> Vec<MerkleRoot<D>> {
        assert_eq!(alphas.len(), self.rounds - 1);
        let mut oracles = Vec::new();

        for (i, _) in (1..self.rounds).enumerate() {
            let previous_round = &self.commits[i];
            let previous_poly = previous_round.poly.clone();
            // FIXME: alpha should be different in each round
            let alpha = alphas[i];
            let folded_poly = FriRound::<W, D, _>::split_and_fold(&previous_poly, alpha);
            let domain_size = self.domain_size(&folded_poly);

            let round = FriRound::<W, D, _>::new(folded_poly, domain_size);
            oracles.push(MerkleRoot(round.commit.root()));
            self.commits.push(round);
        }

        oracles
    }

    pub fn query_phase(&mut self, beta: usize) -> Result<FriProof<D, F>, &str> {
        if self.commits.len() < self.rounds || self.beta.is_some() {
            return Err("wrong time");
        }
        self.beta = Some(beta);
        let mut queries = Vec::new();
        let mut points = Vec::new();
        let mut quotients = Vec::new();

        let mut rounds_iter = self.commits.iter_mut();
        let previous_round = rounds_iter.next().unwrap();
        let mut previous_poly = &previous_round.poly;
        let mut previous_commit = &previous_round.commit;
        let mut previous_domain = &previous_round.domain;
        let mut previous_beta = beta;
        for round in rounds_iter {
            assert_eq!(previous_domain.size() / W, round.domain.size());

            let x1 = previous_domain.element(previous_beta);
            let x2 = previous_domain.element(round.domain.size() + previous_beta);
            let x3 = round.domain.element(previous_beta);
            let y1 = previous_poly.evaluate(&x1);
            let y2 = previous_poly.evaluate(&x2);
            let y3 = round.poly.evaluate(&x3);
            points.push([(x1, y1), (x2, y2), (x3, y3)]);
            assert_eq!(x3, previous_domain.element(2 * previous_beta));

            // quotienting
            // g(x) = ax + b
            let a = (y2 - y1) / (x2 - x1);
            let b = y1 - a * x1;
            let g = DensePolynomial::from_coefficients_vec(vec![b, a]);

            // q(x) = f(x) - g(x) / Z(x)
            let numerator = previous_poly.clone() - g;
            let vanishing_poly = Fri::<W, D, F>::calculate_vanishing_poly(&[x1, x2]);
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
            previous_beta %= round.domain.size();
            println!("another round achieved");
        }

        Ok(FriProof {
            points,
            queries,
            quotients,
        })
    }

    fn calculate_vanishing_poly(roots: &[F]) -> DensePolynomial<F> {
        roots
            .iter()
            .map(|i| DensePolynomial::from_coefficients_slice(&[-*i, F::ONE]))
            .reduce(|acc, e| acc * e)
            .unwrap()
    }
}

struct FriRound<const TREE_WIDH: usize, D: Digest, F: PrimeField> {
    poly: DensePolynomial<F>,
    commit: MerkleTree<TREE_WIDH, D, F>,
    domain: Radix2EvaluationDomain<F>,
}

impl<const W: usize, D: Digest, F: PrimeField> FriRound<W, D, F> {
    fn new(poly: DensePolynomial<F>, domain_size: usize) -> Self {
        let domain = Radix2EvaluationDomain::<F>::new(domain_size).unwrap();
        let commit = FriRound::<W, D, F>::codeword_commit(&poly, domain);

        Self {
            commit,
            poly,
            domain,
        }
    }

    fn split_and_fold(poly: &DensePolynomial<F>, alpha: F) -> DensePolynomial<F> {
        let mut coeffs_vectors = Vec::<Vec<F>>::new();
        for _ in 0..W {
            coeffs_vectors.push(Vec::<F>::new());
        }
        for (i, element) in poly.coeffs().iter().enumerate() {
            let index = i % W;
            coeffs_vectors[index].push(*element);
        }

        coeffs_vectors
            .iter()
            .map(|coeffs: &Vec<F>| DensePolynomial::from_coefficients_slice(coeffs))
            .enumerate()
            .map(|(i, poly)| {
                poly.naive_mul(&DensePolynomial::from_coefficients_slice(&[
                    alpha.pow([i as u64])
                ]))
            })
            .reduce(|acc, e| acc + e)
            .unwrap()
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
    use ark_ff::UniformRand;
    use ark_std::test_rng;
    use sha2::Sha256;

    const TWO: usize = 2;

    #[test]
    fn test_fri_new() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let mut fri = Fri::<TWO, Sha256, _>::new(poly, blowup_factor);
        assert_eq!(fri.rounds, 3);

        let mut alphas = Vec::new();
        for _ in 1..fri.rounds {
            alphas.push(Goldilocks::rand(&mut test_rng()));
        }
        let commitment = fri.commit_phase(&alphas);
        assert_eq!(commitment.len(), 2);

        let beta = 2usize;
        let proof_result = fri.query_phase(beta);
        assert!(proof_result.is_ok());
    }
}
