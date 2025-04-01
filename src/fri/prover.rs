use super::FriProof;
use crate::merkle::{MerkleTree, Tree};
use crate::Hash;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use ark_poly::EvaluationDomain;
use ark_poly::{DenseUVPolynomial, Polynomial};
use digest::core_api::BlockSizeUser;
use digest::{Digest, FixedOutputReset};
use nimue::plugins::ark::FieldChallenges;
use nimue::DigestBridge;
use nimue::{ByteChallenges, ByteWriter, Merlin};

pub struct FriProver<'a, const TREE_WIDH: usize, D: Digest, F: PrimeField>
where
    F: PrimeField,
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
{
    round_num: usize,
    transcript: &'a mut Merlin<DigestBridge<D>>,
    rounds: Vec<FriRound<TREE_WIDH, D, F>>,
}

impl<'a, const W: usize, D, F> FriProver<'a, W, D, F>
where
    D: Digest + FixedOutputReset + BlockSizeUser + Clone,
    F: PrimeField,
{
    pub fn new(
        transcript: &'a mut Merlin<DigestBridge<D>>,
        poly: DensePolynomial<F>,
        blowup_factor: usize,
    ) -> Self {
        let d = poly.degree();
        let domain = Radix2EvaluationDomain::<F>::new(d * blowup_factor).unwrap();
        let domain_size = domain.size as usize;
        let round_num = domain.log_size_of_group as usize;

        // degree padding
        let power_offset = domain_size - d - 1;
        let x_power = SparsePolynomial::<F>::from_coefficients_vec(vec![(power_offset, F::ONE)]);
        let poly_offset = DensePolynomial::from(x_power).naive_mul(&poly) + poly;

        let mut rounds = Vec::<FriRound<W, D, F>>::with_capacity(round_num);
        let first_round = FriRound::new(poly_offset, domain_size);
        rounds.push(first_round);

        Self {
            round_num,
            rounds,
            transcript,
        }
    }

    pub fn get_initial_commit(&self) -> Hash<D> {
        self.rounds[0].commit.root()
    }

    pub fn prove(mut self) -> (FriProof<D, F>, Vec<u8>) {
        self.commit_phase();
        self.query_phase().unwrap()
    }

    pub fn commit_phase(&mut self) -> Vec<Hash<D>> {
        assert_eq!(self.rounds.len(), 1);

        let mut commits = Vec::new();
        for i in 1..self.round_num {
            let previous_round = &self.rounds[i - 1];
            let commit = previous_round.commit.root();
            self.transcript.add_bytes(&commit).unwrap();
            commits.push(commit);

            let alpha: [F; 1] = self.transcript.challenge_scalars().unwrap();
            let folded_poly =
                FriRound::<W, D, _>::split_and_fold(&previous_round.poly.clone(), alpha[0]);
            let domain_size = folded_poly.degree() + 1;

            let round = FriRound::<W, D, _>::new(folded_poly, domain_size);
            self.rounds.push(round);
        }

        let previous_round = &self.rounds.last().unwrap();
        let commit = previous_round.commit.root();
        self.transcript.add_bytes(&commit).unwrap();
        commits.push(commit);

        commits
    }

    pub fn query_phase(&mut self) -> Result<(FriProof<D, F>, Vec<u8>), &str> {
        let mut beta = [0u8; 8]; // usize is 64-bits
        self.transcript.fill_challenge_bytes(&mut beta).unwrap();
        let mut beta = usize::from_le_bytes(beta);
        let domain_size = self.rounds[0].domain.size();
        if beta > domain_size {
            beta %= domain_size;
        }

        let mut queries = Vec::new();
        let mut points = Vec::new();
        let mut quotients = Vec::new();

        let mut rounds_iter = self.rounds.iter_mut();
        let previous_round = rounds_iter.next().unwrap();
        let mut previous_poly = &previous_round.poly;
        let mut previous_commit = &previous_round.commit;
        let mut previous_domain = &previous_round.domain;
        for round in rounds_iter {
            assert_eq!(previous_domain.size() / W, round.domain.size());

            let x1 = previous_domain.element(beta);
            let x2 = previous_domain.element(round.domain.size() + beta);
            let x3 = round.domain.element(beta);
            let y1 = previous_poly.evaluate(&x1);
            let y2 = previous_poly.evaluate(&x2);
            let y3 = round.poly.evaluate(&x3);
            points.push([(x1, y1), (x2, y2), (x3, y3)]);
            assert_eq!(x3, previous_domain.element(2 * beta));

            // quotienting
            // g(x) = ax + b
            let a = (y2 - y1) / (x2 - x1);
            let b = y1 - a * x1;
            let g = DensePolynomial::from_coefficients_vec(vec![b, a]);

            // q(x) = f(x) - g(x) / Z(x)
            let numerator = previous_poly.clone() - g;
            let vanishing_poly = FriProver::<W, D, F>::calculate_vanishing_poly(&[x1, x2]);
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
            beta %= round.domain.size();
            println!("another round achieved");
        }

        Ok((
            FriProof {
                points,
                queries,
                quotients,
            },
            self.transcript.transcript().to_vec(),
        ))
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
    use crate::fri::fiatshamir::FriIOPattern;
    use nimue::IOPattern;
    use sha2::{self, Sha256};

    const TWO: usize = 2;

    #[test]
    fn test_fri_prover_new() {
        let blowup_factor = 2usize;
        let coeffs = (0..4).map(Goldilocks::from).collect::<Vec<_>>();
        let poly = DensePolynomial::from_coefficients_vec(coeffs);
        let io: IOPattern<DigestBridge<Sha256>> = FriIOPattern::<_, Goldilocks>::new_fri("🍟", 3);
        let mut transcript = io.to_merlin();

        let fri = FriProver::<TWO, Sha256, _>::new(&mut transcript, poly, blowup_factor);
        assert_eq!(fri.round_num, 3);

        let _proof = fri.prove();
    }
}
