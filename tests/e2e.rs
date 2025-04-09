use ark_ff::{AdditiveGroup, Field};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial};
use mini_starks::air::{Provable, TraceTable};
use mini_starks::fiatshamir::StarkIOPattern;
use mini_starks::field::Goldilocks;
use mini_starks::starks::Stark;
use nimue::{DigestBridge, IOPattern};
use sha2::Sha256;

const ONE: Goldilocks = Goldilocks::ONE;
const ZERO: Goldilocks = Goldilocks::ZERO;

const TWO: usize = 2usize;

struct FibonacciClaim {
    step: usize, // nth fibonacci number
    output: Goldilocks,
}

struct Witness {
    secret_b: Goldilocks,
}

impl Provable<Witness, Goldilocks> for FibonacciClaim {
    fn trace(&self, witness: &Witness) -> TraceTable<Goldilocks> {
        let trace_width = 2usize;
        let mut trace = TraceTable::new(self.step, trace_width);

        // initial state
        let mut a = ONE;
        let mut b = witness.secret_b;

        // set initial state constrains
        trace.add_boundary_constrain(0, 0);
        trace.add_boundary_constrain(0, 1);

        // trace
        for i in 0..trace.step_number() {
            let c = a + b;
            trace.add_row(i, vec![a, b]);
            a = b;
            b = c;
        }

        // set output constrains
        trace.add_boundary_constrain(self.step - 1, 1);

        // add transition constrains
        trace.add_transition_constrain(Box::new(move |trace_polys| {
            trace_polys[0].clone() * DensePolynomial::from_coefficients_vec(vec![trace.omega])
                - trace_polys[1].clone()
        }));
        trace.add_transition_constrain(Box::new(move |trace_polys| {
            trace_polys[1].clone() * trace.omega - (trace_polys[0].clone() + trace_polys[1].clone())
        }));

        trace
    }
}

fn test_setup() -> (Witness, FibonacciClaim) {
    let witness = Witness {
        secret_b: Goldilocks::from(2),
    };
    let claim = FibonacciClaim {
        step: 5,
        output: Goldilocks::from(13),
    };
    (witness, claim)
}

#[test]
fn test_fibonacci_air_constrains() {
    let (witness, claim) = test_setup();
    let trace = claim.trace(&witness);
    let constrains = trace.derive_constrains();
    let domain = trace.get_domain();

    // check output constrain
    let carry_over_constrain = constrains
        .get_constrain_poly(2)
        .mul_by_vanishing_poly(domain);
    let sum_constrain = constrains
        .get_constrain_poly(3)
        .mul_by_vanishing_poly(domain);
    for i in 0..trace.step_number() - 1 {
        let w_i = domain.element(i);
        assert_eq!(carry_over_constrain.evaluate(&w_i), ZERO);
        assert_eq!(sum_constrain.evaluate(&w_i), ZERO);
    }
}

#[test]
fn test_stark_prover() {
    let (witness, claim) = test_setup();
    let trace = claim.trace(&witness);
    let constrains = trace.derive_constrains();
    let domain = trace.get_domain();

    let io: IOPattern<DigestBridge<Sha256>> =
        StarkIOPattern::<_, Goldilocks>::new_stark(4, 80, "🐺");
    let transcript = io.to_merlin();

    let proof_system = Stark::<TWO, Sha256, Goldilocks>::new(2, 80, domain.size());
    let proof = proof_system.prove(transcript, claim, witness).unwrap();

    let is_alright = proof_system.verify(io, constrains, proof);
    assert!(is_alright);
}
