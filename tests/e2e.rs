use ark_ff::{AdditiveGroup, Field};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial};
use mini_starks::air::{Provable, TraceTable};
use mini_starks::field::Goldilocks;
use mini_starks::starks::{Stark, StarkConfig};
use sha2::Sha256;

const ONE: Goldilocks = Goldilocks::ONE;
const ZERO: Goldilocks = Goldilocks::ZERO;

struct FibonacciClaim {
    step: usize, // nth fibonacci number
    output: Goldilocks,
}

struct Witness {
    secret_b: Goldilocks,
}

impl Provable<Witness, Goldilocks> for FibonacciClaim {
    fn trace(&self, witness: &Witness) -> TraceTable<Goldilocks> {
        let trace_width = 3usize;
        let mut trace = TraceTable::new(self.step, trace_width);

        // initial state
        let mut a = ONE;
        let mut b = witness.secret_b;
        let mut c = a + b;

        // set initial state constrains
        trace.add_boundary_constrain(0, 0);
        trace.add_boundary_constrain(0, 1);
        trace.add_boundary_constrain(0, 2);

        // trace
        for i in 0..trace.step_number() {
            trace.add_row(i, vec![a, b, c]);
            a = b;
            b = c;
            c = a + b;
        }

        // set output constrains
        trace.add_boundary_constrain(self.step - 1, 2);

        // add transition constrains
        // a[1] == b[0]
        trace.add_transition_constrain(Box::new(move |trace_polys| {
            trace_polys[0].clone() * DensePolynomial::from_coefficients_vec(vec![trace.omega])
                - trace_polys[1].clone()
        }));
        // b[1] == c[0]
        trace.add_transition_constrain(Box::new(move |trace_polys| {
            trace_polys[0].clone() * DensePolynomial::from_coefficients_vec(vec![trace.omega])
                - trace_polys[1].clone()
        }));
        trace.add_transition_constrain(Box::new(move |trace_polys| {
            trace_polys[2].clone() - trace_polys[0].clone() - trace_polys[1].clone()
        }));

        trace
    }
}

fn test_setup() -> (Witness, FibonacciClaim) {
    let witness = Witness {
        secret_b: Goldilocks::from(2),
    };
    let claim = FibonacciClaim {
        step: 9,
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

    let blowup_factor = 2;
    let columns = trace.constrain_number();

    let config =
        StarkConfig::<Sha256, Goldilocks>::new(20, blowup_factor, trace.step_number(), columns);
    let proof_system = Stark::new(config);
    let proof = proof_system.prove(claim, witness).unwrap();

    let is_alright = proof_system.verify(constrains, proof);
    assert!(is_alright);
}
