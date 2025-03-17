use ark_ff::{AdditiveGroup, Field};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial};
use mini_starks::air::*;
use mini_starks::field::Goldilocks;

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
    fn trace(&self, witness: Witness) -> TraceTable<Goldilocks> {
        let trace_width = 2usize;
        let mut trace = TraceTable::new(self.step, trace_width);

        // initial state
        let mut a = ONE;
        let mut b = witness.secret_b;
        // trace
        for _ in 0..trace.len() {
            let c = a + b;
            trace.add_row(vec![a, b]);
            a = b;
            b = c;
        }

        trace
    }
}

impl Verifiable<Goldilocks> for FibonacciClaim {
    fn derive_constrains(&self, trace: &TraceTable<Goldilocks>) -> Constrains<Goldilocks> {
        let mut constrains = trace.interpolate_col_polys();

        // boundary polynomials
        constrains.add_boundary_constrain(0, 0, Goldilocks::ONE);
        constrains.add_boundary_constrain(self.step - 1, 1, self.output);

        // transition polynomials
        let col_a = constrains.get_trace_poly(0);
        let col_b = constrains.get_trace_poly(1);
        let omega = DensePolynomial::from_coefficients_slice(&[constrains.get_domain().element(1)]);

        let carry_over_poly = col_a.clone() * omega.clone() - col_b.clone();
        let sum_poly = col_b.clone() * omega - (col_a + col_b);
        constrains.add_transition_constrain(carry_over_poly);
        constrains.add_transition_constrain(sum_poly);

        constrains
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
    let trace = claim.trace(witness);
    let constrains = claim.derive_constrains(&trace);
    let domain = constrains.get_domain();

    // check output constrain
    let omega_4 = domain.element(claim.step - 1);
    let z = DensePolynomial::from_coefficients_slice(&[-omega_4, ONE]);
    let boundary = constrains.get_boundary_constrain(1).clone() * z;
    assert_eq!(boundary.evaluate(&omega_4), ZERO);

    let carry_over_constrain = constrains
        .get_transition_constrain(0)
        .mul_by_vanishing_poly(domain);
    let sum_constrain = constrains
        .get_transition_constrain(1)
        .mul_by_vanishing_poly(domain);
    for i in 0..trace.len() - 1 {
        let w_i = domain.element(i);
        assert_eq!(carry_over_constrain.evaluate(&w_i), ZERO);
        assert_eq!(sum_constrain.evaluate(&w_i), ZERO);
    }
}
