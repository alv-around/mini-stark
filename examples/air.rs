use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use mini_starks::air::{Constrains, Provable, TraceTable, Verifiable};
use mini_starks::field::Goldilocks;
use std::error::Error;

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
        let mut a = Goldilocks::ONE;
        let mut b = witness.secret_b;
        // trace
        for _ in 1..trace.len() {
            let c = a + b;
            trace.add_row(vec![a, b]);
            a = b;
            b = c;
        }

        trace
    }
}

impl Verifiable<Goldilocks> for FibonacciClaim {
    fn constrains(&self, trace: &TraceTable<Goldilocks>) -> Constrains<Goldilocks> {
        let mut constrains = Constrains::new(trace);

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

fn main() -> Result<(), Box<dyn Error>> {
    Ok(())
}
