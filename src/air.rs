use crate::util::ceil_log2_k;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;

const TWO: usize = 2;

pub trait Provable<T, F: PrimeField> {
    fn trace(&self, witness: T) -> TraceTable<F>;
}

pub trait Verifiable<F: PrimeField> {
    fn constrains(&self, trace: &TraceTable<F>) -> Constrains<F>;
}

pub struct TraceTable<F: PrimeField> {
    depth: usize,
    width: usize,
    data: Vec<F>,
}

impl<F: PrimeField> TraceTable<F> {
    pub fn new(steps: usize, width: usize) -> Self {
        let depth = 1usize << ceil_log2_k::<TWO>(steps);
        let data = Vec::<F>::with_capacity(depth * width);
        Self { depth, width, data }
    }

    pub fn add_row(&mut self, row: Vec<F>) {
        assert_eq!(row.len(), self.width);
        self.data.extend(row);
    }

    pub fn get_value(&self, row: usize, col: usize) -> &F {
        assert!(row < self.len() && col < self.width());
        &self.data[row * self.width + col]
    }

    pub fn len(&self) -> usize {
        self.depth
    }

    pub fn is_empty(&self) -> bool {
        self.depth == 0 || self.width == 0
    }

    pub fn width(&self) -> usize {
        self.width
    }
}

pub struct Constrains<F: PrimeField> {
    domain: Radix2EvaluationDomain<F>,
    trace_polys: Vec<DensePolynomial<F>>,
    constrains: Vec<DensePolynomial<F>>,
}

impl<F: PrimeField> Constrains<F> {
    pub fn new(trace: &TraceTable<F>) -> Self {
        let domain = Radix2EvaluationDomain::new(trace.len()).unwrap();
        let trace_polys = Self::interpolate_col_polys(trace);
        let constrains = Vec::<DensePolynomial<F>>::new();
        Self {
            domain,
            trace_polys,
            constrains,
        }
    }

    pub fn get_trace_poly(&self, col: usize) -> DensePolynomial<F> {
        self.trace_polys[col].clone()
    }

    pub fn get_domain(&self) -> Radix2EvaluationDomain<F> {
        self.domain
    }

    pub fn constrains(self) -> Vec<DensePolynomial<F>> {
        self.constrains
    }

    pub fn add_boundary_constrain(&mut self, row: usize, col: usize, value: F) {
        let omega_i = self.domain.element(row);
        let root = DensePolynomial::from_coefficients_vec(vec![-omega_i, F::ONE]);
        let y = DensePolynomial::from_coefficients_slice(&[value]);
        let poly = (self.trace_polys[col].clone() - y) / root;
        self.constrains.push(poly);
    }

    pub fn add_transition_constrain(&mut self, poly: DensePolynomial<F>) {
        let domain_size = self.domain.size();
        let root = DensePolynomial::from(SparsePolynomial::from_coefficients_slice(&[
            (0, -F::ONE),
            (domain_size, F::ONE),
        ]));
        let omega_n_minus_one = self.domain.element(domain_size - 1);
        let last_row_poly =
            DensePolynomial::from_coefficients_vec(vec![-omega_n_minus_one, F::ONE]);
        self.constrains.push(poly * last_row_poly / root);
    }

    fn interpolate_col_polys(trace: &TraceTable<F>) -> Vec<DensePolynomial<F>> {
        let trace_depth = trace.len();
        let trace_width = trace.width();
        let mut trace_polys = Vec::with_capacity(trace_width);

        for i in 0..trace_width {
            let coeffs = (0..trace_depth)
                .map(|j| *trace.get_value(j, i))
                .collect::<Vec<_>>();
            let poly = DensePolynomial::from_coefficients_vec(coeffs);
            trace_polys.push(poly);
        }

        trace_polys
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::Goldilocks;
    use ark_ff::{AdditiveGroup, Field};
    use ark_poly::Polynomial;

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
        fn constrains(&self, trace: &TraceTable<Goldilocks>) -> Constrains<Goldilocks> {
            let mut constrains = Constrains::new(trace);

            // boundary polynomials
            constrains.add_boundary_constrain(0, 0, Goldilocks::ONE);
            constrains.add_boundary_constrain(self.step - 1, 1, self.output);

            // transition polynomials
            let col_a = constrains.get_trace_poly(0);
            let col_b = constrains.get_trace_poly(1);
            let omega =
                DensePolynomial::from_coefficients_slice(&[constrains.get_domain().element(1)]);

            let carry_over_poly = col_a.clone() * omega.clone() - col_b.clone();
            let sum_poly = col_b.clone() * omega - (col_a + col_b);
            constrains.add_transition_constrain(carry_over_poly);
            constrains.add_transition_constrain(sum_poly);

            constrains
        }
    }

    #[test]
    fn test_fibonacci_trace_and_constrains() {
        let witness = Witness {
            secret_b: Goldilocks::from(2),
        };
        let claim = FibonacciClaim {
            step: 5,
            output: Goldilocks::from(13),
        };

        let trace = claim.trace(witness);
        assert_eq!(trace.len(), 8);
        assert_eq!(trace.width, 2);
        assert_eq!(trace.data.capacity(), trace.len() * trace.width());
        assert_eq!(*trace.get_value(claim.step - 1, 1), claim.output);
        assert_eq!(*trace.get_value(claim.step, 0), claim.output);

        let one = Goldilocks::ONE;
        let zero = Goldilocks::ZERO;
        let constrains = claim.constrains(&trace);
        let domain = constrains.get_domain();
        let constrains = constrains.constrains();

        // check constrains
        let boundary1 =
            constrains[0].clone() * DensePolynomial::from_coefficients_slice(&[-one, one]);
        assert_eq!(boundary1.evaluate(&one), zero);

        let omega_4 = domain.element(claim.step - 1);
        let z = DensePolynomial::from_coefficients_slice(&[-omega_4, one]);
        let boundary2 = constrains[1].clone() * z;
        assert_eq!(boundary2.evaluate(&claim.output), zero);
    }
}
