use crate::util::is_power_of_two;
use ark_ff::PrimeField;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_std::test_rng;

pub trait Provable<W, F: PrimeField> {
    fn trace(&self, witness: &W) -> TraceTable<F>;
    // fn derive_constrians(&self) -> Constrains<F>;
}

// Matrix is the basic struct for the trace, and ldes
pub(crate) struct Matrix<F: PrimeField> {
    length: usize,
    width: usize,
    data: Vec<F>,
}

impl<F: PrimeField> Matrix<F> {
    pub(crate) fn new(length: usize, width: usize, entries: Option<Vec<F>>) -> Self {
        assert!(is_power_of_two(length * width));
        let data = match entries {
            Some(data) => {
                assert_eq!(data.len(), length * width);
                data
            }
            None => vec![F::ZERO; length * width],
        };
        Self {
            length,
            width,
            data,
        }
    }

    pub(crate) fn get_data(&self) -> &[F] {
        &self.data[..]
    }

    pub(crate) fn get_value(&self, row: usize, col: usize) -> &F {
        assert!(row < self.length && col < self.width);
        &self.data[row * self.width + col]
    }

    pub(crate) fn len(&self) -> usize {
        self.length
    }

    #[allow(dead_code)]
    pub(crate) fn is_empty(&self) -> bool {
        self.length == 0 || self.width == 0
    }

    pub(crate) fn add_col(&mut self, index: usize, col: Vec<F>) {
        assert_eq!(col.len(), self.length);
        assert!(index < self.width);
        for (i, val) in col.into_iter().enumerate() {
            self.data[i * self.width + index] = val;
        }
    }
}

pub struct TraceTable<F: PrimeField> {
    pub(super) trace: Matrix<F>,
    steps: usize,
    domain: Radix2EvaluationDomain<F>,
    pub omega: F,
    boundaries: Vec<(usize, usize)>,
    transition_constrains: Vec<Box<dyn Fn(&Vec<DensePolynomial<F>>) -> DensePolynomial<F>>>,
}

impl<F: PrimeField> TraceTable<F> {
    pub fn new(steps: usize, registers: usize) -> Self {
        let domain = Radix2EvaluationDomain::<F>::new(steps + 1).unwrap();
        let omega = domain.group_gen();

        let mut uninitialized_trace = vec![F::ZERO; steps * registers];
        // add random padding to make zk
        let padding_length = (domain.size() - steps) * registers;
        let random_padding = (0..padding_length)
            .map(|_| F::rand(&mut test_rng()))
            .collect::<Vec<F>>();
        uninitialized_trace.extend(random_padding);

        let trace = Matrix::<F>::new(domain.size(), registers, Some(uninitialized_trace));

        let boundaries = Vec::new();
        Self {
            steps,
            domain,
            omega,
            boundaries,
            transition_constrains: Vec::new(),
            trace,
        }
    }

    pub fn step_number(&self) -> usize {
        self.steps
    }

    pub fn get_domain(&self) -> Radix2EvaluationDomain<F> {
        self.domain
    }

    pub fn add_row(&mut self, index: usize, row: Vec<F>) {
        assert_eq!(row.len(), self.trace.width);
        assert!(index < self.steps);
        for (j, val) in row.into_iter().enumerate() {
            self.trace.data[index * self.trace.width + j] = val;
        }
    }

    pub fn add_boundary_constrain(&mut self, row: usize, col: usize) {
        assert!(row < self.steps && col < self.trace.width);
        self.boundaries.push((row, col));
    }

    pub fn add_transition_constrain(
        &mut self,
        f: Box<dyn Fn(&Vec<DensePolynomial<F>>) -> DensePolynomial<F>>,
    ) {
        self.transition_constrains.push(f);
    }

    pub fn derive_constrains(&self) -> Constrains<F> {
        let mut constrains = self.get_trace_polys();

        let transition_evals = self
            .transition_constrains
            .iter()
            .map(|f| f(&constrains))
            .collect::<Vec<_>>();

        let trace_constrains_num = self.trace.width;
        let transition_constrains_num = transition_evals.len();
        constrains.extend(transition_evals);
        Constrains {
            trace_constrains_num,
            transition_constrains_num,
            constrains,
        }
    }

    // compute trace polynomials
    pub(crate) fn get_trace_polys(&self) -> Vec<DensePolynomial<F>> {
        let mut trace_polys = Vec::with_capacity(self.trace.width);
        for i in 0..self.trace.width {
            // derive trace polys
            let evaluations = (0..self.trace.length)
                .map(|j| *self.trace.get_value(j, i))
                .collect::<Vec<_>>();
            let coeffs = self.domain.ifft(&evaluations);
            let column_poly = DensePolynomial::from_coefficients_vec(coeffs);
            trace_polys.push(column_poly);
        }

        trace_polys
    }
}

pub struct Constrains<F: PrimeField> {
    trace_constrains_num: usize,
    transition_constrains_num: usize,
    constrains: Vec<DensePolynomial<F>>,
}

impl<F: PrimeField> Constrains<F> {
    pub fn len(&self) -> usize {
        self.constrains.len()
    }

    pub fn get_constrain_poly(&self, col: usize) -> DensePolynomial<F> {
        assert!(col < self.trace_constrains_num + self.transition_constrains_num);
        self.constrains[col].clone()
    }

    pub fn get_polynomials(&self) -> Vec<DensePolynomial<F>> {
        self.constrains.clone()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::Goldilocks;
    use ark_ff::{AdditiveGroup, Field};
    use ark_poly::Polynomial;

    const ONE: Goldilocks = Goldilocks::ONE;
    const ZERO: Goldilocks = Goldilocks::ZERO;

    struct FibonacciClaim {
        step: usize, // nth fibonacci number
        output: Goldilocks,
    }

    struct Witness;

    impl Provable<Witness, Goldilocks> for FibonacciClaim {
        fn trace(&self, _witness: &Witness) -> TraceTable<Goldilocks> {
            let trace_width = 2usize;
            let mut trace = TraceTable::new(self.step, trace_width);

            // initial state
            let mut a = ONE;
            let mut b = ONE;

            // set initial state constrains
            trace.add_boundary_constrain(0, 0);
            trace.add_boundary_constrain(0, 1);

            // trace
            for i in 0..self.step {
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
                trace_polys[1].clone() * trace.omega
                    - (trace_polys[0].clone() + trace_polys[1].clone())
            }));

            trace
        }
    }

    /// Test multiple trace lengths, assert step number = output, and random padding after steps
    #[test]
    fn test_air_trace() {
        let third_fibonacci = FibonacciClaim {
            step: 3,
            output: Goldilocks::from(3),
        };
        let trace = third_fibonacci.trace(&Witness);
        assert_eq!(trace.trace.len(), 4);
        assert_eq!(trace.trace.width, 2);
        assert_eq!(*trace.trace.get_value(0, 0), ONE);
        assert_eq!(
            *trace.trace.get_value(third_fibonacci.step - 1, 1),
            third_fibonacci.output
        );
        // assert that trace has random panding from after steps
        assert_ne!(
            *trace.trace.get_value(third_fibonacci.step, 0),
            third_fibonacci.output
        );
        assert_ne!(
            *trace.trace.get_value(third_fibonacci.step, 0),
            Goldilocks::ZERO
        );

        let fourth_fib = FibonacciClaim {
            step: 4,
            output: Goldilocks::from(5),
        };
        let trace = fourth_fib.trace(&Witness);
        assert_eq!(trace.trace.length, 8);
        assert_eq!(trace.trace.width, 2);
        assert_eq!(*trace.trace.get_value(0, 1), ONE);
        assert_eq!(
            *trace.trace.get_value(fourth_fib.step - 1, 1),
            fourth_fib.output
        );
        assert_ne!(
            *trace.trace.get_value(fourth_fib.step, 0),
            fourth_fib.output
        );
        assert_ne!(*trace.trace.get_value(fourth_fib.step, 0), Goldilocks::ZERO);

        let fith_fib = FibonacciClaim {
            step: 5,
            output: Goldilocks::from(8),
        };
        let trace = fith_fib.trace(&Witness);
        assert_eq!(trace.trace.length, 8);
        assert_eq!(trace.trace.width, 2);
        assert_eq!(*trace.trace.get_value(0, 1), ONE);
        assert_eq!(
            *trace.trace.get_value(fith_fib.step - 1, 1),
            fith_fib.output
        );
        assert_ne!(*trace.trace.get_value(fith_fib.step, 0), fith_fib.output);
        assert_ne!(*trace.trace.get_value(fith_fib.step, 0), Goldilocks::ZERO);
    }

    #[test]
    fn test_air_trace_polynomials() {
        let claim = FibonacciClaim {
            step: 3,
            output: Goldilocks::from(3),
        };
        let trace = claim.trace(&Witness);
        let trace_polys = trace.get_trace_polys();
        let domain = Radix2EvaluationDomain::<Goldilocks>::new(trace.trace.length).unwrap();

        // check trace polynomials are computed correctly
        for i in 0..claim.step {
            let row = domain.element(i);
            assert_eq!(*trace.trace.get_value(i, 0), trace_polys[0].evaluate(&row));
            assert_eq!(*trace.trace.get_value(i, 1), trace_polys[1].evaluate(&row));
        }
    }

    #[test]
    fn test_air_constrains() {
        let claim = FibonacciClaim {
            step: 3,
            output: Goldilocks::from(3),
        };
        let trace = claim.trace(&Witness);
        let domain = trace.domain;
        let constrains = trace.derive_constrains();
        assert_eq!(constrains.transition_constrains_num, 2);

        // boundary_constrains
        let w0 = domain.element(0);
        let root = DensePolynomial::from_coefficients_vec(vec![-w0, ONE]);
        let boundary1 = constrains.get_constrain_poly(0);
        assert_eq!((boundary1 * root).evaluate(&ONE), ZERO);

        let w2 = domain.element(claim.step - 1);
        let root = DensePolynomial::from_coefficients_vec(vec![-w2, ONE]);
        let boundary3 = constrains.get_constrain_poly(1);
        assert_eq!((boundary3 * root).evaluate(&w2), ZERO);

        // transition_constrains
        let carry_over_constrain = constrains
            .get_constrain_poly(2)
            .mul_by_vanishing_poly(domain);
        let sum_constrain = constrains
            .get_constrain_poly(3)
            .mul_by_vanishing_poly(domain);
        for i in 0..trace.trace.length - 1 {
            let w_i = domain.element(i);
            assert_eq!(carry_over_constrain.evaluate(&w_i), ZERO);
            assert_eq!(sum_constrain.evaluate(&w_i), ZERO);
        }
    }
}
