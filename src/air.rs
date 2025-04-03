use crate::util::ceil_log2_k;
use ark_ff::PrimeField;
use ark_ff::Zero;
use ark_poly::domain::Radix2EvaluationDomain;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;

const TWO: usize = 2;

pub trait Provable<W, F: PrimeField> {
    fn trace(&self, witness: &W) -> TraceTable<F>;
    // fn derive_constrians(&self) -> Constrains<F>;
}

pub struct TraceTable<F: PrimeField> {
    depth: usize,
    width: usize,
    pub omega: F,
    boundaries: Vec<(usize, usize)>,
    transition_constrains: Vec<Box<dyn Fn(&Vec<DensePolynomial<F>>) -> DensePolynomial<F>>>,
    data: Vec<F>,
}

impl<F: PrimeField> TraceTable<F> {
    pub fn new(steps: usize, width: usize) -> Self {
        // TODO: add random padding for [steps, depth) rows
        // TODO: replace ceil with domain::new
        let boundaries = Vec::new();
        let depth = 1usize << ceil_log2_k::<TWO>(steps);
        let capacity = depth * width;
        let mut data = Vec::<F>::with_capacity(capacity);
        let omega = Radix2EvaluationDomain::<F>::new(depth).unwrap().group_gen();
        data.resize(capacity, F::ZERO);
        Self {
            depth,
            width,
            omega,
            boundaries,
            transition_constrains: Vec::new(),
            data,
        }
    }

    pub fn get_data(&self) -> &[F] {
        &self.data[..]
    }

    // compute trace polynomials
    pub fn get_trace_polys(&self) -> Vec<DensePolynomial<F>> {
        let domain = Radix2EvaluationDomain::<F>::new(self.len()).unwrap();
        let trace_depth = self.len();
        let trace_width = self.width();

        let mut trace_polys = Vec::with_capacity(trace_width);
        for i in 0..trace_width {
            // derive trace polys
            let evaluations = (0..trace_depth)
                .map(|j| *self.get_value(j, i))
                .collect::<Vec<_>>();
            let coeffs = domain.ifft(&evaluations);
            let column_poly = DensePolynomial::from_coefficients_vec(coeffs);
            trace_polys.push(column_poly);
        }

        trace_polys
    }

    pub fn add_row(&mut self, index: usize, row: Vec<F>) {
        assert_eq!(row.len(), self.width);
        assert!(index < self.depth);
        for (j, val) in row.into_iter().enumerate() {
            self.data[index * self.width + j] = val;
        }
    }

    pub fn add_col(&mut self, index: usize, col: Vec<F>) {
        assert_eq!(col.len(), self.depth);
        assert!(index < self.width);
        for (i, val) in col.into_iter().enumerate() {
            self.data[i * self.width + index] = val;
        }
    }

    pub fn add_boundary_constrain(&mut self, row: usize, col: usize) {
        assert!(row < self.len() && col < self.width());
        self.boundaries.push((row, col));
    }

    pub fn add_transition_constrain(
        &mut self,
        f: Box<dyn Fn(&Vec<DensePolynomial<F>>) -> DensePolynomial<F>>,
    ) {
        self.transition_constrains.push(f);
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

    pub fn derive_constrains(&self) -> Constrains<F> {
        let domain = Radix2EvaluationDomain::<F>::new(self.len()).unwrap();
        let trace_depth = self.len();
        let trace_polys = self.get_trace_polys();

        // derive the transition quotients
        let omega_n_minus_one = domain.element(trace_depth - 1);
        let last_row_root =
            DensePolynomial::from_coefficients_vec(vec![-omega_n_minus_one, F::ONE]);
        let transition_quotients = self
            .transition_constrains
            .iter()
            .map(|f| f(&trace_polys))
            .map(|poly| {
                let (rest, quotient) = poly.divide_by_vanishing_poly(domain);
                assert_eq!(rest, DensePolynomial::zero());
                quotient * last_row_root.clone()
            })
            .collect::<Vec<_>>();

        // add boundary quotients to trace polys
        let mut trace_quotients = trace_polys
            .into_iter()
            .enumerate()
            .map(|(i, poly)| {
                let boundaries = self
                    .boundaries
                    .iter()
                    .filter(|(_, col)| *col == i)
                    .map(|(row, col)| (row, self.get_value(*row, *col)));

                let (indexes, values): (Vec<&usize>, Vec<F>) = boundaries.unzip();
                let interpolant_domain = Radix2EvaluationDomain::<F>::new(values.len()).unwrap();
                let interpolant = interpolant_domain.ifft(&values);
                let numerator = poly - DensePolynomial::from_coefficients_vec(interpolant);
                let zerofier = indexes
                    .iter()
                    .map(|i| domain.element(**i))
                    .map(|i| DensePolynomial::from_coefficients_vec(vec![-i, F::ONE]))
                    .fold(
                        DensePolynomial::from_coefficients_vec(vec![F::ONE]),
                        |acc, root| acc * root,
                    );
                numerator / zerofier
            })
            .collect::<Vec<_>>();

        let trace_constrains = self.width();
        let transition_constrains = transition_quotients.len();
        trace_quotients.extend(transition_quotients);
        Constrains {
            domain,
            trace_constrains,
            transition_constrains,
            constrains: trace_quotients,
        }
    }
}

pub struct Constrains<F: PrimeField> {
    pub domain: Radix2EvaluationDomain<F>,
    trace_constrains: usize,
    transition_constrains: usize,
    constrains: Vec<DensePolynomial<F>>,
}

impl<F: PrimeField> Constrains<F> {
    pub fn len(&self) -> usize {
        self.constrains.len()
    }

    pub fn get_constrain_poly(&self, col: usize) -> DensePolynomial<F> {
        assert!(col < self.trace_constrains + self.transition_constrains);
        self.constrains[col].clone()
    }

    pub fn get_polynomials(&self) -> Vec<DensePolynomial<F>> {
        self.constrains.clone()
    }

    pub fn get_domain(&self) -> Radix2EvaluationDomain<F> {
        self.domain
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
            for i in 0..trace.len() {
                let c = a + b;
                trace.add_row(i, vec![a, b]);
                a = b;
                b = c;
            }

            // set output constrains
            trace.add_boundary_constrain(self.step - 1, 1);

            // add transition constrains
            trace.add_transition_constrain(Box::new(move |trace_polys| {
                trace_polys[0].clone()
                    * DensePolynomial::from_coefficients_vec(vec![trace.omega.clone()])
                    - trace_polys[1].clone()
            }));
            trace.add_transition_constrain(Box::new(move |trace_polys| {
                trace_polys[1].clone() * trace.omega.clone()
                    - (trace_polys[0].clone() + trace_polys[1].clone())
            }));

            trace
        }
    }

    #[test]
    fn test_air_trace() {
        let third_fibonacci = FibonacciClaim {
            step: 3,
            output: Goldilocks::from(3),
        };
        let trace = third_fibonacci.trace(&Witness);
        assert_eq!(trace.len(), 4);
        assert_eq!(trace.width(), 2);
        assert_eq!(*trace.get_value(0, 0), ONE);
        assert_eq!(
            *trace.get_value(third_fibonacci.step - 1, 1),
            third_fibonacci.output
        );
        assert_eq!(
            *trace.get_value(third_fibonacci.step, 0),
            third_fibonacci.output
        );

        let fourth_fib = FibonacciClaim {
            step: 4,
            output: Goldilocks::from(5),
        };
        let trace = fourth_fib.trace(&Witness);
        assert_eq!(trace.len(), 4);
        assert_eq!(trace.width(), 2);
        assert_eq!(*trace.get_value(0, 1), ONE);
        assert_eq!(
            *trace.get_value(third_fibonacci.step - 1, 1),
            third_fibonacci.output
        );
    }

    #[test]
    fn test_air_trace_polynomials() {
        let claim = FibonacciClaim {
            step: 3,
            output: Goldilocks::from(3),
        };
        let trace = claim.trace(&Witness);
        let trace_polys = trace.get_trace_polys();
        let domain = Radix2EvaluationDomain::<Goldilocks>::new(trace.len()).unwrap();

        // check trace polynomials are computed correctly
        for i in 0..claim.step {
            let row = domain.element(i);
            assert_eq!(*trace.get_value(i, 0), trace_polys[0].evaluate(&row));
            assert_eq!(*trace.get_value(i, 1), trace_polys[1].evaluate(&row));
        }
    }

    #[test]
    fn test_air_constrains() {
        let claim = FibonacciClaim {
            step: 3,
            output: Goldilocks::from(3),
        };
        let trace = claim.trace(&Witness);
        let constrains = trace.derive_constrains();
        let domain = constrains.get_domain();
        assert_eq!(constrains.transition_constrains, 2);

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
        for i in 0..trace.len() - 1 {
            let w_i = domain.element(i);
            assert_eq!(carry_over_constrain.evaluate(&w_i), ZERO);
            assert_eq!(sum_constrain.evaluate(&w_i), ZERO);
        }
    }
}
