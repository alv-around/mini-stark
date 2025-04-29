use ark_ff::fields::{Fp2Config, MontFp};
use ark_ff::fields::{MontBackend, MontConfig};
use ark_ff::{FftField, Field, Fp, Fp2ConfigWrapper, PrimeField, QuadExtField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};

// ========== Stark Field Config ==========
pub trait StarkField {
    type Base: PrimeField;
    type Extension: Field<BasePrimeField = Self::Base> + Clone + FftField;

    // Compile-time safety check
    fn SOUNDNESS_CHECK(&self) {
        // Verify extension field size > 2^100 (Goldilocks² has 128-bit modulus)
        assert!(
            <Self::Base as PrimeField>::MODULUS_BIT_SIZE as u64
                * <Self::Extension as Field>::extension_degree()
                > 100
        );
    }

    fn extend_poly(poly: &DensePolynomial<Self::Base>) -> DensePolynomial<Self::Extension> {
        // Convert the base field coefficients to the extension field
        let extension_coefficients = poly
            .coeffs()
            .iter()
            .map(|&coeff| Self::Extension::from_base_prime_field(coeff))
            .collect();

        DensePolynomial::from_coefficients_vec(extension_coefficients)
    }
}

// Goldilocks implementation
pub struct Goldilocks;

impl StarkField for Goldilocks {
    type Base = GoldilocksFp;
    type Extension = GoldilocksFp2;
}

#[derive(MontConfig)]
#[modulus = "18446744069414584321"]
#[generator = "7"]
pub struct GoldilocksConfig;
pub type GoldilocksFp = Fp<MontBackend<GoldilocksConfig, 1>, 1>;

// field extension implementation taken from: https://github.com/WizardOfMenlo/whir/blob/main/src/crypto/fields.rs
pub type GoldilocksFp2 = QuadExtField<GoldilocksQuadraticExtensionConfig>;
pub type GoldilocksQuadraticExtensionConfig = Fp2ConfigWrapper<GoldilocksFp2Config>;
pub struct GoldilocksFp2Config;
impl Fp2Config for GoldilocksFp2Config {
    type Fp = GoldilocksFp;
    const NONRESIDUE: Self::Fp = MontFp!("7");
    const FROBENIUS_COEFF_FP2_C1: &'static [Self::Fp] = &[
        // Fq(7)**(((q^0) - 1) / 2)
        MontFp!("1"),
        // Fq(7)**(((q^1) - 1) / 2)
        MontFp!("18446744069414584320"),
    ];
}

#[derive(MontConfig)]
#[modulus = "2013265921"]
#[generator = "7"]
pub struct BabyBearConfig;
pub type BabyBear = Fp<MontBackend<BabyBearConfig, 1>, 1>;

#[cfg(test)]
mod test {
    // super::*;
    //
    // #[test]
    // fn test_evaluate_on_extension() {
    //     use ark_test_curves::bls12_381::{Fq as ExtensionField, Fr as BaseField};
    //
    //     // Define a polynomial with coefficients in the base field (Fr)
    //     let coefficients = vec![
    //         BaseField::from(1u64),
    //         BaseField::from(2u64),
    //         BaseField::from(3u64),
    //     ];
    //     // The polynomial is conceptually 1 + 2x + 3x^2
    //
    //     // Choose a point in the extension field (Fq) to evaluate at
    //     let evaluation_point = ExtensionField::from_random(&mut ark_std::test_rng());
    //
    //     // Evaluate the polynomial at the extension field point
    //     let result = evaluate_base_poly_on_extension_point(&coefficients, evaluation_point);
    //
    //     println!("Polynomial (base coefficients): 1 + 2x + 3x^2");
    //     println!("Evaluation Point (in Fq): {}", evaluation_point);
    //     println!("Result of Evaluation (in Fq): {}", result);
    //
    //     // Manually compute the result by lifting coefficients
    //     let one_in_fq = ExtensionField::from_base_prime_field(BaseField::from(1u64));
    //     let two_in_fq = ExtensionField::from_base_prime_field(BaseField::from(2u64));
    //     let three_in_fq = ExtensionField::from_base_prime_field(BaseField::from(3u64));
    //
    //     let manual_result =
    //         one_in_fq + two_in_fq * evaluation_point + three_in_fq * evaluation_point.square();
    //
    //     assert_eq!(result, manual_result);
    // }
}
