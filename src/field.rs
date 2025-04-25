use ark_ff::fields::{Fp2Config, MontFp};
use ark_ff::fields::{MontBackend, MontConfig};
use ark_ff::{FftField, Fp, Fp2ConfigWrapper, PrimeField, QuadExtConfig};

#[derive(MontConfig)]
#[modulus = "2013265921"]
#[generator = "7"]
pub struct BabyBearConfig;
pub type BabyBear = Fp<MontBackend<BabyBearConfig, 1>, 1>;

#[derive(MontConfig)]
#[modulus = "18446744069414584321"]
#[generator = "7"]
pub struct GoldilocksConfig;
pub type GoldilocksFp = Fp<MontBackend<GoldilocksConfig, 1>, 1>;

// field extension implementation taken from: https://github.com/WizardOfMenlo/whir/blob/main/src/crypto/fields.rs
pub type GoldilocksQuadraticExtension = Fp2ConfigWrapper<GoldilocksFp2Config>;
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
// ========== Stark Field Config ==========
pub trait StarkField {
    type Base: PrimeField + FftField;
    type Extension: QuadExtConfig;

    // Compile-time safety check
    fn SOUNDNESS_CHECK(&self) {
        // Verify extension field size > 2^100 (Goldilocks² has 128-bit modulus)
        assert!(
            <Self::Base as PrimeField>::MODULUS_BIT_SIZE as usize
                * <Self::Extension as QuadExtConfig>::DEGREE_OVER_BASE_PRIME_FIELD
                > 100
        );
    }
}

// Goldilocks implementation
pub struct Goldilocks;

impl StarkField for Goldilocks {
    type Base = GoldilocksFp;
    type Extension = GoldilocksQuadraticExtension;
}
