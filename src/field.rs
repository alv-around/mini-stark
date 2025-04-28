use ark_ff::fields::{Fp2Config, MontFp};
use ark_ff::fields::{MontBackend, MontConfig};
use ark_ff::{FftField, Field, Fp, Fp2ConfigWrapper, PrimeField, QuadExtField};
use std::marker::PhantomData;

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
    type Extension: Field<BasePrimeField = Self::Base>;

    // Compile-time safety check
    fn SOUNDNESS_CHECK(&self) {
        // Verify extension field size > 2^100 (Goldilocks² has 128-bit modulus)
        assert!(
            <Self::Base as PrimeField>::MODULUS_BIT_SIZE as u64
                * <Self::Extension as Field>::extension_degree()
                > 100
        );
    }

    // fn from_base_prime_field_elems(
    //     elems: impl IntoIterator<Item = Self::Base>,
    // ) -> Option<Self::Extension>;
}

// Goldilocks implementation
pub struct Goldilocks {
    base: PhantomData<GoldilocksFp>,
    pub(crate) extension: PhantomData<QuadExtField<GoldilocksQuadraticExtension>>,
}

impl Goldilocks {
    pub fn new() -> Self {
        Self {
            base: PhantomData::<GoldilocksFp>,
            extension: PhantomData::<QuadExtField<GoldilocksQuadraticExtension>>,
        }
    }
}

impl StarkField for Goldilocks {
    type Base = GoldilocksFp;
    type Extension = QuadExtField<GoldilocksQuadraticExtension>;

    // fn from_base_prime_field_elems(
    //     elems: impl IntoIterator<Item = Self::Base>,
    // ) -> Option<Self::Extension> {
    //     QuadExtField::<Self::Extension>::from_base_prime_field_elems(elems)
    // }
}
