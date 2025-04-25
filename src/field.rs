use ark_ff::fields::{Fp2, Fp2Config, MontFp};
use ark_ff::fields::{MontBackend, MontConfig};
use ark_ff::{FftField, Fp, PrimeField};

#[derive(MontConfig)]
#[modulus = "2013265921"]
#[generator = "7"]
pub struct BabyBearConfig;
pub type BabyBear = Fp<MontBackend<BabyBearConfig, 1>, 1>;

#[derive(MontConfig)]
#[modulus = "18446744069414584321"]
#[generator = "7"]
pub struct GoldilocksConfig;
pub type Goldilocks = Fp<MontBackend<GoldilocksConfig, 1>, 1>;

// field extension implementation taken from: https://github.com/WizardOfMenlo/whir/blob/main/src/crypto/fields.rs
pub type GoldilocksQuadraticExtension = Fp2<GoldilocksQuadraticExtensionConfig>;
pub struct GoldilocksQuadraticExtensionConfig;
impl Fp2Config for GoldilocksQuadraticExtensionConfig {
    type Fp = Goldilocks;
    const NONRESIDUE: Self::Fp = MontFp!("7");
    const FROBENIUS_COEFF_FP2_C1: &'static [Self::Fp] = &[
        // Fq(7)**(((q^0) - 1) / 2)
        MontFp!("1"),
        // Fq(7)**(((q^1) - 1) / 2)
        MontFp!("18446744069414584320"),
    ];
}

// ---- Helper Trait ----
pub trait ExtensionFieldBase {
    type Base: PrimeField;
}

impl<Config: Fp2Config> ExtensionFieldBase for Fp2<Config> {
    type Base = Config::Fp;
}

// ---- StarkField Trait ----
pub trait StarkField<F: FftField>
where
    Self::BaseField: FftField,
    Self::ExtensionField: ExtensionFieldBase<Base = Self::BaseField>,
{
    type BaseField;
    type ExtensionField;
}

struct GoldilocksStark;

impl StarkField<Goldilocks> for GoldilocksStark {
    type BaseField = Goldilocks;
    type ExtensionField = GoldilocksQuadraticExtension;
}
