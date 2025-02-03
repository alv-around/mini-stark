use ark_ff::fields::{MontBackend, MontConfig};
use ark_ff::{FftField, Fp};

#[derive(MontConfig)]
#[modulus = "2147483647"]
#[generator = "7"]
pub struct Mersenne31Config;
pub type M31 = Fp<MontBackend<Mersenne31Config, 1>, 1>;

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

pub enum SmallPrimeField {
    M31,
    BabyBear,
    Goldilocks,
}
