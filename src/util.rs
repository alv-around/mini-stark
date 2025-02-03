use ark_ff::FftField;
use core::ops::{BitAnd, Sub};
use num_traits::{One, Zero};
use std::cmp::PartialEq;

pub fn is_power_of_two<T>(number: T) -> bool
where
    T: BitAnd<Output = T> + Sub<Output = T> + PartialEq + Copy + One + Zero,
{
    number & (number - T::one()) == T::zero()
}

pub fn get_nth_root<F: FftField>(nth_root: u64) -> F {
    assert!(is_power_of_two(nth_root));
    FftField::get_root_of_unity(nth_root).unwrap()
}
