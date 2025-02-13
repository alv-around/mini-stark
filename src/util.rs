use core::ops::{BitAnd, Sub};
use num_traits::{One, Zero};

pub fn is_power_of_two<T>(number: T) -> bool
where
    T: BitAnd<Output = T> + Sub<Output = T> + PartialOrd + Copy + One + Zero,
{
    if number.gt(&T::one()) {
        return number & (number - T::one()) == T::zero();
    }
    false
}

pub fn logarithm_of_two_k<'a, const N: usize>(number: usize) -> Result<usize, &'a str> {
    if !is_power_of_two(number) {
        return Err("number if not a power of 2");
    }
    let power_of_two = number.trailing_zeros() as usize;
    if power_of_two % N != 0 {
        return Err("number if not a power of base");
    }
    Ok(power_of_two / N)
}
