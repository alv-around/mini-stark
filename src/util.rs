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

pub fn logarithm_of_two_k<const N: usize>(number: usize) -> Result<usize, &'static str> {
    assert!(is_power_of_two(N));
    let log_n = N.trailing_zeros() as usize;

    if !is_power_of_two(number) {
        return Err("number if not a power of 2");
    }
    let power_of_two = number.trailing_zeros() as usize;
    if power_of_two % log_n != 0 {
        return Err("number if not a power of base");
    }
    Ok(power_of_two / log_n)
}

#[cfg(test)]
mod test {
    use super::*;

    const TWO: usize = 2;
    const FOUR: usize = 4;
    const EIGHT: usize = 8;
    const SIXTEEN: usize = 16;

    #[test]
    fn test_is_power_of_two() {
        assert!(is_power_of_two(2u8));
        assert!(is_power_of_two(32u16));
        assert!(is_power_of_two(128u32));
        assert!(is_power_of_two(512u64));
        assert!(is_power_of_two(1024usize));

        assert!(!is_power_of_two(1u8));
        assert!(!is_power_of_two(0u16));
        assert!(!is_power_of_two(24u32));
        assert!(!is_power_of_two(48usize));
    }

    #[test]
    fn test_logarithm_of_two_k() {
        let not_power_of_2 = Err("number if not a power of 2");
        let not_power_of_2k = Err("number if not a power of base");

        assert_eq!(Ok(5), logarithm_of_two_k::<TWO>(32));
        assert_eq!(not_power_of_2, logarithm_of_two_k::<TWO>(6));

        assert_eq!(Ok(4), logarithm_of_two_k::<FOUR>(256));
        assert_eq!(not_power_of_2, logarithm_of_two_k::<FOUR>(12));
        assert_eq!(not_power_of_2k, logarithm_of_two_k::<FOUR>(32));

        assert_eq!(Ok(3), logarithm_of_two_k::<EIGHT>(512));
        assert_eq!(not_power_of_2, logarithm_of_two_k::<EIGHT>(15));
        assert_eq!(not_power_of_2k, logarithm_of_two_k::<EIGHT>(16));

        assert_eq!(Ok(2), logarithm_of_two_k::<SIXTEEN>(256));
        assert_eq!(not_power_of_2, logarithm_of_two_k::<SIXTEEN>(48));
        assert_eq!(not_power_of_2k, logarithm_of_two_k::<SIXTEEN>(64));
    }
}
