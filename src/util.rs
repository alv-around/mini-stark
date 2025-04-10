use core::ops::{BitAnd, Sub};
use num_traits::{One, Zero};

pub fn is_power_of_two<T>(number: T) -> bool
where
    T: BitAnd<Output = T> + Sub<Output = T> + PartialOrd + Copy + One + Zero,
{
    if number.gt(&T::zero()) {
        return number & (number - T::one()) == T::zero();
    }
    false
}

pub fn logarithm_of_two_k(number: usize, base: usize) -> Result<usize, &'static str> {
    assert!(is_power_of_two(base));
    let log_n = base.trailing_zeros() as usize;

    if !is_power_of_two(number) {
        return Err("number if not a power of 2");
    }
    let power_of_two = number.trailing_zeros() as usize;
    if power_of_two % log_n != 0 {
        return Err("number if not a power of base");
    }
    Ok(power_of_two / log_n)
}

pub fn ceil_log2_k(number: usize, base: usize) -> usize {
    assert!(is_power_of_two(base));
    assert!(number != 0);
    if number == 1 {
        return 1;
    }
    let log2_base = base.trailing_zeros() as usize;
    let log2_number = number.trailing_zeros() as usize;
    if is_power_of_two(number) && log2_number % log2_base == 0 {
        log2_number
    } else {
        let next_power_2 = usize::BITS - number.leading_zeros();
        next_power_2.div_ceil(log2_base as u32) as usize * log2_base
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_is_power_of_two() {
        assert!(is_power_of_two(1u8));
        assert!(is_power_of_two(2u8));
        assert!(is_power_of_two(32u16));
        assert!(is_power_of_two(128u32));
        assert!(is_power_of_two(512u64));
        assert!(is_power_of_two(1024usize));

        assert!(!is_power_of_two(0u16));
        assert!(!is_power_of_two(24u32));
        assert!(!is_power_of_two(48usize));
    }

    #[test]
    fn test_logarithm_of_two_k() {
        let not_power_of_2 = Err("number if not a power of 2");
        let not_power_of_2k = Err("number if not a power of base");

        assert_eq!(Ok(5), logarithm_of_two_k(32, 2));
        assert_eq!(not_power_of_2, logarithm_of_two_k(6, 2));

        assert_eq!(Ok(4), logarithm_of_two_k(256, 4));
        assert_eq!(not_power_of_2, logarithm_of_two_k(12, 4));
        assert_eq!(not_power_of_2k, logarithm_of_two_k(32, 4));

        assert_eq!(Ok(3), logarithm_of_two_k(512, 8));
        assert_eq!(not_power_of_2, logarithm_of_two_k(15, 8));
        assert_eq!(not_power_of_2k, logarithm_of_two_k(16, 8));

        assert_eq!(Ok(2), logarithm_of_two_k(256, 16));
        assert_eq!(not_power_of_2, logarithm_of_two_k(48, 16));
        assert_eq!(not_power_of_2k, logarithm_of_two_k(64, 16));
    }

    #[test]
    fn test_ceil_log_power_two() {
        // assert_eq!(1, ceil_log2_k(1, 2));
        assert_eq!(1, ceil_log2_k(2, 2));
        assert_eq!(5, ceil_log2_k(21, 2));
        assert_eq!(5, ceil_log2_k(32, 2));

        assert_eq!(2, ceil_log2_k(4, 4));
        assert_eq!(2, ceil_log2_k(3, 4));
        assert_eq!(4, ceil_log2_k(13, 4));
        assert_eq!(6, ceil_log2_k(21, 4));
    }
}
