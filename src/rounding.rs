use crate::constants::{D, Q};

/// Given a finite field element (named 'val' here), this function computes (and returns) two
/// numbers (a0, a1) such that a mod ^{+Q} = a1*2^D + a0 with -2^{D-1} < a0 <= 2^{D-1}.
///
/// This function also expects a to be standard representative
///
/// # Arguments:
/// * `val` - A finite field element
pub fn power2round(val: i32) -> (i32, i32) {
    let a_one = (val + (1 << (D - 1)) - 1) >> D;
    let a_zero = val - (a_one << D);
    return (a_zero, a_one);
}

/// Given a finite field element (named 'val' here), this function computes high and low bits a0,
/// a1 such that a mod^+ Q = a1*ALPHA + a0 with -ALPHA/2 < a0 <= ALPHA/2 except
/// if a1 = (Q-1)/ALPHA where we set a1 = 0 and -ALPHA/2 <= a0 = a mod^+ Q - Q < 0.
///
/// # Arguments:
/// * `val` - A finite field element

pub fn decompose<const GAMMA_TWO: usize>(val: i32) -> (i32, i32) {
    let mut a_one: i32 = (val + 127) >> 7;

    if GAMMA_TWO == ((Q - 1) / 32) as usize {
        a_one = (a_one * 1025 + (1 << 21)) >> 22;
        a_one &= 15;
    } else {
        a_one = (a_one * 11275 + (1 << 23)) >> 24;
        a_one ^= ((43 - a_one) >> 31) & a_one;
    }

    let mut a_zero = val - a_one * 2 * GAMMA_TWO as i32;
    a_zero -= (((Q as i32 - 1) / 2 - a_zero) >> 31) & Q as i32;

    return (a_zero, a_one);
}

/// This function computes hint bit indicating whether the low bits of the
///  input element overflow into the high bits. If the function returns 1, it means that there is
/// an overflow.
///
/// # Arguments
///
/// * `a0` - Low bits of input element
/// * `a1` - High bits of input element
pub fn make_hint<const GAMMA_TWO: usize>(a0: i32, a1: i32) -> i32 {
    return if a0 > GAMMA_TWO as i32
        || a0 < -1 * GAMMA_TWO as i32
        || (a0 == -1 * GAMMA_TWO as i32 && a1 != 0)
    {
        1
    } else {
        0
    };
}

#[cfg(test)]
mod tests {
    use crate::rounding::power2round;

    #[test]
    fn test_it_should_output_correct_values() {
        let values = [[100, 0, 100], [2i32.pow(14), 2, 0]];

        for item in values {
            let (a_zero, a_one) = power2round(item[0]);
            assert_eq!(a_zero, item[2]);
            assert_eq!(a_one, item[1])
        }
    }
}
