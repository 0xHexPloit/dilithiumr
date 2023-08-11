use crate::constants::D;

/// Given a finite field element (named 'val' here), this function computes (and returns) two
/// numbers (a0, a1) such that a mod ^{+Q} = a1*2^D + a0 with -2^{D-1} < a0 <= 2^{D-1}.
///
/// This function also expects a to be standard representative
///
/// # Arguments:
/// * `val` - A finite field element
pub fn power2round(val: i32) -> (i32, i32) {
    let a_one = (val + (1 << (D-1)) - 1) >> D;
    let a_zero = val - (a_one << D);
    return (a_zero, a_one)
}

#[cfg(test)]
mod tests {
    use crate::round::power2round;

    #[test]
    fn test_it_should_output_correct_values() {
        let values = [
            [100, 0, 100],
            [2i32.pow(14), 2, 0]
        ];

        for item in values {
            let (a_zero, a_one) = power2round(item[0]);
            assert_eq!(a_zero, item[2]);
            assert_eq!(a_one, item[1])
        }
    }
}