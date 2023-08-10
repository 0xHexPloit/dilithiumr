use crate::constants::{Q, Q_INVERSE};

pub fn montgomery_reduce(value: i64) -> i32 {
    let reduction = (value as i32).wrapping_mul(Q_INVERSE);
    ((value - reduction as i64 * Q as i64) >> 32) as i32
}


#[cfg(test)]
mod tests {
    use crate::constants::Q;
    use crate::reduce::montgomery_reduce;

    #[test]
    fn test_should_output_the_correct_residue() {
        let z_value: i32 = Q as i32 + 10;
        let r_inverse = 8265825;
        let residue = (z_value as i64 * r_inverse as i64).rem_euclid(Q as i64);
        let output = montgomery_reduce(z_value as i64);

        assert_eq!(residue as i32, output.rem_euclid(Q as i32));
        assert!(output.unsigned_abs() <= (z_value as u64 >> 32) as u32  + (Q >> 2u32))
    }
}