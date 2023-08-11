use crate::constants::{Q, Q_INVERSE};

/// For a finite field element (named 'val' here) such that -2^31 * Q <= val <= 2^31 * Q, this
/// function computes r such that r = a * (2^32)^-1 mod Q and -Q < r < Q.
pub fn montgomery_reduce(val: i64) -> i32 {
    let reduction = (val as i32).wrapping_mul(Q_INVERSE);
    ((val - reduction as i64 * Q as i64) >> 32) as i32
}

/// For a finite field element (named 'val' here) such that val < Q, this function computes
///r such that r = a [Q] and -6283009 <= r <= 6283007.
///
/// # Arguments
///
/// * `val` - A finite field element
pub fn reduce32(val: i32) -> i32 {
    let t = (val + (1 << 22)) >> 23;
    val - t * Q as i32
}

///
///
pub fn caddq(val: i32) -> i32 {
    val + ((val >> 31) & Q as i32)
}

#[cfg(test)]
mod tests {
    use rand::distributions::{Distribution, Uniform};
    use rand::thread_rng;
    use crate::constants::Q;
    use crate::reduce::{caddq, montgomery_reduce, reduce32};

    #[test]
    fn test_should_output_the_correct_residue() {
        let z_value: i32 = Q as i32 + 10;
        let r_inverse = 8265825;
        let residue = (z_value as i64 * r_inverse as i64).rem_euclid(Q as i64);
        let output = montgomery_reduce(z_value as i64);

        assert_eq!(residue as i32, output.rem_euclid(Q as i32));
        assert!(output.unsigned_abs() <= (z_value as u64 >> 32) as u32  + (Q >> 2u32))
    }


    #[test]
    fn test_should_output_a_value_that_respects_the_reduce_constraints() {
        let mut rng = thread_rng();
        let dist = Uniform::from((-1i32 * Q as i32 + 1)..Q as i32);
        let mut buff = [0i32; 10];
        buff.iter_mut().for_each(|val| *val = dist.sample(&mut rng));
        buff.iter().for_each(|val| {
            let output = reduce32(*val);
            assert_eq!(val.rem_euclid(Q as i32), output.rem_euclid(Q as i32));
            assert!(-6283009 <= output && output <= 6283007)
        })

    }

    #[test]
    fn test_should_output_the_same_value_as_input_as_it_is_positive() {
        let input = 5;
        let output = caddq(input);
        assert_eq!(input, output);
    }

    #[test]
    fn test_should_add_q_to_input_as_it_is_negative() {
        let input = -5;
        let output = caddq(input);
        assert_eq!(output, input + Q as i32);
    }
}