use crate::constants::{Q, Q_INVERSE};

pub fn montgomery_reduce(value: i64) -> i32 {
    let reduction = (value as i32).wrapping_mul(Q_INVERSE) as i64;
    ((value  - reduction * Q as i64) >> 32) as i32
}