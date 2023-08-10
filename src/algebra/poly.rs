use std::ops::Index;
use crate::constants::N;

#[derive(Debug, Clone, Copy)]
pub struct Polynomial {
    coefficients: [i32; N],
}

impl Polynomial {
    pub fn zeroes() -> Self {
        Self {
            coefficients: [0; N]
        }
    }

    pub fn set(&mut self, idx: usize, value: i32) {
        self.coefficients[idx] = value;
    }
}


impl Default for Polynomial {
    fn default() -> Self {
        Self::zeroes()
    }
}

impl Index<usize> for Polynomial {
    type Output = i32;
    fn index(&self, index: usize) -> &Self::Output {
        &self.coefficients[index]
    }
}