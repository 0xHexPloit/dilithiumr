use crate::constants::N;

#[derive(Debug, Clone, Copy)]
pub struct Polynomial {
    coefficients: [usize; N],
}

impl Polynomial {
    pub fn zeroes() -> Self {
        Self {
            coefficients: [0; N]
        }
    }

    pub fn set(&mut self, idx: usize, value: usize) {
        self.coefficients[idx] = value;
    }
}


impl Default for Polynomial {
    fn default() -> Self {
        Self::zeroes()
    }
}