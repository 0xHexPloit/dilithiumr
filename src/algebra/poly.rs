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

    pub fn uniform_random(rho: &[u8], nonce: usize) -> Self {
        Self::zeroes()
    }
}


impl Default for Polynomial {
    fn default() -> Self {
        Self::zeroes()
    }
}