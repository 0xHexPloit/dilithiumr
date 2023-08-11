use std::ops::{AddAssign, Index, IndexMut, Mul};
use std::slice::IterMut;
use crate::algebra::poly::Polynomial;

#[derive(Debug, Copy, Clone)]
pub struct PolyVec<const N: usize> {
    data: [Polynomial; N]
}

impl <const N: usize> PolyVec<N> {
    pub fn iter_mut(&mut self) -> IterMut<Polynomial> {
        self.data.iter_mut()
    }

    pub fn to_ntt(&self) -> Self {
        let mut clone = self.clone();
        for i in 0..N {
            clone[i].ntt();
        }
        clone
    }

    pub fn caddq(&mut self) {
        for i in 0..N {
            self[i].caddq();
        }
    }

    pub fn power2round(&self) -> (PolyVec<N>, PolyVec<N>) {
        let mut poly_vec_zero = PolyVec::<N>::default();
        let mut poly_vec_one = PolyVec::<N>::default();

        for i in 0..N {
            let (poly_zero, poly_one) = self[i].power2round();
            poly_vec_zero[i] = poly_zero;
            poly_vec_one[i] = poly_one;
        }

        (poly_vec_zero, poly_vec_one)
    }
}

impl <const N: usize>Default for PolyVec<N> {
    fn default() -> Self {
        let data = [Polynomial::default(); N];
        Self {
            data
        }
    }
}

impl <const N: usize>Index<usize> for PolyVec<N> {
    type Output = Polynomial;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl <const N: usize> IndexMut<usize> for PolyVec<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}


impl <const N: usize> AddAssign<&PolyVec<N>> for PolyVec<N> {
    fn add_assign(&mut self, rhs: &PolyVec<N>) {
        for i in 0..N {
            self[i] += &rhs[i]
        }
    }
}

impl <const N: usize> Mul<&PolyVec<N>> for PolyVec<N> {
    type Output = Polynomial;
    fn mul(self, rhs: &PolyVec<N>) -> Self::Output {
        let mut poly = Polynomial::default();

        for i in 0..N {
            poly += &(self[i] * &rhs[i]);
        }

        poly.reduce();
        poly.ntt_to_montgomery();
        poly
    }
}