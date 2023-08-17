use crate::algebra::poly::Polynomial;
use crate::algebra::vec::PolyVec;
use std::ops::{Index, IndexMut, Mul};

#[derive(Debug)]
pub struct PolyMatrix<const K: usize, const L: usize> {
    data: [PolyVec<L>; K],
}

impl<const K: usize, const L: usize> PolyMatrix<K, L> {
    pub fn zeroes() -> Self {
        let data = [PolyVec::<L>::default(); K];
        Self { data }
    }

    pub fn expand_a(rho: &[u8]) -> PolyMatrix<K, L> {
        let mut matrix = PolyMatrix::<K, L>::default();
        for i in 0..K {
            for j in 0..L {
                let poly = &mut matrix[i][j];
                Polynomial::fill_uniform_random(rho, ((i << 8) + j) as u16, poly);
            }
        }

        matrix
    }
}

impl<const K: usize, const L: usize> Default for PolyMatrix<K, L> {
    fn default() -> Self {
        Self::zeroes()
    }
}

impl<const K: usize, const L: usize> Index<usize> for PolyMatrix<K, L> {
    type Output = PolyVec<L>;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<const K: usize, const L: usize> IndexMut<usize> for PolyMatrix<K, L> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl<const K: usize, const L: usize> Mul<&PolyVec<L>> for &PolyMatrix<K, L> {
    type Output = PolyVec<K>;
    fn mul(self, rhs: &PolyVec<L>) -> Self::Output {
        let mut poly_vec = PolyVec::<K>::default();

        for i in 0..K {
            let row = self.data[i];
            poly_vec[i] = row * rhs;
        }

        poly_vec
    }
}
