use crate::algebra::poly::Polynomial;
use std::ops::{AddAssign, Index, IndexMut, Mul, Sub};
use std::slice::IterMut;

#[derive(Debug, Copy, Clone)]
pub struct PolyVec<const N: usize> {
    data: [Polynomial; N],
}

impl<const N: usize> PolyVec<N> {
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

    pub fn reduce(&mut self) {
        for i in 0..N {
            self[i].reduce()
        }
    }

    /// This function returns a Boolean indicating whether the infinite norm of the vector of
    /// polynomials is greater than or equal to a specific bound.
    ///
    /// # Arguments:
    /// * `bound` - The upper bound
    pub fn check_infinite_norm(&self, bound: usize) -> bool {
        for poly in &self.data {
            if poly.check_infinite_norm(bound) {
                return true;
            }
        }

        false
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

    pub fn expand_mask<const ARR_SIZE: usize, const GAMMA_ONE: usize>(
        rho_prime: &[u8],
        kappa: u16,
    ) -> PolyVec<N> {
        let mut vec = PolyVec::default();
        let mut nonce = kappa;

        for i in 0..N {
            Polynomial::fill_uniform_gamma::<ARR_SIZE, GAMMA_ONE>(rho_prime, nonce, &mut vec[i]);
            nonce += 1;
        }

        vec
    }

    pub fn decompose<const GAMMA_TWO: usize>(vec: &PolyVec<N>) -> (PolyVec<N>, PolyVec<N>) {
        let mut r_zero = PolyVec::<N>::default();
        let mut r_one = PolyVec::<N>::default();

        for i in 0..N {
            Polynomial::decompose::<GAMMA_TWO>(&mut r_zero[i], &mut r_one[i], &vec[i]);
        }

        return (r_one, r_zero);
    }

    pub fn make_hint<const GAMMA_TWO: usize>(
        output_vec: &mut PolyVec<N>,
        v_one: &PolyVec<N>,
        v_zero: &PolyVec<N>,
    ) -> usize {
        let mut number_ones = 0;

        for i in 0..N {
            number_ones +=
                Polynomial::make_hint::<GAMMA_TWO>(&mut output_vec[i], &v_one[i], &v_zero[i]);
        }

        number_ones
    }

    pub fn expand_s<const ETA: usize>(rho_prime: &[u8], nonce_init_value: u16) -> PolyVec<N> {
        let mut vec = PolyVec::<N>::default();
        let mut nonce = nonce_init_value;
        vec.iter_mut().for_each(|poly| {
            Polynomial::fill_uniform_eta::<ETA>(rho_prime, nonce, poly);
            nonce += 1
        });
        vec
    }
}

impl<const N: usize> Default for PolyVec<N> {
    fn default() -> Self {
        let data = [Polynomial::default(); N];
        Self { data }
    }
}

impl<const N: usize> Index<usize> for PolyVec<N> {
    type Output = Polynomial;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<const N: usize> IndexMut<usize> for PolyVec<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl<const N: usize> AddAssign<&PolyVec<N>> for PolyVec<N> {
    fn add_assign(&mut self, rhs: &PolyVec<N>) {
        for i in 0..N {
            self[i] += &rhs[i]
        }
    }
}

impl<const N: usize> Mul<&PolyVec<N>> for PolyVec<N> {
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

impl<const N: usize> Mul<&Polynomial> for PolyVec<N> {
    type Output = PolyVec<N>;

    fn mul(self, rhs: &Polynomial) -> Self::Output {
        let mut vec = PolyVec::<N>::default();

        for i in 0..N {
            vec[i] = self[i] * rhs;
            vec[i].ntt_to_montgomery();
            vec[i].reduce();
        }

        vec
    }
}

impl<const N: usize> Sub<&PolyVec<N>> for PolyVec<N> {
    type Output = PolyVec<N>;
    fn sub(self, rhs: &PolyVec<N>) -> Self::Output {
        let mut vec = PolyVec::default();

        for i in 0..N {
            vec[i] = self[i] - &rhs[i]
        }

        vec
    }
}
