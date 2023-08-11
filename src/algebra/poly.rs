use sha3::digest::{XofReader};
use std::ops::{AddAssign, Index, IndexMut, Mul};
use crate::constants::{D, N, POLY_UNIFORM_ETA_2_BUFLEN, POLY_UNIFORM_ETA_4_BUFLEN, REJ_UNIFORM_BUFLEN, STREAM128_BLOCKBYTES, STREAM256_BLOCKBYTES, ZETAS};
use crate::helper::{shake_128_reader, shake_256_reader};
use crate::reduce::{caddq, montgomery_reduce, reduce32};
use crate::rejection::{rejection_eta, rejection_uniform};
use crate::round::power2round;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Polynomial {
    pub (crate) coefficients: [i32; N],
}

impl Polynomial {
    pub fn zeroes() -> Self {
        Self {
            coefficients: [0; N]
        }
    }

    pub fn caddq(&mut self) {
        for i in 0..N {
            self[i] = caddq(self[i]);
        }
    }

    pub fn power2round(&self) -> (Polynomial, Polynomial) {
        let mut poly_zero = Polynomial::default();
        let mut poly_one = Polynomial::default();

        for i in 0..N {
            let (a_zero, a_one) = power2round(self[i]);
            poly_zero[i] = a_zero;
            poly_one[i] = a_one;
        }

        (poly_zero, poly_one)
    }


    pub fn reduce(&mut self) {
        for i in 0..N {
            self[i] = reduce32(self[i])
        }
    }

    pub fn pack_t1(buf: &mut [u8], poly: &Polynomial) {
        for i in 0..N/4 {
            buf[5*i+0] = (poly[4*i+0] >> 0) as u8;
            buf[5*i+1] = (poly[4*i+0] >> 8) as u8 | (poly[4*i+1] << 2) as u8;
            buf[5*i+2] = (poly[4*i+1] >> 6) as u8 | (poly[4*i+2] << 4) as u8;
            buf[5*i+3] = (poly[4*i+2] >> 4) as u8 | (poly[4*i+3] << 6) as u8;
            buf[5*i+4] = (poly[4*i+3] >> 2) as u8;
        }
    }

    pub fn pack_eta<const ETA: usize>(buf: &mut [u8], poly: &Polynomial) {
        let mut arr = [0u8; 8];

        if ETA == 2 {
            for i in 0..N/8 {
                arr[0] = (ETA as i32- poly[8*i+0]) as u8;
                arr[1] = (ETA as i32 - poly[8*i+1]) as u8;
                arr[2] = (ETA as i32 - poly[8*i+2]) as u8;
                arr[3] = (ETA as i32 - poly[8*i+3]) as u8;
                arr[4] = (ETA as i32 - poly[8*i+4]) as u8;
                arr[5] = (ETA as i32 - poly[8*i+5]) as u8;
                arr[6] = (ETA as i32 - poly[8*i+6]) as u8;
                arr[7] = (ETA as i32- poly[8*i+7]) as u8;

                buf[3*i+0]  = (arr[0] >> 0) | (arr[1] << 3) | (arr[2] << 6);
                buf[3*i+1]  = (arr[2] >> 2) | (arr[3] << 1) | (arr[4] << 4) | (arr[5] << 7);
                buf[3*i+2]  = (arr[5] >> 1) | (arr[6] << 2) | (arr[7] << 5);
            }
        } else {
            for i in 0..N/2 {
                arr[0] = (ETA as i32 - poly[2*i+0]) as u8;
                arr[1] = (ETA as i32 - poly[2*i+1]) as u8;
                buf[i] = arr[0] | (arr[1] << 4);
            }
        }
    }

    pub fn pack_t0(buf: &mut [u8], poly: &Polynomial) {
        let mut arr = [0i32; 8];
        const UPPER_BOUND: i32 = 1 << (D-1);

        for i in 0..N/8 {
            arr[0] = UPPER_BOUND - poly[8*i+0];
            arr[1] = UPPER_BOUND - poly[8*i+1];
            arr[2] = UPPER_BOUND - poly[8 * i + 2];
            arr[3] = UPPER_BOUND - poly[8 * i + 3];
            arr[4] = UPPER_BOUND - poly[8 * i + 4];
            arr[5] = UPPER_BOUND - poly[8 * i + 5];
            arr[6] = UPPER_BOUND - poly[8 * i + 6];
            arr[7] = UPPER_BOUND - poly[8 * i + 7];

            buf[13 * i + 0] = (arr[0]) as u8;
            buf[13 * i + 1] = (arr[0] >> 8) as u8;
            buf[13 * i + 1] |= (arr[1] << 5) as u8;
            buf[13 * i + 2] = (arr[1] >> 3) as u8;
            buf[13 * i + 3] = (arr[1] >> 11) as u8;
            buf[13 * i + 3] |= (arr[2] << 2) as u8;
            buf[13 * i + 4] = (arr[2] >> 6) as u8;
            buf[13 * i + 4] |= (arr[3] << 7) as u8;
            buf[13 * i + 5] = (arr[3] >> 1) as u8;
            buf[13 * i + 6] = (arr[3] >> 9) as u8;
            buf[13 * i + 6] |= (arr[4] << 4) as u8;
            buf[13 * i + 7] = (arr[4] >> 4) as u8;
            buf[13 * i + 8] = (arr[4] >> 12) as u8;
            buf[13 * i + 8] |= (arr[5] << 1) as u8;
            buf[13 * i + 9] = (arr[5] >> 7) as u8;
            buf[13 * i + 9] |= (arr[6] << 6) as u8;
            buf[13 * i + 10] = (arr[6] >> 2) as u8;
            buf[13 * i + 11] = (arr[6] >> 10) as u8;
            buf[13 * i + 11] |= (arr[7] << 3) as u8;
            buf[13 * i + 12] = (arr[7] >> 5) as u8;
        }
    }

    /// This function converts the coefficients of the polynomial from the NTT domain to
    /// the montgomery domain.
    ///
    /// Be aware that due to the Montgomery reduction optimization, we have:
    ///  p != poly_from_ntt(poly_to_ntt(p))
    ///
    pub fn ntt_to_montgomery(&mut self) {
        let mut layer = 1;
        let mut zeta_index = 256;
        const F: i64 = 41978; // mont**2/256

        while layer < N {
            let mut offset = 0;

            while offset < N {
                zeta_index -= 1;
                let zeta = -ZETAS[zeta_index] as i64;
                let mut j = offset;

                while j < offset + layer {
                    let temp_value = self[j];
                    self[j] = temp_value + self[j + layer];
                    self[j + layer] = temp_value - self[j + layer];
                    self[j + layer] = montgomery_reduce(zeta * self[j + layer] as i64);
                    j  += 1;
                }

                offset = j + layer
            }
            layer <<= 1;
        }

        for j in 0..N {
            self[j] = montgomery_reduce(F * self[j] as i64);
        }
    }

    /// This function converts the coefficients of the polynomial into the NTT domain.
    pub fn ntt(&mut self) {
        let mut zeta_index = 0usize;
        let mut layer = 128u8;

        while layer != 0 {
            let mut offset = 0;
            while offset < N {
                zeta_index += 1;
                let zeta = ZETAS[zeta_index] as i64;

                let mut j = offset;
                while j < offset + layer as usize {
                    let t = montgomery_reduce(zeta * self[j + layer as usize] as i64);

                    self[j + layer as usize] = self[j] - t;
                    self[j] = self[j] + t;

                    j += 1;
                }

                offset = j + layer as usize;
            }
            layer >>= 1;
        }
    }

    fn fill_uniform_random_base<
        const XOF_INPUT_BUFFER_SIZE:usize,
        const ARR_SIZE: usize,
        const STREAM_BLOCK_BYTES: usize
    >(
        rho: &[u8],
        nonce: u16,
        poly:
        &mut Polynomial,
        rejection_fn: impl Fn(&mut Polynomial, usize, usize, &[u8], usize) -> usize,
        get_reader: impl Fn(&[u8]) -> Box<dyn XofReader>
    )
    {
        let mut xof_input = [0u8; XOF_INPUT_BUFFER_SIZE];
        xof_input[..XOF_INPUT_BUFFER_SIZE - 2].copy_from_slice(rho);
        xof_input[XOF_INPUT_BUFFER_SIZE - 2] = (nonce & 0xFF) as u8;
        xof_input[XOF_INPUT_BUFFER_SIZE - 1] = ((nonce >> 8)  & 0xFF) as u8;

        let mut reader = get_reader(&xof_input);
        let mut buffer = [0u8; ARR_SIZE];
        let mut buffer_len = ARR_SIZE;
        reader.read(&mut buffer);

        let mut counter = rejection_fn(
            poly,
            N,
            0,
            &buffer,
            buffer_len
        );

        while counter < N {
            let offset = buffer_len % 3;
            for i in 0..offset {
                buffer[i] = buffer[buffer_len - offset + 1];
            }
            reader.read(&mut buffer[offset..STREAM_BLOCK_BYTES]);
            buffer_len = STREAM_BLOCK_BYTES + offset;
            counter += rejection_fn(
                poly,
                N - counter,
                counter,
                &buffer,
                buffer_len
            );
        }
    }

    pub fn fill_uniform_random_using_eta<const ETA: usize>(
        rho: &[u8],
        nonce: u16,
        poly: &mut Polynomial
    ) {
        if ETA == 2 {
            Polynomial::fill_uniform_random_base::<
                66,
                POLY_UNIFORM_ETA_2_BUFLEN,
                STREAM256_BLOCKBYTES
            >(
                rho,
                nonce,
                poly,
                rejection_eta::<ETA>,
                |input| Box::new(shake_256_reader(input))
            );
        } else {
            Polynomial::fill_uniform_random_base::<
                66,
                POLY_UNIFORM_ETA_4_BUFLEN,
                STREAM256_BLOCKBYTES
            >(
                rho,
                nonce,
                poly,
                rejection_eta::<ETA>,
                |input| Box::new(shake_256_reader(input))
            );
        }
    }


    /// This function fills the coefficients of a polynomial with values sampled from [0, Q-1].
    ///
    /// # Arguments
    /// * `rho` -  An array containing 32 random bytes
    /// * `nonce` - A random integer
    /// * `poly` - A polynomial from which the coefficients will be converted into the NTT domain.
    pub fn fill_uniform_random(rho: &[u8], nonce: u16, poly: &mut Polynomial) {
        Polynomial::fill_uniform_random_base::<
            34,
            REJ_UNIFORM_BUFLEN,
            STREAM128_BLOCKBYTES
        >(
            rho,
            nonce,
            poly,
            rejection_uniform,
            |input| Box::new(shake_128_reader(input))
        );
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

impl IndexMut<usize> for Polynomial {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.coefficients[index]
    }
}

impl AddAssign<&Polynomial> for Polynomial {
    fn add_assign(&mut self, rhs: &Polynomial) {
        for i in 0..N {
            self[i] += rhs[i];
        }
    }
}

impl Mul<&Polynomial> for Polynomial {
    type Output = Self;
    fn mul(self, rhs: &Polynomial) -> Self::Output {
        let mut poly = Polynomial::default();

        for i in 0..N {
            poly[i] = montgomery_reduce(self[i] as i64 * rhs[i] as i64)
        }

        poly
    }
}