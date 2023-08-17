use crate::constants::{
    D, N, POLY_UNIFORM_ETA_2_BUFLEN, POLY_UNIFORM_ETA_4_BUFLEN, Q, REJ_UNIFORM_BUFLEN,
    SHAKE256_RATE, STREAM128_BLOCKBYTES, STREAM256_BLOCKBYTES, ZETAS,
};
use crate::conversion::u16_to_bytes_le;
use crate::helper::{shake_128_reader, shake_256_reader};
use crate::reduce::{caddq, montgomery_reduce, reduce32};
use crate::rejection::{rejection_eta, rejection_uniform};
use crate::rounding::{decompose, make_hint, power2round};
use sha3::digest::XofReader;
use std::ops::{AddAssign, Index, IndexMut, Mul, Sub};

const D_MINUS_ONE: u8 = D - 1;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Polynomial {
    pub(crate) coefficients: [i32; N],
}

impl Polynomial {
    pub fn zeroes() -> Self {
        Self {
            coefficients: [0; N],
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
        for i in 0..N / 4 {
            buf[5 * i + 0] = (poly[4 * i + 0] >> 0) as u8;
            buf[5 * i + 1] = (poly[4 * i + 0] >> 8) as u8 | (poly[4 * i + 1] << 2) as u8;
            buf[5 * i + 2] = (poly[4 * i + 1] >> 6) as u8 | (poly[4 * i + 2] << 4) as u8;
            buf[5 * i + 3] = (poly[4 * i + 2] >> 4) as u8 | (poly[4 * i + 3] << 6) as u8;
            buf[5 * i + 4] = (poly[4 * i + 3] >> 2) as u8;
        }
    }

    pub fn pack_eta<const ETA: usize>(buf: &mut [u8], poly: &Polynomial) {
        let mut arr = [0u8; 8];

        if ETA == 2 {
            for i in 0..N / 8 {
                arr[0] = (ETA as i32 - poly[8 * i + 0]) as u8;
                arr[1] = (ETA as i32 - poly[8 * i + 1]) as u8;
                arr[2] = (ETA as i32 - poly[8 * i + 2]) as u8;
                arr[3] = (ETA as i32 - poly[8 * i + 3]) as u8;
                arr[4] = (ETA as i32 - poly[8 * i + 4]) as u8;
                arr[5] = (ETA as i32 - poly[8 * i + 5]) as u8;
                arr[6] = (ETA as i32 - poly[8 * i + 6]) as u8;
                arr[7] = (ETA as i32 - poly[8 * i + 7]) as u8;

                buf[3 * i + 0] = (arr[0] >> 0) | (arr[1] << 3) | (arr[2] << 6);
                buf[3 * i + 1] = (arr[2] >> 2) | (arr[3] << 1) | (arr[4] << 4) | (arr[5] << 7);
                buf[3 * i + 2] = (arr[5] >> 1) | (arr[6] << 2) | (arr[7] << 5);
            }
        } else {
            for i in 0..N / 2 {
                arr[0] = (ETA as i32 - poly[2 * i + 0]) as u8;
                arr[1] = (ETA as i32 - poly[2 * i + 1]) as u8;
                buf[i] = arr[0] | (arr[1] << 4);
            }
        }
    }

    pub fn pack_t0(buf: &mut [u8], poly: &Polynomial) {
        let mut arr = [0i32; 8];
        const UPPER_BOUND: i32 = 1 << (D - 1);

        for i in 0..N / 8 {
            arr[0] = UPPER_BOUND - poly[8 * i + 0];
            arr[1] = UPPER_BOUND - poly[8 * i + 1];
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

    pub fn pack_w1<const GAMMA_TWO: usize>(buf: &mut [u8], poly: &Polynomial) {
        if GAMMA_TWO == ((Q - 1) / 88) as usize {
            for i in 0..N / 4 {
                buf[3 * i + 0] = poly[4 * i + 0] as u8;
                buf[3 * i + 0] |= (poly[4 * i + 1] << 6) as u8;
                buf[3 * i + 1] = (poly[4 * i + 1] >> 2) as u8;
                buf[3 * i + 1] |= (poly[4 * i + 2] << 4) as u8;
                buf[3 * i + 2] = (poly[4 * i + 2] >> 4) as u8;
                buf[3 * i + 2] |= (poly[4 * i + 3] << 2) as u8;
            }
        } else {
            for i in 0..N / 2 {
                buf[i] = (poly[2 * i + 0] | poly[2 * i + 1] << 4) as u8;
            }
        }
    }

    pub fn pack_z<const GAMMA_ONE: usize>(buf: &mut [u8], poly: &Polynomial) {
        let mut t = [0u32; 4];

        if GAMMA_ONE == 1 << 17 {
            for i in 0..N / 4 {
                t[0] = (GAMMA_ONE as i32 - poly[4 * i + 0]) as u32;
                t[1] = (GAMMA_ONE as i32 - poly[4 * i + 1]) as u32;
                t[2] = (GAMMA_ONE as i32 - poly[4 * i + 2]) as u32;
                t[3] = (GAMMA_ONE as i32 - poly[4 * i + 3]) as u32;

                buf[9 * i + 0] = t[0] as u8;
                buf[9 * i + 1] = (t[0] >> 8) as u8;
                buf[9 * i + 2] = (t[0] >> 16) as u8;
                buf[9 * i + 2] |= (t[1] << 2) as u8;
                buf[9 * i + 3] = (t[1] >> 6) as u8;
                buf[9 * i + 4] = (t[1] >> 14) as u8;
                buf[9 * i + 4] |= (t[2] << 4) as u8;
                buf[9 * i + 5] = (t[2] >> 4) as u8;
                buf[9 * i + 6] = (t[2] >> 12) as u8;
                buf[9 * i + 6] |= (t[3] << 6) as u8;
                buf[9 * i + 7] = (t[3] >> 2) as u8;
                buf[9 * i + 8] = (t[3] >> 10) as u8;
            }
        } else {
            for i in 0..N / 2 {
                t[0] = (GAMMA_ONE as i32 - poly[2 * i + 0]) as u32;
                t[1] = (GAMMA_ONE as i32 - poly[2 * i + 1]) as u32;

                buf[5 * i + 0] = t[0] as u8;
                buf[5 * i + 1] = (t[0] >> 8) as u8;
                buf[5 * i + 2] = (t[0] >> 16) as u8;
                buf[5 * i + 2] |= (t[1] << 4) as u8;
                buf[5 * i + 3] = (t[1] >> 4) as u8;
                buf[5 * i + 4] = (t[1] >> 12) as u8;
            }
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
                    j += 1;
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
        const XOF_INPUT_BUFFER_SIZE: usize,
        const ARR_SIZE: usize,
        const STREAM_BLOCK_BYTES: usize,
    >(
        poly: &mut Polynomial,
        rejection_fn: impl Fn(&mut Polynomial, usize, usize, &[u8], usize) -> usize,
        reader: &mut impl XofReader,
    ) {
        let mut buffer = [0u8; ARR_SIZE];
        let mut buffer_len = ARR_SIZE;
        reader.read(&mut buffer);

        let mut counter = rejection_fn(poly, N, 0, &buffer, buffer_len);

        while counter < N {
            let offset = buffer_len % 3;
            for i in 0..offset {
                buffer[i] = buffer[buffer_len - offset + 1];
            }
            reader.read(&mut buffer[offset..STREAM_BLOCK_BYTES]);
            buffer_len = STREAM_BLOCK_BYTES + offset;
            counter += rejection_fn(poly, N - counter, counter, &buffer, buffer_len);
        }
    }

    /// This function fills the coefficients of a polynomial with values sampled from
    /// [-ETA, ETA].
    ///
    /// # Arguments
    /// * `rho` -  An array containing 32 random bytes
    /// * `nonce` - A random integer
    /// * `poly` - A polynomial
    pub fn fill_uniform_eta<const ETA: usize>(rho: &[u8], nonce: u16, poly: &mut Polynomial) {
        if ETA == 2 {
            Polynomial::fill_uniform_random_base::<
                66,
                POLY_UNIFORM_ETA_2_BUFLEN,
                STREAM256_BLOCKBYTES,
            >(
                poly,
                rejection_eta::<ETA>,
                &mut shake_256_reader(&[rho, &u16_to_bytes_le(nonce)]),
            );
        } else {
            Polynomial::fill_uniform_random_base::<
                66,
                POLY_UNIFORM_ETA_4_BUFLEN,
                STREAM256_BLOCKBYTES,
            >(
                poly,
                rejection_eta::<ETA>,
                &mut shake_256_reader(&[rho, &u16_to_bytes_le(nonce)]),
            );
        }
    }

    /// This function fills the coefficients of a polynomial with values sampled from
    /// [-(GAMMA1 - 1), GAMMA1].
    ///
    /// # Arguments
    /// * `rho` -  An array containing 32 random bytes
    /// * `nonce` - A random integer
    /// * `poly` - A polynomial
    pub fn fill_uniform_gamma<const N: usize, const GAMMA_ONE: usize>(
        rho: &[u8],
        nonce: u16,
        poly: &mut Polynomial,
    ) {
        let mut buffer = [0u8; N];
        let bytes = u16_to_bytes_le(nonce);
        let mut reader = shake_256_reader(&[&rho, &bytes]);
        reader.read(&mut buffer);
        Polynomial::unpack_z::<GAMMA_ONE>(&buffer, poly);
    }

    /// This function fills the coefficients of a polynomial with values sampled from [0, Q-1].
    ///
    /// # Arguments
    /// * `rho` -  An array containing 32 random bytes
    /// * `nonce` - A random integer
    /// * `poly` - A polynomial
    pub fn fill_uniform_random(rho: &[u8], nonce: u16, poly: &mut Polynomial) {
        Polynomial::fill_uniform_random_base::<34, REJ_UNIFORM_BUFLEN, STREAM128_BLOCKBYTES>(
            poly,
            rejection_uniform,
            &mut shake_128_reader(&[rho, &u16_to_bytes_le(nonce)]),
        );
    }

    pub fn unpack_eta<const ETA: usize>(poly: &mut Polynomial, buf: &[u8]) {
        if ETA == 2 {
            for i in 0..N / 8 {
                poly[8 * i + 0] = ((buf[3 * i + 0] >> 0) & 7) as i32;
                poly[8 * i + 1] = ((buf[3 * i + 0] >> 3) & 7) as i32;
                poly[8 * i + 2] = (((buf[3 * i + 0] >> 6) | (buf[3 * i + 1] << 2)) & 7) as i32;
                poly[8 * i + 3] = ((buf[3 * i + 1] >> 1) & 7) as i32;
                poly[8 * i + 4] = ((buf[3 * i + 1] >> 4) & 7) as i32;
                poly[8 * i + 5] = (((buf[3 * i + 1] >> 7) | (buf[3 * i + 2] << 1)) & 7) as i32;
                poly[8 * i + 6] = ((buf[3 * i + 2] >> 2) & 7) as i32;
                poly[8 * i + 7] = ((buf[3 * i + 2] >> 5) & 7) as i32;

                poly[8 * i + 0] = ETA as i32 - poly[8 * i + 0];
                poly[8 * i + 1] = ETA as i32 - poly[8 * i + 1];
                poly[8 * i + 2] = ETA as i32 - poly[8 * i + 2];
                poly[8 * i + 3] = ETA as i32 - poly[8 * i + 3];
                poly[8 * i + 4] = ETA as i32 - poly[8 * i + 4];
                poly[8 * i + 5] = ETA as i32 - poly[8 * i + 5];
                poly[8 * i + 6] = ETA as i32 - poly[8 * i + 6];
                poly[8 * i + 7] = ETA as i32 - poly[8 * i + 7];
            }
        } else {
            for i in 0..N / 2 {
                poly[2 * i + 0] = (buf[i] & 0x0F) as i32;
                poly[2 * i + 1] = (buf[i] >> 4) as i32;
                poly[2 * i + 0] = ETA as i32 - poly[2 * i + 0];
                poly[2 * i + 1] = ETA as i32 - poly[2 * i + 1];
            }
        }
    }

    pub fn unpack_t0(poly: &mut Polynomial, buf: &[u8]) {
        for i in 0..N / 8 {
            poly[8 * i + 0] = buf[13 * i + 0] as i32;
            poly[8 * i + 0] |= (buf[13 * i + 1] as i32) << 8;
            poly[8 * i + 0] &= 0x1FFF;

            poly[8 * i + 1] = (buf[13 * i + 1] >> 5) as i32;
            poly[8 * i + 1] |= (buf[13 * i + 2] as i32) << 3;
            poly[8 * i + 1] |= (buf[13 * i + 3] as i32) << 11;
            poly[8 * i + 1] &= 0x1FFF;

            poly[8 * i + 2] = (buf[13 * i + 3] >> 2) as i32;
            poly[8 * i + 2] |= (buf[13 * i + 4] as i32) << 6;
            poly[8 * i + 2] &= 0x1FFF;

            poly[8 * i + 3] = (buf[13 * i + 4] as u32 >> 7) as i32;
            poly[8 * i + 3] |= (buf[13 * i + 5] as i32) << 1;
            poly[8 * i + 3] |= (buf[13 * i + 6] as i32) << 9;
            poly[8 * i + 3] &= 0x1FFF;

            poly[8 * i + 4] = (buf[13 * i + 6] >> 4) as i32;
            poly[8 * i + 4] |= (buf[13 * i + 7] as i32) << 4;
            poly[8 * i + 4] |= (buf[13 * i + 8] as i32) << 12;
            poly[8 * i + 4] &= 0x1FFF;

            poly[8 * i + 5] = (buf[13 * i + 8] >> 1) as i32;
            poly[8 * i + 5] |= (buf[13 * i + 9] as i32) << 7;
            poly[8 * i + 5] &= 0x1FFF;

            poly[8 * i + 6] = (buf[13 * i + 9] >> 6) as i32;
            poly[8 * i + 6] |= (buf[13 * i + 10] as i32) << 2;
            poly[8 * i + 6] |= (buf[13 * i + 11] as i32) << 10;
            poly[8 * i + 6] &= 0x1FFF;

            poly[8 * i + 7] = (buf[13 * i + 11] >> 3) as i32;
            poly[8 * i + 7] |= (buf[13 * i + 12] as i32) << 5;
            poly[8 * i + 7] &= 0x1FFF;

            poly[8 * i + 0] = (1 << (D_MINUS_ONE)) - poly[8 * i + 0];
            poly[8 * i + 1] = (1 << (D_MINUS_ONE)) - poly[8 * i + 1];
            poly[8 * i + 2] = (1 << (D_MINUS_ONE)) - poly[8 * i + 2];
            poly[8 * i + 3] = (1 << (D_MINUS_ONE)) - poly[8 * i + 3];
            poly[8 * i + 4] = (1 << (D_MINUS_ONE)) - poly[8 * i + 4];
            poly[8 * i + 5] = (1 << (D_MINUS_ONE)) - poly[8 * i + 5];
            poly[8 * i + 6] = (1 << (D_MINUS_ONE)) - poly[8 * i + 6];
            poly[8 * i + 7] = (1 << (D_MINUS_ONE)) - poly[8 * i + 7];
        }
    }

    fn unpack_z<const GAMMA_ONE: usize>(buf: &[u8], poly: &mut Polynomial) {
        if GAMMA_ONE == 1 << 17 {
            for i in 0..N / 4 {
                poly[4 * i + 0] = buf[9 * i + 0] as i32;
                poly[4 * i + 0] |= (buf[9 * i + 1] as i32) << 8;
                poly[4 * i + 0] |= (buf[9 * i + 2] as i32) << 16;
                poly[4 * i + 0] &= 0x3FFFF;

                poly[4 * i + 1] = buf[9 * i + 2] as i32 >> 2;
                poly[4 * i + 1] |= (buf[9 * i + 3] as i32) << 6;
                poly[4 * i + 1] |= (buf[9 * i + 4] as i32) << 14;
                poly[4 * i + 1] &= 0x3FFFF;

                poly[4 * i + 2] = buf[9 * i + 4] as i32 >> 4;
                poly[4 * i + 2] |= (buf[9 * i + 5] as i32) << 4;
                poly[4 * i + 2] |= (buf[9 * i + 6] as i32) << 12;
                poly[4 * i + 2] &= 0x3FFFF;

                poly[4 * i + 3] = (buf[9 * i + 6] as i32) >> 6;
                poly[4 * i + 3] |= (buf[9 * i + 7] as i32) << 2;
                poly[4 * i + 3] |= (buf[9 * i + 8] as i32) << 10;
                poly[4 * i + 3] &= 0x3FFFF;

                poly[4 * i + 0] = GAMMA_ONE as i32 - poly[4 * i + 0];
                poly[4 * i + 1] = GAMMA_ONE as i32 - poly[4 * i + 1];
                poly[4 * i + 2] = GAMMA_ONE as i32 - poly[4 * i + 2];
                poly[4 * i + 3] = GAMMA_ONE as i32 - poly[4 * i + 3];
            }
        } else {
            for i in 0..N / 2 {
                poly[2 * i + 0] = buf[5 * i + 0] as i32;
                poly[2 * i + 0] |= (buf[5 * i + 1] as i32) << 8;
                poly[2 * i + 0] |= (buf[5 * i + 2] as i32) << 16;
                poly[2 * i + 0] &= 0xFFFFF;

                poly[2 * i + 1] = (buf[5 * i + 2] as i32) >> 4;
                poly[2 * i + 1] |= (buf[5 * i + 3] as i32) << 4;
                poly[2 * i + 1] |= (buf[5 * i + 4] as i32) << 12;
                poly[2 * i + 0] &= 0xFFFFF;

                poly[2 * i + 0] = GAMMA_ONE as i32 - poly[2 * i + 0];
                poly[2 * i + 1] = GAMMA_ONE as i32 - poly[2 * i + 1];
            }
        }
    }

    pub fn decompose<const GAMMA_TWO: usize>(
        r_zero: &mut Polynomial,
        r_one: &mut Polynomial,
        poly: &Polynomial,
    ) {
        for i in 0..N {
            let (a_zero, a_one) = decompose::<GAMMA_TWO>(poly[i]);
            r_zero[i] = a_zero;
            r_one[i] = a_one
        }
    }

    /// This function returns a Boolean indicating whether the infinite norm of the polynomial
    /// is greater than or equal to a specific bound.
    ///
    /// # Arguments:
    /// * `bound` - The upper bound
    pub fn check_infinite_norm(&self, bound: usize) -> bool {
        if bound > ((Q - 1) / 8) as usize {
            return true;
        }

        for coeff in self.coefficients {
            let coeff_abs = coeff.unsigned_abs() as usize;

            if coeff_abs >= bound {
                return true;
            }
        }

        false
    }

    pub fn make_hint<const GAMMA_TWO: usize>(
        output_poly: &mut Polynomial,
        poly_one: &Polynomial,
        poly_zero: &Polynomial,
    ) -> usize {
        let mut number_ones = 0;

        for i in 0..N {
            output_poly[i] = make_hint::<GAMMA_TWO>(poly_zero[i], poly_one[i]);
            number_ones += output_poly[i]
        }

        number_ones as usize
    }

    pub fn sample_in_ball<const TAU: usize>(c_tilde: &[u8]) -> Polynomial {
        let mut poly = Polynomial::default();
        let mut signs = 0u64;
        let mut buffer = [0u8; SHAKE256_RATE];

        let mut reader = shake_256_reader(&[c_tilde]);
        reader.read(&mut buffer);

        for i in 0..8 {
            signs |= (buffer[i] as u64) << 8 * i;
        }

        let mut pos: usize = 8;
        let mut byte_value;

        for i in N - TAU..N {
            loop {
                if pos >= SHAKE256_RATE {
                    reader.read(&mut buffer);
                    pos = 0;
                }

                byte_value = buffer[pos];
                pos += 1;

                if byte_value <= i as u8 {
                    break;
                }
            }

            poly[i] = poly[byte_value as usize];
            poly[byte_value as usize] = 1i32 - 2 * (signs & 1) as i32;
            signs >>= 1;
        }

        poly
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

impl Sub<&Polynomial> for Polynomial {
    type Output = Polynomial;

    fn sub(self, rhs: &Polynomial) -> Self::Output {
        let mut poly = Polynomial::default();

        for i in 0..N {
            poly[i] = self[i] - rhs[i]
        }

        poly
    }
}
