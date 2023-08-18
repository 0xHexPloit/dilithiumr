use crate::constants::{
    D, N, POLY_UNIFORM_ETA_2_BUFLEN, POLY_UNIFORM_ETA_4_BUFLEN, Q, REJ_UNIFORM_BUFLEN,
    SHAKE256_RATE, STREAM128_BLOCKBYTES, STREAM256_BLOCKBYTES, ZETAS,
};
use crate::conversion::u16_to_bytes_le;
use crate::helper::{shake_128_reader, shake_256_reader};
use crate::reduce::{caddq, montgomery_reduce, reduce32};
use crate::rejection::{rejection_eta, rejection_uniform};
use crate::rounding::{decompose, make_hint, power2round, use_hint};
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

    pub fn unpack_z<const GAMMA_ONE: usize>(buf: &[u8], poly: &mut Polynomial) {
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

    pub fn unpack_t1(buf: &[u8], poly: &mut Polynomial) {
        for i in 0..N / 4 {
            poly[4 * i + 0] =
                (((buf[5 * i + 0] as u32 >> 0) | ((buf[5 * i + 1] as u32) << 8)) & 0x3FF) as i32;
            poly[4 * i + 1] =
                (((buf[5 * i + 1] as u32 >> 2) | ((buf[5 * i + 2] as u32) << 6)) & 0x3FF) as i32;
            poly[4 * i + 2] =
                (((buf[5 * i + 2] as u32 >> 4) | ((buf[5 * i + 3] as u32) << 4)) & 0x3FF) as i32;
            poly[4 * i + 3] =
                (((buf[5 * i + 3] as u32 >> 6) | ((buf[5 * i + 4] as u32) << 2)) & 0x3FF) as i32;
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

    /// This function fills a hint polynomial
    ///
    /// # Arguments
    /// * `poly` - The hint polynomial to be filled
    /// * `poly_zero` - The low part of an input polynomial
    /// * `poly_one` - The high part of an input polynomial
    pub fn make_hint<const GAMMA_TWO: usize>(
        poly: &mut Polynomial,
        poly_zero: &Polynomial,
        poly_one: &Polynomial,
    ) -> usize {
        let mut number_ones = 0;

        for i in 0..N {
            poly[i] = make_hint::<GAMMA_TWO>(poly_zero[i], poly_one[i]);
            number_ones += poly[i]
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

    /// This functions multiplies polynomial by 2^D without modular reduction. The latter assumes
    /// input coefficients to be less than 2^{31-D} in absolute value.
    pub fn shiftl(&mut self) {
        for i in 0..N {
            self[i] <<= D;
        }
    }

    /// This function uses the hint polynomial to correct the high bits of input polynomial and
    /// fills the output polynomial with the correction.
    ///
    /// # Arguments
    ///
    /// * `h` - A hint polynomial
    /// * `vec` - A polynomial for which we should correct the high bits.
    /// * `poly` - A poly that will be filled with the correction
    pub fn use_hint<const GAMMA_TWO: usize>(
        h: &Polynomial,
        vec: &Polynomial,
        poly: &mut Polynomial,
    ) {
        for i in 0..N {
            poly[i] = use_hint::<GAMMA_TWO>(vec[i], h[i] as u32);
        }
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

#[cfg(test)]
mod tests {
    use crate::algebra::poly::Polynomial;
    use crate::algebra::vec::PolyVec;
    use crate::constants::POLY_T0_PACKED_BYTES;

    #[test]
    fn test_should_unpack_t0_correctly() {
        let t0_packed_bytes = [
            246, 244, 86, 109, 35, 71, 226, 135, 131, 140, 75, 207, 56, 56, 92, 66, 119, 128, 105,
            142, 11, 254, 171, 32, 213, 225, 172, 220, 226, 70, 232, 26, 105, 76, 245, 178, 69,
            254, 186, 31, 139, 176, 69, 159, 198, 110, 186, 43, 139, 66, 13, 218, 109, 83, 192, 84,
            91, 51, 123, 240, 250, 49, 230, 230, 22, 49, 114, 134, 23, 153, 162, 127, 150, 151, 30,
            119, 230, 176, 217, 110, 184, 100, 52, 196, 185, 188, 74, 134, 56, 108, 53, 193, 169,
            28, 69, 87, 23, 205, 12, 40, 10, 82, 186, 47, 54, 63, 102, 31, 168, 9, 227, 152, 24,
            24, 205, 98, 188, 205, 119, 59, 108, 240, 251, 11, 23, 196, 243, 37, 143, 149, 53, 127,
            179, 50, 245, 144, 150, 143, 19, 94, 174, 37, 13, 186, 255, 232, 157, 195, 103, 80, 91,
            55, 83, 75, 120, 119, 70, 189, 60, 161, 175, 155, 115, 243, 163, 173, 215, 39, 206, 2,
            211, 240, 112, 105, 98, 68, 119, 113, 198, 107, 152, 158, 211, 149, 84, 151, 221, 144,
            176, 80, 155, 130, 182, 120, 239, 170, 144, 145, 252, 193, 200, 87, 29, 8, 7, 88, 110,
            70, 233, 216, 69, 179, 140, 107, 245, 117, 201, 3, 116, 201, 3, 203, 108, 249, 149,
            255, 57, 34, 91, 104, 139, 27, 139, 238, 54, 254, 48, 36, 242, 136, 118, 211, 15, 248,
            123, 64, 239, 133, 93, 150, 27, 53, 83, 140, 18, 101, 6, 243, 33, 206, 204, 234, 41, 3,
            77, 23, 47, 63, 192, 2, 249, 243, 36, 241, 225, 245, 37, 179, 95, 12, 150, 61, 209,
            166, 33, 161, 65, 103, 163, 74, 254, 158, 204, 76, 62, 150, 229, 249, 8, 151, 183, 197,
            45, 9, 193, 149, 113, 41, 12, 198, 26, 25, 67, 194, 53, 101, 19, 76, 115, 178, 58, 71,
            241, 96, 36, 222, 217, 103, 167, 112, 80, 200, 39, 33, 104, 50, 96, 80, 33, 214, 162,
            164, 139, 100, 90, 40, 168, 170, 236, 158, 2, 48, 226, 232, 72, 35, 197, 228, 209, 164,
            237, 134, 164, 68, 187, 153, 24, 148, 2, 130, 211, 66, 141, 170, 176, 174, 111, 73,
            214, 16, 63, 197, 84, 210, 95, 135, 206, 219, 189, 84, 242, 38, 120, 32, 213, 126, 172,
            253, 75, 169, 76, 207, 92, 173, 235, 206, 59, 94, 197, 59, 208, 123, 195, 81, 231, 135,
            150, 49, 180, 111, 248, 93, 145, 219, 132, 142, 158, 94, 198, 61, 242, 177, 149, 234,
            213, 21, 153, 77, 79, 147, 232, 65, 32, 232, 247, 10, 225, 167, 63, 38, 149, 122, 130,
            31, 235, 207, 49, 184, 180, 234, 223, 41, 145, 59, 9, 81, 76, 65, 170, 119, 218, 17,
            12, 173, 119, 71, 204, 22, 238, 39, 111, 209, 215, 117, 216, 82, 70, 24, 253, 184, 79,
            171, 202, 134, 118, 118, 34, 114, 164, 176, 14, 105, 23, 29, 71, 75, 102, 35, 156, 43,
            165, 17, 201, 74, 183, 5, 2, 236, 53, 0, 162, 255, 73, 96, 54, 35, 212, 177, 192, 19,
            38, 189, 20, 216, 135, 224, 76, 187, 12, 116, 79, 148, 182, 155, 138, 150, 27, 197, 88,
            112, 143, 184, 40, 141, 209, 174, 1, 43, 17, 12, 134, 170, 63, 187, 80, 42, 90, 136,
            58, 92, 23, 138, 237, 149, 247, 231, 126, 30, 32, 210, 173, 249, 34, 202, 136, 250,
            230, 177, 129, 27, 144, 230, 96, 27, 76, 170, 52, 240, 120, 180, 185, 71, 2, 83, 111,
            33, 149, 147, 246, 62, 206, 218, 146, 220, 64, 240, 157, 238, 82, 205, 224, 85, 243,
            202, 176, 126, 232, 146, 67, 229, 250, 185, 196, 116, 159, 143, 75, 10, 203, 17, 163,
            116, 80, 44, 159, 141, 138, 148, 225, 237, 16, 186, 106, 205, 45, 98, 123, 65, 103,
            247, 188, 183, 204, 169, 222, 53, 177, 116, 45, 176, 204, 126, 252, 116, 114, 44, 59,
            201, 51, 50, 49, 162, 128, 89, 153, 126, 122, 137, 16, 251, 227, 153, 206, 104, 111,
            179, 107, 127, 65, 86, 246, 149, 179, 30, 171, 109, 95, 233, 47, 86, 172, 43, 146, 214,
            238, 140, 181, 174, 173, 89, 193, 242, 96, 234, 209, 44, 171, 122, 59, 175, 53, 176,
            173, 124, 235, 213, 44, 147, 11, 182, 99, 39, 98, 48, 194, 51, 96, 81, 10, 146, 210,
            20, 64, 152, 228, 110, 140, 134, 239, 38, 131, 174, 212, 111, 4, 157, 235, 195, 139,
            191, 7, 246, 125, 114, 125, 204, 61, 95, 219, 55, 20, 156, 142, 38, 245, 71, 170, 33,
            140, 161, 12, 204, 36, 231, 142, 105, 120, 34, 97, 211, 226, 67, 193, 53, 179, 105, 78,
            130, 31, 134, 195, 118, 136, 101, 167, 227, 100, 82, 119, 99, 35, 59, 46, 175, 92, 200,
            101, 76, 200, 254, 240, 254, 184, 154, 81, 147, 234, 216, 42, 155, 45, 214, 198, 69, 8,
            255, 147, 192, 108, 86, 202, 176, 27, 144, 203, 194, 49, 187, 80, 20, 54, 146, 208, 60,
            237, 151, 22, 205, 6, 89, 175, 117, 250, 87, 78, 224, 100, 92, 247, 176, 50, 6, 49,
            175, 61, 8, 205, 30, 231, 64, 233, 115, 206, 112, 71, 60, 162, 230, 61, 212, 109, 239,
            232, 170, 42, 145, 73, 71, 122, 231, 186, 227, 207, 163, 228, 172, 113, 149, 140, 53,
            145, 55, 86, 51, 167, 129, 7, 106, 244, 57, 136, 212, 50, 243, 166, 101, 181, 107, 27,
            107, 166, 161, 4, 102, 77, 66, 2, 9, 52, 194, 209, 46, 204, 44, 131, 226, 6, 102, 233,
            193, 61, 100, 99, 181, 167, 229, 225, 217, 115, 101, 155, 115, 235, 120, 236, 221, 181,
            234, 111, 64, 85, 145, 55, 251, 157, 69, 211, 108, 91, 96, 213, 185, 160, 38, 68, 107,
            151, 13, 104, 34, 228, 10, 120, 250, 32, 229, 138, 105, 231, 231, 80, 236, 191, 185,
            53, 2, 171, 141, 107, 178, 154, 28, 143, 183, 140, 52, 82, 81, 11, 205, 136, 73, 252,
            77, 200, 92, 164, 166, 217, 60, 6, 189, 96, 231, 143, 56, 204, 46, 9, 83, 181, 31, 133,
            179, 9, 98, 62, 230, 242, 151, 122, 90, 75, 100, 53, 181, 223, 124, 25, 50, 157, 47,
            120, 198, 103, 108, 112, 48, 224, 59, 160, 253, 201, 119, 198, 91, 176, 20, 199, 202,
            8, 131, 112, 71, 68, 29, 126, 241, 186, 214, 22, 28, 112, 117, 41, 94, 221, 239, 25,
            115, 117, 249, 246, 88, 203, 225, 220, 255, 97, 229, 132, 217, 139, 152, 39, 84, 162,
            100, 153, 57, 79, 250, 24, 87, 39, 226, 58, 167, 180, 253, 189, 123, 102, 165, 58, 170,
            19, 221, 30, 242, 89, 17, 138, 235, 190, 248, 95, 227, 172, 220, 103, 228, 220, 120,
            118, 238, 62, 71, 210, 252, 228, 0, 198, 127, 168, 112, 168, 167, 167, 99, 190, 127,
            165, 247, 178, 139, 239, 224, 74, 198, 62, 242, 184, 214, 138, 11, 229, 160, 90, 153,
            86, 76, 164, 243, 116, 188, 202, 170, 12, 57, 154, 240, 58, 235, 128, 218, 136, 35, 27,
            142, 83, 7, 255, 177, 96, 152, 91, 184, 19, 136, 10, 131, 50, 195, 137, 127, 69, 252,
            90, 49, 244, 191, 236, 151, 213, 111, 201, 50, 32, 121, 156, 125, 52, 191, 54, 139, 8,
            144, 117, 19, 255, 82, 48, 42, 6, 62, 62, 202, 78, 133, 56, 217, 135, 195, 189, 5, 114,
            205, 35, 85, 74, 173, 222, 170, 27, 239, 6, 98, 120, 128, 195, 129, 23, 45, 214, 18,
            199, 3, 125, 207, 174, 168, 133, 111, 163, 103, 247, 98, 206, 84, 25, 169, 150, 162,
            129, 62, 107, 156, 119, 143, 81, 222, 226, 252, 234, 174, 106, 53, 64, 137, 236, 247,
            110, 164, 89, 227, 8, 97, 126, 84, 202, 153, 64, 95, 178, 211, 57, 166, 177, 238, 127,
            239, 163, 36, 167, 216, 144, 124, 117, 219, 205, 9, 218, 74, 113, 225, 195, 54, 156,
            74, 93, 229, 173, 169, 255, 73, 94, 169, 209, 213, 107, 242, 1, 66, 182, 73, 184, 221,
            218, 157, 6, 45, 228, 107, 54, 152, 64, 238, 64, 121, 14, 12, 164, 5, 227, 124, 232,
            252, 36, 59, 120, 69, 63, 159, 128, 239, 95, 38, 67, 155, 213, 229, 176, 84, 97, 134,
            202, 57, 64, 215, 139, 55, 87, 100, 248, 158, 57, 209, 117, 17, 206, 242, 183, 26, 242,
            190, 144, 231, 112, 17, 132, 242, 15, 59, 150, 232, 30, 76, 95, 76, 188, 101, 187, 219,
            203, 207, 104, 204, 82, 120, 53, 164, 211, 230, 253, 214, 61, 230, 214, 14, 177, 248,
            107, 124, 130, 16, 255, 182, 181, 76, 203, 160, 29, 218, 78, 146, 190, 138, 43, 115,
            57, 182, 76, 177, 204, 84, 246, 232, 83, 218, 235, 35, 189, 252, 177, 164, 26, 83, 242,
            12, 122, 216, 6, 198, 25, 24, 18, 50, 83, 217, 99, 103, 95, 246, 144, 109, 179, 102,
            171, 139, 186, 171, 209, 247, 129, 192, 233, 77, 179, 78, 23, 110, 36, 70, 11, 73, 33,
            57, 238, 246, 254, 3, 113, 12, 38, 129, 165, 132, 58, 180, 84, 94, 81, 6, 255, 138,
            158, 179, 164, 69, 41, 149, 220, 144, 107, 254, 165, 42, 233, 231, 73, 15, 5, 109, 105,
            87, 60, 81, 187, 81, 203, 240, 121, 32, 140, 146, 181, 179, 40, 84, 95, 121, 213, 151,
            51, 249, 116, 85, 142, 66, 25, 24, 199, 8, 133, 21, 142, 187, 144, 230, 41, 253, 120,
            45, 143, 149, 46, 13, 7, 10, 124, 20, 88, 244, 161, 39, 49, 45, 219, 87, 163, 77, 102,
            12, 139, 217, 130, 67, 106, 71, 196, 57, 119, 199,
        ];
        let mut t_0 = PolyVec::<4>::default();

        for i in 0..4 {
            Polynomial::unpack_t0(&mut t_0[i], &t0_packed_bytes[i * POLY_T0_PACKED_BYTES..]);
        }

        let expected_t0 = PolyVec::<4> {
            data: [
                Polynomial {
                    coefficients: [
                        -1270, 1353, 1829, 2930, -2174, 2495, -3374, 2279, -3128, -2578, 4067,
                        -3283, 3912, -1535, -1154, -3130, -3244, -1814, -2577, -565, -1222, -2426,
                        -2326, -1887, 1249, 636, 2095, -3469, -2982, 2667, -1290, -2881, -877,
                        2558, -1749, -1638, 249, -2301, -2968, 3364, -561, -3123, 2491, -3909,
                        -2407, 181, -2524, -1564, 295, 2621, 743, -904, 1077, 3291, -226, 2387,
                        1599, 1819, -1489, -2606, 3892, 2796, 1720, 2569, -3894, -2865, 1529, 2541,
                        1650, 1012, 1228, -1932, -1997, 3621, -3099, -2039, 3728, -2530, -3223,
                        -689, -3893, -1435, -3404, 735, -2297, 247, -1721, 3676, -4026, 185, -231,
                        3889, -1461, 1637, 3795, 273, -3398, 1563, 1048, 2249, -3895, -1745, -3934,
                        -2500, -770, 2170, -2652, 1852, -1908, 3272, 3665, -979, -1491, -2724,
                        -1893, 3807, -1291, 3763, 3366, -3567, -170, 2932, -127, 111, 3627, 3196,
                        -2400, 1843, -2281, -2606, 3284, 1321, -1887, 3612, 2608, 3975, 821, 53,
                        -4069, 2957, 2638, 2636, 914, -3537, -3638, 3705, -3209, 751, -3383, -3079,
                        3601, -3560, -3461, -3250, 698, -2214, -296, 3278, 2100, -2500, 1332,
                        -2383, -832, -3630, 3086, 3744, 28, 2914, 3599, 81, 823, -2239, -2400,
                        1890, 2405, -1060, 2239, -1307, -3986, -2365, 2868, 1249, 2154, 3809,
                        -1943, 467, 3509, 1150, -1817, 2540, 1256, 3293, 3517, 1618, 2855, 2408,
                        1241, 3171, 3131, 2932, -2526, -2878, -3113, -160, -636, -1040, 3895, 1524,
                        -1569, 2794, 3351, -1225, 3451, -1364, -2994, 4013, 3536, 2233, 1838, 1654,
                        738, -1746, -539, 1900, -2491, 3900, 3931, 2300, -1069, -1350, -2754, 523,
                        -1609, -2182, -335, 2903, -1533, 2237, -1903, 1385, 2318, 3135, -1352,
                        -2301, -4058, -1189, -3378, 1127, 1107, -3703, -1934, -1930, -3331, 3651,
                        -3399, -252,
                    ],
                },
                Polynomial {
                    coefficients: [
                        -406, -3489, -3611, 3397, 583, 2238, -2682, -2251, -573, 625, -2725, 1109,
                        -2449, 2138, 3507, 1987, 2016, -1983, -2114, -3919, -611, -3402, -3593,
                        -3427, -463, 2623, -2733, -959, -2322, 2915, -324, 2007, -1962, 301, 3324,
                        166, 2953, 1178, -4024, 540, -2001, 3154, -1206, -140, 47, 2084, 1363,
                        -217, -1654, -275, 1764, -3425, -1680, 373, 740, 823, -3107, 1700, 2967,
                        -1426, -2932, 3838, -1968, 4090, -4002, 3505, 616, 1978, 1251, 1568, -1176,
                        3433, 2088, 2300, 301, 2023, 2825, -2890, 1426, -721, 2789, 3386, 3108,
                        -369, -2258, -1896, 1018, 3547, 2548, -3412, 305, -1185, 2654, -3396,
                        -3440, -322, -1517, -4028, -4025, 4036, -3362, -3286, 1909, -281, 2310,
                        625, 2336, 736, -1550, 2547, -681, -3590, -1144, -3533, 3951, -3750, -534,
                        1590, -3034, -2503, -730, 2332, -3088, -3387, -1326, -102, 681, -2398,
                        -3760, -1859, -228, -1482, 1121, -2658, -3709, 1679, 1270, -2190, -3368,
                        -2208, -2546, 2746, 2478, -3516, -2576, 1195, 1165, -1732, -1047, -2995,
                        -3827, -2454, -3753, 1618, -3372, 4006, 821, -3647, 1581, 2674, 1733, -414,
                        948, 3772, -1432, -3916, 2583, 3567, 3077, -1231, -2611, 2338, -1723, 3905,
                        -2393, -702, -3763, 680, -2011, -4050, 2718, -1494, -2632, -3546, -1420,
                        651, -1643, 2686, 2545, 1803, 845, 171, 197, 3667, 1172, -1785, 674, 1642,
                        -2094, 906, 3545, -387, 784, 3392, 3931, 1719, 4013, -776, 284, -1123,
                        -3041, 2483, 1304, -2026, -1041, -3443, 1085, -3580, -3457, 2821, 2089,
                        -3814, 644, 2309, -3092, -1140, -3401, -1167, 3558, -198, -50, 2919, 281,
                        3252, 1890, 2366, -3629, 3935, 809, 714, 3506, -252, -225, -237, -1624,
                        -467, 1645, 278, 3229, -473, 1077, -185, 2468, 3034, 3077, -4062,
                    ],
                },
                Polynomial {
                    coefficients: [
                        -2744, -2700, -2724, -1457, -2482, 1258, -1819, 3832, -1023, 2556, -1435,
                        3692, 3653, 2616, 2293, -1894, -1104, -432, -1060, -2681, 1666, 2421, 3045,
                        -1515, -2677, -703, -2067, -2249, 139, -2392, 3048, -1510, 1987, -1640,
                        -2503, -641, 2242, -2151, -285, -1095, -3558, 351, -3035, -1489, -682,
                        2872, 1763, -3311, 3142, -3711, -2344, 3239, 1705, -2758, -3652, 1338,
                        2253, -3085, -2689, -1000, 1917, -2410, -3020, 844, 1099, -2267, 1638,
                        1725, -1632, 3802, 3063, 2431, -450, 3722, 1229, 2810, 3986, -1203, -1799,
                        889, -1379, 707, -2169, 2125, -1623, -2509, 3155, -3471, -1501, -3925, -27,
                        3414, -889, 259, 746, 614, 4005, 341, 2002, 1971, -1716, 2357, 1632, -3204,
                        -2058, 2093, -2376, -789, -3702, 1933, -4017, -1847, 3531, 680, -2787,
                        -1380, -457, -3015, -562, 1466, 1199, 2456, -610, -3064, 892, -558, 2406,
                        2149, -3334, -2821, 3079, -2161, -748, 1660, -3797, -163, 1613, -784,
                        -2447, 27, 2135, 2643, -1425, -1702, -3295, -203, 2228, -95, -3175, -1587,
                        3647, -3078, 4037, 19, -3570, -1932, 1275, 3190, 3285, -97, 2192, 1502,
                        -3975, -1506, 661, 501, -1472, 2770, -3422, 130, -3270, -746, 145, 2644,
                        -903, -4091, 2719, 985, 3338, 207, 2750, -593, 2459, 1561, -2298, -2744,
                        -2185, 395, -2890, -3838, -2542, -1196, 1478, 1891, 2121, -996, 3819,
                        -1477, 3333, 1025, 797, -3813, -2329, -441, 2201, -3959, 1764, -3994, 3868,
                        -3632, 1505, -225, -2682, -467, -3833, -1199, -759, -3165, -2107, 875,
                        3092, -3193, 1190, 3727, 3867, 1323, -1446, 1896, 198, -3642, 1238, 3691,
                        -2617, -1924, -2766, -1281, -2189, 623, 456, 3862, -511, 3323, -1766, 2192,
                        1919, 3707, 822, -312, 2689, -2018, 938, -4072, -3787, 1333, 2625, 2471,
                    ],
                },
                Polynomial {
                    coefficients: [
                        -2336, 797, 737, 642, 1869, 2044, 554, -4066, -82, -337, 127, -1148, -1260,
                        -3138, -3940, -2160, 2627, 1136, 1805, -1194, 1324, -1391, -3182, 3875,
                        -2146, -3075, 3984, -2607, 670, 3191, -1039, -2543, 1874, -3117, 1829, 305,
                        2513, 1433, 2971, -725, 3678, -2548, 2278, -3823, 2792, -367, 1037, -1501,
                        -1386, 1535, -2850, -3567, -2630, -428, 3037, 52, 1452, 2866, -2000, 2204,
                        3171, -2259, -4026, -3567, 2909, 2759, 3018, 1287, -3511, 2842, 1176, 471,
                        3103, 3658, -679, 1350, -2782, -4052, -2343, -1323, -1489, -862, 3972, 892,
                        2917, 292, -1899, 3885, 3027, -863, 2547, -3201, -1038, 2244, -48, 3916,
                        -3299, 2237, 1729, -118, -1111, 97, -3586, 1027, 3290, 806, -2421, 1695,
                        2539, 2749, 3865, -2792, -1931, 3399, -3609, -829, -3347, 1862, 1224,
                        -1790, -538, 2569, -2532, 3359, 1983, 2055, -2284, -3346, 994, 3334, 237,
                        -1739, -3515, 2075, -419, 1447, -1400, -3361, -2484, 517, 3107, 1165, 3013,
                        -3862, -3179, 3053, -4036, 1171, -1227, -101, 1930, 1573, -3730, -3157,
                        -3274, 910, -1227, 2472, -2387, -3358, -2643, -3934, 184, 3079, 1461, 1651,
                        -969, 191, 2344, 464, 2554, 3036, -1331, -492, -3485, -3787, 624, -1435,
                        1319, -1303, -2747, -3048, 3577, -3384, -845, -2677, -2949, 952, -180,
                        -164, -2276, -3805, 3074, 3192, 1661, 1278, 1974, -2589, -2386, 1493,
                        -3846, -1111, 793, 1207, -660, 438, 445, -4045, 1371, -3913, -633, 1506,
                        -1744, 1100, 2831, -1898, 1199, 122, 2018, 2792, -2905, -1113, -3408, 213,
                        -2005, 1636, -3390, -3242, -1064, 1012, 3300, -161, 491, 2596, -2468,
                        -2643, -1935, 2154, -2646, 3675, 1529, 3104, -1541, 3096, -634, -1688,
                        -3948, -1130, 2483, -2147, -1634, 2299, -1700, 3549, -3303, -2286,
                    ],
                },
            ],
        };

        for i in 0..4 {
            assert_eq!(t_0[i], expected_t0[i])
        }
    }
}
