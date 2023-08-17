use crate::algebra::matrix::PolyMatrix;
use crate::algebra::poly::Polynomial;
use crate::algebra::vec::PolyVec;
use crate::constants::{
    C_TILDE_BYTES, MU_BYTES, RHO_PRIME_BYTES, SEED_BYTES, TR_BYTES, ZETA_BYTES,
};
use crate::helper::{random_bytes, shake_256};
use crate::packing::{
    pack_public_key, pack_signature, pack_signing_key, pack_w1, unpack_public_key,
    unpack_signature, unpack_signing_key,
};

pub struct Dilithium<
    const K: usize,
    const L: usize,
    const ETA: usize,
    const GAMMA_ONE: usize,
    const GAMMA_TWO: usize,
    const TAU: usize,
    const BETA: usize,
    const OMEGA: usize,
    const POLY_Z_PACKED_BYTES: usize,
    const POLY_W_ONE_PACKED_BYTES: usize,
    const SIGNING_KEY_SIZE: usize,
    const PUBLIC_KEY_SIZE: usize,
    const SIGNATURE_SIZE: usize,
>;

impl<
        const K: usize,
        const L: usize,
        const ETA: usize,
        const GAMMA_ONE: usize,
        const GAMMA_TWO: usize,
        const TAU: usize,
        const BETA: usize,
        const OMEGA: usize,
        const POLY_Z_PACKED_BYTES: usize,
        const POLY_W_ONE_PACKED_BYTES: usize,
        const SIGNING_KEY_SIZE: usize,
        const PUBLIC_KEY_SIZE: usize,
        const SIGNATURE_SIZE: usize,
    >
    Dilithium<
        K,
        L,
        ETA,
        GAMMA_ONE,
        GAMMA_TWO,
        TAU,
        BETA,
        OMEGA,
        POLY_Z_PACKED_BYTES,
        POLY_W_ONE_PACKED_BYTES,
        SIGNING_KEY_SIZE,
        PUBLIC_KEY_SIZE,
        SIGNATURE_SIZE,
    >
{
    pub fn keygen(
        &self,
        custom_zeta: Option<&[u8; ZETA_BYTES]>,
    ) -> ([u8; SIGNING_KEY_SIZE], [u8; PUBLIC_KEY_SIZE]) {
        let zeta;

        if let Some(custom_zeta) = custom_zeta {
            zeta = *custom_zeta;
        } else {
            zeta = random_bytes::<ZETA_BYTES>();
        }

        let seed_bytes = shake_256::<SEED_BYTES>(&[&zeta]);

        let (rho, seed_bytes) = seed_bytes.split_at(32);

        let (rho_prime, big_k) = seed_bytes.split_at(64);

        let a_matrix = PolyMatrix::<K, L>::expand_a(rho);

        let s_one = PolyVec::<L>::expand_s::<ETA>(&rho_prime, 0);
        let s_two = PolyVec::<K>::expand_s::<ETA>(&rho_prime, L as u16);

        let s_one_hat = s_one.to_ntt();
        let mut t = &a_matrix * &s_one_hat;
        t += &s_two;
        t.caddq();

        let (t_zero, t_one) = t.power2round();

        let public_key = pack_public_key::<PUBLIC_KEY_SIZE, K>(rho, &t_one);

        let tr = shake_256::<TR_BYTES>(&[&public_key]);
        let signing_key = pack_signing_key::<SIGNING_KEY_SIZE, K, L, ETA>(
            rho, big_k, &tr, &s_one, &s_two, &t_zero,
        );

        (signing_key, public_key)
    }

    pub fn sign(
        &self,
        signing_key: &[u8],
        message: &[u8],
        randomized_signing: bool,
    ) -> [u8; SIGNATURE_SIZE] {
        // Checking signing key
        if signing_key.len() != SIGNING_KEY_SIZE {
            panic!("Invalid length for signing key!")
        }

        // Extracting data from signing key
        let mut rho: &[u8] = &[];
        let mut big_k: &[u8] = &[];
        let mut tr: &[u8] = &[];

        let mut s_one = PolyVec::<L>::default();
        let mut s_two = PolyVec::<K>::default();
        let mut t_zero = PolyVec::<K>::default();

        unpack_signing_key::<K, L, ETA>(
            signing_key,
            &mut rho,
            &mut big_k,
            &mut tr,
            &mut s_one,
            &mut s_two,
            &mut t_zero,
        );

        let s_one_hat = s_one.to_ntt();
        let s_two_hat = s_two.to_ntt();
        let t_zero_hat = t_zero.to_ntt();

        // Starting the signature process
        let mut kappa = 0;
        let a_matrix = PolyMatrix::<K, L>::expand_a(&rho);
        let mu = shake_256::<MU_BYTES>(&[&tr, message]);

        let rho_prime;
        if randomized_signing {
            rho_prime = random_bytes::<RHO_PRIME_BYTES>();
        } else {
            rho_prime = shake_256::<RHO_PRIME_BYTES>(&[&big_k, &mu]);
        }

        loop {
            // Sample intermediate vector y
            let y = PolyVec::<L>::expand_mask::<POLY_Z_PACKED_BYTES, GAMMA_ONE>(&rho_prime, kappa);
            kappa += 1;

            let y_hat = y.to_ntt();
            let mut w = &a_matrix * &y_hat;
            w.caddq();

            // Decompose w and call the random oracle
            let (w_one, w_zero) = PolyVec::<K>::decompose::<GAMMA_TWO>(&w);
            let w_one_bytes = pack_w1::<POLY_W_ONE_PACKED_BYTES, K, GAMMA_TWO>(&w_one);

            let c_tilde = shake_256::<C_TILDE_BYTES>(&[&mu, &w_one_bytes]);
            let mut c = Polynomial::sample_in_ball::<TAU>(&c_tilde);
            c.ntt();

            // Computing z vector
            let mut z = s_one_hat * &c;
            z += &y;

            if z.check_infinite_norm(GAMMA_ONE - BETA) {
                continue;
            }

            // Computing r0 vector to verify that +subtracting c * s2 does not change high bits of
            // w and low bits don't reveal secret information.
            let mut r_zero = s_two_hat * &c;
            r_zero = w_zero - &r_zero;

            if r_zero.check_infinite_norm(GAMMA_TWO - BETA) {
                continue;
            }

            // Computing hints for w1
            let mut vec = t_zero_hat * &c;
            if vec.check_infinite_norm(GAMMA_TWO) {
                continue;
            }

            vec += &r_zero;

            let mut h = PolyVec::<K>::default();

            let number_ones = PolyVec::<K>::make_hint::<GAMMA_TWO>(&mut h, &w_one, &vec);

            if number_ones > OMEGA {
                continue;
            }

            let signature =
                pack_signature::<SIGNATURE_SIZE, K, L, POLY_Z_PACKED_BYTES, OMEGA, GAMMA_ONE>(
                    &c_tilde, &z, &h,
                );
            return signature;
        }
    }

    pub fn verify(&self, public_key: &[u8], message: &[u8], signature: &[u8]) -> bool {
        // Checking public key
        if public_key.len() != PUBLIC_KEY_SIZE {
            panic!("Invalid length for public key !");
        }

        // Checking signature
        if signature.len() != SIGNATURE_SIZE {
            panic!("Invalid length for the signature !");
        }

        // Unpacking public key
        let mut rho: &[u8] = &[];
        let mut t_one: PolyVec<K> = PolyVec::<K>::default();

        unpack_public_key(&public_key, &mut rho, &mut t_one);

        let mut c_tilde: &[u8] = &[];
        let mut z = PolyVec::<L>::default();
        let mut h = PolyVec::<K>::default();

        let malformed_signature = unpack_signature::<K, L, GAMMA_ONE, POLY_Z_PACKED_BYTES, OMEGA>(
            &signature,
            &mut c_tilde,
            &mut z,
            &mut h,
        );

        if malformed_signature {
            return false;
        }

        // Starting the verification process
        let a_matrix = PolyMatrix::<K, L>::expand_a(&rho);
        let digest = shake_256::<32>(&[&public_key]);
        let mu = shake_256::<MU_BYTES>(&[&digest, &message]);

        true
    }
}
