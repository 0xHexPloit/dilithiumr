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

#[cfg(test)]
mod tests {
    use crate::constants::SIGNING_KEY_TEST;
    use crate::DILITHIUM2;

    #[test]
    fn test_should_produce_the_correct_signing_and_public_keys() {
        let seed = [0u8; 32];
        let (sk, pk) = DILITHIUM2.keygen(Some(&seed));

        let expected_public_key = [
            245, 151, 124, 130, 131, 84, 106, 99, 114, 59, 195, 29, 38, 25, 18, 79, 17, 219, 70,
            88, 100, 51, 54, 116, 29, 248, 23, 87, 213, 173, 48, 98, 142, 244, 52, 176, 225, 204,
            234, 77, 254, 222, 183, 162, 40, 111, 3, 204, 221, 134, 109, 116, 153, 14, 137, 134,
            101, 201, 48, 24, 85, 231, 139, 145, 171, 234, 10, 76, 152, 45, 200, 206, 57, 13, 203,
            77, 249, 97, 45, 148, 20, 25, 69, 137, 251, 169, 36, 88, 177, 21, 162, 220, 75, 0, 142,
            151, 170, 157, 66, 244, 237, 234, 238, 52, 20, 31, 211, 171, 152, 89, 163, 254, 206,
            234, 73, 128, 95, 186, 15, 75, 184, 38, 200, 45, 205, 136, 42, 187, 184, 207, 180, 220,
            181, 241, 160, 14, 214, 152, 197, 59, 98, 55, 15, 58, 13, 104, 46, 11, 87, 59, 98, 104,
            211, 251, 139, 239, 219, 132, 205, 251, 199, 194, 77, 130, 238, 80, 150, 232, 73, 120,
            44, 154, 132, 107, 6, 88, 139, 144, 235, 177, 189, 98, 161, 108, 184, 24, 69, 2, 81,
            125, 40, 119, 54, 123, 106, 36, 81, 40, 174, 23, 96, 143, 50, 4, 111, 166, 35, 5, 96,
            185, 135, 86, 169, 199, 243, 235, 117, 185, 81, 240, 254, 41, 175, 19, 240, 92, 88,
            143, 247, 14, 214, 229, 0, 173, 172, 4, 87, 185, 85, 119, 168, 37, 84, 81, 1, 81, 144,
            89, 223, 64, 24, 184, 233, 78, 22, 33, 37, 130, 19, 127, 161, 101, 94, 38, 159, 250,
            186, 177, 204, 69, 146, 25, 76, 200, 1, 250, 40, 253, 171, 86, 163, 58, 204, 249, 207,
            236, 195, 73, 62, 105, 161, 177, 232, 92, 253, 44, 169, 90, 158, 105, 53, 106, 154, 57,
            249, 147, 61, 8, 249, 228, 66, 70, 11, 184, 160, 131, 152, 186, 73, 205, 243, 160, 189,
            37, 124, 75, 160, 219, 4, 108, 212, 117, 234, 27, 82, 16, 181, 107, 35, 13, 20, 27,
            133, 255, 51, 166, 145, 68, 225, 151, 60, 132, 66, 145, 40, 119, 70, 35, 117, 144, 219,
            27, 246, 95, 143, 224, 102, 40, 55, 0, 108, 140, 199, 147, 116, 176, 141, 53, 91, 35,
            84, 86, 133, 25, 41, 153, 154, 122, 95, 110, 186, 126, 92, 222, 216, 245, 210, 167,
            112, 98, 89, 88, 126, 11, 219, 209, 228, 86, 193, 201, 5, 151, 74, 169, 76, 32, 252,
            187, 70, 125, 107, 168, 243, 137, 160, 115, 75, 85, 31, 139, 235, 199, 239, 146, 222,
            189, 101, 80, 213, 80, 216, 100, 150, 144, 139, 102, 194, 240, 220, 221, 89, 151, 99,
            93, 126, 124, 120, 191, 206, 62, 105, 92, 131, 176, 241, 79, 61, 95, 148, 205, 42, 231,
            196, 163, 216, 197, 144, 9, 211, 65, 82, 182, 216, 5, 61, 62, 8, 50, 227, 0, 88, 199,
            63, 219, 107, 167, 165, 207, 73, 190, 147, 131, 59, 29, 107, 7, 201, 118, 132, 63, 158,
            235, 115, 104, 1, 187, 164, 235, 80, 126, 46, 194, 94, 3, 154, 37, 166, 14, 230, 199,
            73, 239, 91, 55, 95, 114, 225, 88, 218, 89, 73, 160, 149, 79, 143, 239, 67, 114, 181,
            217, 6, 174, 12, 169, 97, 90, 209, 197, 91, 51, 243, 112, 202, 81, 221, 27, 5, 148,
            107, 189, 63, 53, 19, 7, 72, 66, 12, 166, 56, 224, 106, 60, 150, 144, 214, 198, 35, 10,
            216, 107, 245, 206, 167, 137, 246, 49, 114, 151, 51, 107, 6, 185, 196, 50, 202, 118, 2,
            18, 235, 188, 66, 206, 251, 147, 48, 153, 11, 169, 30, 118, 121, 106, 49, 245, 222,
            115, 234, 140, 87, 188, 31, 232, 146, 220, 140, 137, 222, 8, 202, 156, 12, 225, 206,
            128, 168, 71, 164, 130, 235, 0, 56, 104, 97, 34, 101, 124, 120, 63, 53, 103, 156, 201,
            8, 51, 231, 44, 34, 29, 19, 155, 124, 147, 83, 78, 32, 33, 71, 182, 124, 11, 246, 213,
            82, 37, 157, 99, 242, 248, 156, 219, 182, 208, 114, 51, 236, 104, 120, 58, 132, 90, 28,
            173, 79, 182, 166, 28, 133, 47, 41, 119, 242, 176, 110, 203, 111, 141, 73, 184, 103,
            193, 214, 196, 169, 23, 139, 196, 246, 102, 121, 22, 210, 47, 254, 187, 110, 42, 211,
            208, 88, 102, 229, 165, 182, 120, 87, 135, 229, 246, 138, 222, 33, 160, 72, 85, 142,
            38, 162, 109, 131, 133, 37, 93, 115, 187, 171, 83, 232, 63, 59, 156, 5, 111, 227, 125,
            166, 203, 133, 109, 222, 109, 192, 133, 43, 96, 188, 53, 100, 207, 92, 230, 135, 19,
            68, 203, 96, 202, 169, 200, 118, 2, 251, 51, 37, 184, 96, 141, 111, 127, 167, 249, 165,
            173, 48, 61, 86, 231, 0, 36, 210, 55, 86, 64, 182, 75, 124, 158, 139, 248, 119, 75, 0,
            102, 234, 231, 204, 253, 89, 3, 17, 157, 37, 134, 244, 0, 57, 66, 89, 223, 36, 21, 83,
            87, 37, 45, 249, 150, 78, 223, 77, 223, 49, 227, 65, 0, 174, 191, 233, 45, 153, 240,
            220, 174, 18, 236, 98, 145, 62, 116, 233, 206, 158, 52, 98, 112, 181, 22, 235, 149,
            213, 42, 236, 159, 171, 113, 39, 11, 71, 90, 129, 94, 65, 27, 218, 146, 11, 97, 181,
            37, 1, 54, 19, 204, 36, 157, 156, 41, 60, 87, 245, 67, 145, 207, 209, 188, 43, 79, 154,
            33, 175, 75, 103, 91, 169, 118, 166, 84, 219, 93, 82, 141, 253, 44, 174, 96, 76, 246,
            141, 222, 201, 246, 106, 22, 147, 145, 189, 64, 4, 156, 248, 44, 221, 96, 57, 159, 110,
            249, 54, 89, 205, 42, 255, 251, 180, 250, 239, 214, 137, 184, 70, 65, 136, 27, 116,
            212, 159, 42, 28, 236, 99, 181, 177, 4, 166, 158, 161, 167, 2, 195, 247, 55, 197, 63,
            186, 241, 146, 171, 247, 190, 173, 40, 94, 98, 11, 112, 174, 20, 101, 11, 121, 180, 74,
            53, 48, 56, 59, 133, 32, 148, 83, 189, 7, 50, 167, 138, 108, 155, 77, 164, 140, 3, 1,
            160, 96, 66, 52, 18, 148, 126, 184, 184, 127, 208, 128, 15, 117, 128, 46, 1, 83, 169,
            203, 45, 63, 177, 182, 152, 224, 143, 43, 44, 11, 247, 234, 59, 37, 104, 143, 4, 40,
            44, 156, 191, 171, 33, 28, 21, 110, 206, 13, 29, 6, 10, 84, 23, 21, 129, 174, 9, 184,
            72, 194, 29, 182, 69, 183, 249, 21, 100, 5, 165, 203, 249, 26, 105, 246, 126, 221, 72,
            151, 175, 53, 250, 122, 66, 49, 94, 199, 23, 101, 111, 77, 195, 150, 59, 171, 61, 183,
            130, 6, 105, 153, 53, 159, 46, 40, 122, 38, 0, 89, 255, 26, 235, 23, 229, 191, 113,
            247, 125, 226, 27, 175, 147, 197, 70, 200, 73, 21, 104, 104, 120, 218, 4, 19, 114, 140,
            181, 146, 162, 33, 88, 211, 45, 156, 118, 165, 14, 12, 18, 218, 165, 46, 214, 47, 226,
            92, 248, 156, 181, 104, 219, 211, 24, 232, 51, 78, 160, 30, 217, 122, 99, 56, 120, 173,
            45, 128, 156, 134, 173, 57, 27, 121, 43, 23, 194, 17, 73, 178, 223, 254, 69, 42, 139,
            238, 38, 26, 136, 152, 214, 226, 68, 187, 172, 4, 24, 21, 97, 42, 41, 35, 29, 71, 166,
            223, 37, 204, 143, 67, 45, 79, 59, 61, 216, 191, 61, 153, 112, 29, 62, 118, 88, 50,
            177, 207, 35, 18, 245, 147, 106, 110, 63, 130, 145, 138, 30, 87, 160, 148, 212, 171,
            188, 202, 67, 37, 26, 225, 207, 211, 126, 128, 86, 182, 150, 138, 213, 246, 43, 126,
            136, 158, 199, 34, 115,
        ];

        assert_eq!(sk, SIGNING_KEY_TEST);
        assert_eq!(pk, expected_public_key);
    }
}
