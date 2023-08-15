use crate::constants::{SEED_BYTES, TR_BYTES, ZETA_BYTES};
use crate::helper::{expand_a, expand_s, random_bytes, shake_256};
use crate::packing::{pack_public_key, pack_signing_key};

pub struct Dilithium<
    const K: usize,
    const L: usize,
    const ETA: usize,
    const SIGNATURE_KEY_SIZE: usize,
    const PUBLIC_KEY_SIZE: usize,
>;


impl <
    const K: usize,
    const L: usize,
    const ETA: usize,
    const SIGNATURE_KEY_SIZE: usize,
    const PUBLIC_KEY_SIZE: usize,
>Dilithium<
    K,
    L,
    ETA,
    SIGNATURE_KEY_SIZE,
    PUBLIC_KEY_SIZE,
> {

    pub fn keygen(&self, custom_zeta: Option<&[u8; ZETA_BYTES]>) -> ([u8; SIGNATURE_KEY_SIZE], [u8; PUBLIC_KEY_SIZE]) {
        let zeta;

        if let Some(custom_zeta) = custom_zeta {
            zeta = *custom_zeta;
        } else {
            zeta = random_bytes::<ZETA_BYTES>();
        }

        let seed_bytes = shake_256::<SEED_BYTES>(&zeta);

        let (rho, seed_bytes) = seed_bytes.split_at(32);
        let (rho_prime, big_k) = seed_bytes.split_at(64);

        let a_matrix = expand_a::<K, L>(rho);

        let (s_one, s_two) = expand_s::<L, K, ETA>(rho_prime);

        let s_one_hat = s_one.to_ntt();
        let mut t = a_matrix * &s_one_hat;
        t += &s_two;
        t.caddq();

        let (t_zero, t_one) = t.power2round();

        let public_key = pack_public_key::<PUBLIC_KEY_SIZE, K>(rho, &t_one);

        let tr = shake_256::<TR_BYTES>(&public_key);
        let signing_key = pack_signing_key::<SIGNATURE_KEY_SIZE, K, L, ETA>(
            rho,
            big_k,
            &tr,
            &s_one,
            &s_two,
            &t_zero
        );

        (signing_key, public_key)
    }
}