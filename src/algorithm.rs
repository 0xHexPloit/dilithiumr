use crate::helper::{expand_a, expand_s, random_bytes, shake_256};

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

    pub fn keygen(&self) -> ([u8; SIGNATURE_KEY_SIZE], [u8; PUBLIC_KEY_SIZE]) {
        let zeta = random_bytes::<32>();
        let seed_bytes = shake_256::<128>(&zeta);

        let (rho, seed_bytes) = seed_bytes.split_at(32);
        let (rho_prime, k) = seed_bytes.split_at(64);

        let a_matrix = expand_a::<K, L>(rho);

        let (s_one, s_two) = expand_s::<L, K, ETA>(rho_prime);

        let signature_key = [0u8; SIGNATURE_KEY_SIZE];
        let public_key = [0u8; PUBLIC_KEY_SIZE];

        (signature_key, public_key)
    }
}