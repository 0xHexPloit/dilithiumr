use crate::helper::{random_bytes, shake_256};

pub struct Dilithium<
    const K: u8,
    const L: u8,
    const SIGNATURE_KEY_SIZE: usize,
    const PUBLIC_KEY_SIZE: usize,
    const AES_VARIANT: bool
>;


impl <
    const K: u8,
    const L: u8,
    const SIGNATURE_KEY_SIZE: usize,
    const PUBLIC_KEY_SIZE: usize,
    const AES_VARIANT: bool
>Dilithium<
    K,
    L,
    SIGNATURE_KEY_SIZE,
    PUBLIC_KEY_SIZE,
    AES_VARIANT
> {

    pub fn keygen(&self) -> ([u8; SIGNATURE_KEY_SIZE], [u8; PUBLIC_KEY_SIZE]) {
        let zeta = random_bytes::<32>();
        let seed_bytes = shake_256::<128>(&zeta);

        let (rho, seed_bytes) = seed_bytes.split_at(32);
        let (rho_prime, k) = seed_bytes.split_at(64);


        let signature_key = [0u8; SIGNATURE_KEY_SIZE];
        let public_key = [0u8; PUBLIC_KEY_SIZE];

        (signature_key, public_key)
    }
}