pub struct Dilithium<
    const K: u8,
    const L: u8,
    const SIGNATURE_KEY_SIZE: usize,
    const PUBLIC_KEY_SIZE: usize
>;


impl <
    const K: u8,
    const L: u8,
    const SIGNATURE_KEY_SIZE: usize,
    const PUBLIC_KEY_SIZE: usize
>Dilithium<
    K,
    L,
    SIGNATURE_KEY_SIZE,
    PUBLIC_KEY_SIZE
> {

    pub fn keygen(&self) -> ([u8; SIGNATURE_KEY_SIZE], [u8; PUBLIC_KEY_SIZE]) {



        let signature_key = [0u8; SIGNATURE_KEY_SIZE];
        let public_key = [0u8; PUBLIC_KEY_SIZE];

        (signature_key, public_key)
    }
}