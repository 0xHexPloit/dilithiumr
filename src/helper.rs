use rand::{thread_rng, Rng};
use sha3::{Shake256, digest::{Update, ExtendableOutput, XofReader}};
use crate::algebra::{matrix::PolyMatrix, poly::Polynomial};

pub fn random_bytes<const N: usize>() -> [u8; N] {
    let mut rng = thread_rng();
    let mut output = [0u8; N];
    output.iter_mut().for_each(|val| *val = rng.gen::<u8>());
    output
}

pub fn shake_256<const N: usize>(input: &[u8]) -> [u8; N] {
    let mut hasher = Shake256::default();
    hasher.update(input);
    let mut reader = hasher.finalize_xof();
    let mut output = [0u8; N];
    reader.read(&mut output);
    output
}

pub fn expand_a<const K: usize, const L: usize>(rho: &[u8]) {
    let mut matrix = PolyMatrix::<K, L>::default();

    for i in 0..K {
        for j in 0..L {
            matrix.set(i, j, Polynomial::uniform_random(rho, (i << 8) + j));
        }
    }

}