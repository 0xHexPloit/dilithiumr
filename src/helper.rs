use rand::{thread_rng, Rng};
use sha3::{Shake256, digest::{Update, ExtendableOutput, XofReader}, Shake256ReaderCore, Shake128, Shake128ReaderCore};
use sha3::digest::core_api::XofReaderCoreWrapper;
use crate::algebra::{matrix::PolyMatrix};
use crate::algebra::poly::Polynomial;
use crate::algebra::vec::PolyVec;

pub fn random_bytes<const N: usize>() -> [u8; N] {
    let mut rng = thread_rng();
    let mut output = [0u8; N];
    output.iter_mut().for_each(|val| *val = rng.gen::<u8>());
    output
}

pub fn shake_128_reader(input: &[u8]) -> XofReaderCoreWrapper<Shake128ReaderCore> {
    let mut hasher = Shake128::default();
    hasher.update(input);
    hasher.finalize_xof()
}


pub fn shake_256_reader(input: &[u8]) -> XofReaderCoreWrapper<Shake256ReaderCore> {
    let mut hasher = Shake256::default();
    hasher.update(input);
    hasher.finalize_xof()
}

pub fn shake_256<const N: usize>(input: &[u8]) -> [u8; N] {
    let mut output = [0u8; N];
    let mut reader = shake_256_reader(input);
    reader.read(&mut output);
    output
}

pub fn expand_a<const K: usize, const L: usize>(rho: &[u8]) -> PolyMatrix<K, L> {
    let matrix = PolyMatrix::<K, L>::default();


    for i in 0..K {
        for j in 0..L {
            let mut poly = matrix[i][j];
            Polynomial::fill_uniform_random(rho, ((i << 8) + j ) as u16, &mut poly);
            poly.ntt()
        }
    }

    matrix
}

pub fn expand_s<const K: usize, const L: usize, const ETA: usize>(rho_prime: &[u8]) -> (PolyVec<K>, PolyVec<L>) {
    let mut vec_one = PolyVec::<K>::default();
    let mut nonce = 0u16;
    vec_one.iter_mut().for_each(|poly| {
        nonce += 1;
        Polynomial::fill_uniform_random_using_eta::<ETA>(rho_prime, nonce, poly);
    });

    nonce = L as u16;
    let mut vec_two = PolyVec::<L>::default();
    vec_two.iter_mut().for_each(|poly| {
        nonce += 1;
        Polynomial::fill_uniform_random_using_eta::<ETA>(rho_prime, nonce, poly);
    });

    (vec_one, vec_two)
}