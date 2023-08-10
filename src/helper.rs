use rand::{thread_rng, Rng};
use sha3::{Shake256, digest::{Update, ExtendableOutput, XofReader}, Shake256ReaderCore};
use sha3::digest::core_api::XofReaderCoreWrapper;
use crate::algebra::{matrix::PolyMatrix};
use crate::algebra::poly::Polynomial;
use crate::algebra::vec::PolyVec;
use crate::poly::{poly_fill_uniform_random, poly_fill_uniform_random_using_eta, poly_to_ntt};

pub fn random_bytes<const N: usize>() -> [u8; N] {
    let mut rng = thread_rng();
    let mut output = [0u8; N];
    output.iter_mut().for_each(|val| *val = rng.gen::<u8>());
    output
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
    let mut matrix = PolyMatrix::<K, L>::default();


    for i in 0..K {
        for j in 0..L {
            let poly = matrix.get_mut(i, j);
            poly_fill_uniform_random(rho, (i << 8) + j, poly);
            poly_to_ntt(poly);
        }
    }

    matrix
}

pub fn expand_s<const K: usize, const L: usize, const ETA: usize>(rho_prime: &[u8]) -> (PolyVec<K>, PolyVec<L>) {
    let mut vec_one = PolyVec::<K>::default();
    let mut nonce = 0u16;
    vec_one.into_iter().as_mut_slice().iter_mut().for_each(|poly| {
        nonce += 1;
        poly_fill_uniform_random_using_eta::<ETA>(rho_prime, nonce, poly);
    });

    nonce = L as u16;
    let vec_two = PolyVec::<L>::default();

    (vec_one, vec_two)
}