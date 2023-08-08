use rand::{thread_rng, Rng};
use sha3::{Shake256, digest::{Update, ExtendableOutput, XofReader}, Shake256ReaderCore};
use sha3::digest::core_api::XofReaderCoreWrapper;
use crate::algebra::{matrix::PolyMatrix, poly::Polynomial};
use crate::constants::{N, Q, REJ_UNIFORM_BUFLEN, STREAM128_BLOCKBYTES};

pub fn random_bytes<const N: usize>() -> [u8; N] {
    let mut rng = thread_rng();
    let mut output = [0u8; N];
    output.iter_mut().for_each(|val| *val = rng.gen::<u8>());
    output
}


fn shake_256_reader(input: &[u8]) -> XofReaderCoreWrapper<Shake256ReaderCore> {
    let mut hasher = Shake256::default();
    hasher.update(input);
    let mut reader = hasher.finalize_xof();
    reader
}

pub fn shake_256<const N: usize>(input: &[u8]) -> [u8; N] {
    let mut output = [0u8; N];
    let mut reader = shake_256_reader(input);
    reader.read(&mut output);
    output
}


fn rejection_uniform(
    poly: &mut Polynomial,
    remaining_coefficients: usize,
    offset: usize,
    buffer: &[u8],
    buffer_len: usize
) -> usize{
    let mut counter = 0;
    let mut current_index = 0;
    let mut coefficient_value: usize = 0;

    while counter < remaining_coefficients && current_index + 3 <= buffer_len {
        // In the paper, coefficient value is computed like this: b2_prime * 2**16 + b1 * 2**8 + b0
        coefficient_value = buffer[current_index] as usize & 0xFF;
        current_index += 1;
        coefficient_value |= (buffer[current_index] as usize & 0xFF) << 8;
        current_index += 1;
        coefficient_value |= (buffer[current_index] as usize & 0xFF) << 16;
        coefficient_value &= 0x7FFFFF; // Logical and between b2 and 2**128-1

        if coefficient_value < Q {
            poly.set(offset + counter, coefficient_value);
            counter += 1;
        }
    }

    counter
}


pub fn poly_fill_uniform_random(rho: &[u8], nonce: usize, poly: &mut Polynomial) {
    let mut xof_input = [0u8; 32 + 2];
    xof_input[..32].copy_from_slice(rho);
    xof_input[32] = (nonce & 0xFF) as u8;
    xof_input[33] = ((nonce >> 8)  & 0xFF) as u8;

    let mut reader = shake_256_reader(&xof_input);
    let mut buffer = [0u8; REJ_UNIFORM_BUFLEN];
    let mut buffer_len = REJ_UNIFORM_BUFLEN;
    reader.read(&mut buffer);

    let mut counter = 0;
    let mut offset = 0;

    counter = rejection_uniform(
        poly,
        N,
        0,
        &buffer,
        buffer_len
    );

    while counter < N {
        offset = buffer_len % 3;
        for i in 0..offset {
            buffer[i] = buffer[buffer_len - offset + 1];
        }
        reader.read(&mut buffer[offset..STREAM128_BLOCKBYTES]);
        buffer_len = STREAM128_BLOCKBYTES + offset;
        counter += rejection_uniform(
            poly,
            N - counter,
            counter,
            &buffer,
            buffer_len
        );
    }
}



pub fn expand_a<const K: usize, const L: usize>(rho: &[u8]) -> PolyMatrix<K, L> {
    let mut matrix = PolyMatrix::<K, L>::default();

    for i in 0..K {
        for j in 0..L {
            let mut poly = matrix.get_mut(i, j);
            poly_fill_uniform_random(rho, (i << 8) + j, poly);

            // TODO Perform NTT Transformation to the poly
        }
    }

    matrix
}