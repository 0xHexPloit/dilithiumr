use sha3::digest::XofReader;
use crate::algebra::poly::Polynomial;
use crate::constants::{N, Q, REJ_UNIFORM_BUFLEN, STREAM128_BLOCKBYTES, ZETAS};
use crate::helper::shake_256_reader;
use crate::reduce::montgomery_reduce;

fn rejection_uniform(
    poly: &mut Polynomial,
    remaining_coefficients: usize,
    offset: usize,
    buffer: &[u8],
    buffer_len: usize
) -> usize{
    let mut counter = 0;
    let mut current_index = 0;
    let mut coefficient_value: i32 = 0;

    while counter < remaining_coefficients && current_index + 3 <= buffer_len {
        // In the paper, coefficient value is computed like this: b2_prime * 2**16 + b1 * 2**8 + b0
        coefficient_value = buffer[current_index] as i32 & 0xFF;
        current_index += 1;
        coefficient_value |= (buffer[current_index] as i32 & 0xFF) << 8;
        current_index += 1;
        coefficient_value |= (buffer[current_index] as i32 & 0xFF) << 16;
        coefficient_value &= 0x7FFFFF; // Logical and between b2 and 2**128-1

        if coefficient_value < Q as i32 {
            poly.set(offset + counter, coefficient_value );
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

pub fn poly_to_ntt(poly: &mut Polynomial) {
    let mut zeta_index = 0usize;
    let mut layer = 128u8;

    while layer != 0 {
        let mut start = 0;
        while start < N {
            zeta_index += 1;
            let zeta = ZETAS[zeta_index];

            let mut j = 0;
            while j < start + layer as usize {
                let t = montgomery_reduce(zeta as i64 * poly[j + layer as usize] as i64);

                poly.set(j + layer as usize, poly[j] - t);
                poly.set(j, poly[j] + t);

                j += 1;
            }

            start = j + layer as usize;
        }
        layer >>= 1;
    }
}