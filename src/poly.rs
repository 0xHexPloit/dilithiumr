use sha3::digest::XofReader;
use crate::algebra::poly::Polynomial;
use crate::constants::{N, POLY_UNIFORM_ETA_2_NB_BLOCKS, POLY_UNIFORM_ETA_4_NB_BLOCKS, Q, REJ_UNIFORM_BUFLEN, STREAM128_BLOCKBYTES, ZETAS};
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

    while counter < remaining_coefficients && current_index + 3 <= buffer_len {
        // In the paper, coefficient value is computed like this: b2_prime * 2**16 + b1 * 2**8 + b0
        let mut coefficient_value = buffer[current_index] as i32 & 0xFF;
        current_index += 1;
        coefficient_value |= (buffer[current_index] as i32 & 0xFF) << 8;
        current_index += 1;
        coefficient_value |= (buffer[current_index] as i32 & 0xFF) << 16;
        coefficient_value &= 0x7FFFFF; // Logical and between b2 and 2**128-1

        if coefficient_value < Q as i32 {
            poly[offset + counter] = coefficient_value;
            counter += 1;
        }
    }

    counter
}

fn rejection_eta<const ETA: usize>(
    poly: &mut Polynomial,
    remaining_coefficients: usize,
    offset: usize,
    buffer: &[u8],
    buffer_len: usize
) -> usize {
    0
}


fn poly_fill_uniform_random_base<const ARR_SIZE: usize>(
    rho: &[u8],
    nonce: u16,
    poly:
    &mut Polynomial,
    rejection_fn: impl Fn(&mut Polynomial, usize, usize, &[u8], usize) -> usize
)
{
    let mut xof_input = [0u8; 32 + 2];
    xof_input[..32].copy_from_slice(rho);
    xof_input[32] = (nonce & 0xFF) as u8;
    xof_input[33] = ((nonce >> 8)  & 0xFF) as u8;

    let mut reader = shake_256_reader(&xof_input);
    let mut buffer = [0u8; ARR_SIZE];
    let mut buffer_len = ARR_SIZE;
    reader.read(&mut buffer);


    let mut counter = rejection_fn(
        poly,
        N,
        0,
        &buffer,
        buffer_len
    );

    while counter < N {
        let mut offset = buffer_len % 3;
        for i in 0..offset {
            buffer[i] = buffer[buffer_len - offset + 1];
        }
        reader.read(&mut buffer[offset..STREAM128_BLOCKBYTES]);
        buffer_len = STREAM128_BLOCKBYTES + offset;
        counter += rejection_fn(
            poly,
            N - counter,
            counter,
            &buffer,
            buffer_len
        );
    }
}


pub fn poly_fill_uniform_random_using_eta<const ETA: usize>(rho: &[u8], nonce: u16, poly: &mut Polynomial) {
    if ETA == 2 {
        poly_fill_uniform_random_base::<POLY_UNIFORM_ETA_2_NB_BLOCKS>(
            rho,
            nonce,
            poly,
            rejection_eta::<ETA>
        );
    } else {
        poly_fill_uniform_random_base::<POLY_UNIFORM_ETA_4_NB_BLOCKS>(
            rho,
            nonce,
            poly,
            rejection_eta::<ETA>
        );
    }
}


/// This function fills the coefficients of a polynomial with values sampled from [0, Q-1].
///
/// # Arguments
/// * `rho` -  An array containing 32 random bytes
/// * `nonce` - A random integer
/// * `poly` - A polynomial from which the coefficients will be converted into the NTT domain.
pub fn poly_fill_uniform_random(rho: &[u8], nonce: u16, poly: &mut Polynomial) {
    poly_fill_uniform_random_base::<REJ_UNIFORM_BUFLEN>(
        rho,
        nonce,
        poly,
        rejection_uniform
    );
}


/// This function converts the coefficients of a polynomial into the NTT domain.
///
/// # Arguments
/// * `poly` - A polynomial from which the coefficients will be converted into the NTT domain.
pub fn poly_to_ntt(poly: &mut Polynomial) {
    let mut zeta_index = 0usize;
    let mut layer = 128u8;

    while layer != 0 {
        let mut offset = 0;
        while offset < N {
            zeta_index += 1;
            let zeta = ZETAS[zeta_index] as i64;

            let mut j = offset;
            while j < offset + layer as usize {
                let t = montgomery_reduce(zeta * poly[j + layer as usize] as i64);

                poly.set(j + layer as usize, poly[j] - t);
                poly.set(j, poly[j] + t);

                j += 1;
            }

            offset = j + layer as usize;
        }
        layer >>= 1;
    }
}

/// This function converts the coefficients of a polynomial from the NTT domain to the montgomery
/// domain.
///
/// Be aware that due to the Montgomery reduction optimization, we have:
///  p != poly_from_ntt(poly_to_ntt(p))
///
/// # Arguments
/// * `poly` - A polynomial from which the coefficients will be converted into the Montgomery domain.
pub fn poly_from_ntt(poly: &mut Polynomial) {
    let mut layer = 1;
    let mut zeta_index = 256;
    const F: i64 = 41978; // mont**2/256

    while layer < N {
        let mut offset = 0;

        while offset < N {
            zeta_index -= 1;
            let zeta = -ZETAS[zeta_index] as i64;
            let mut j = offset;

            while j < offset + layer {
                let temp_value = poly[j];
                poly[j] = temp_value + poly[j + layer];
                poly[j + layer] = temp_value - poly[j + layer];
                poly[j + layer] = montgomery_reduce(zeta * poly[j + layer] as i64);
                j  += 1;
            }

            offset = j + layer
        }
        layer <<= 1;
    }

    for j in 0..N {
        poly[j] = montgomery_reduce(F * poly[j] as i64);
    }
}


#[cfg(test)]
mod tests {
    use crate::algebra::poly::Polynomial;
    use crate::constants::N;
    use crate::poly::{poly_from_ntt, poly_to_ntt};
    const POLY_COEFFICIENTS: [i32; N] = [0, 2, 1, 0, 1, 2, 1, 2, 1, 2, 1, 1, 0, 2, 0, 2, 0, 0, 2, 2, 1, 0, 2, 0, 2, 2, 1, 0, 2, 0, 0, 0, 2, 1, 2, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 2, 0, 2, 0, 1, 1, 1, 1, 2, 0, 1, 2, 0, 2, 2, 1, 1, 0, 2, 1, 1, 1, 0, 0, 2, 0, 1, 0, 1, 2, 2, 1, 0, 0, 1, 0, 0, 1, 2, 1, 2, 0, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1, 0, 0, 0, 1, 2, 0, 2, 1, 2, 1, 0, 1, 2, 0, 2, 2, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 2, 0, 2, 1, 2, 2, 1, 0, 2, 1, 1, 1, 2, 2, 1, 0, 1, 1, 0, 2, 0, 2, 1, 1, 1, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 2, 2, 2, 1, 1, 2, 0, 0, 0, 2, 2, 1, 0, 1, 2, 1, 2, 0, 0, 2, 2, 2, 1, 0, 1, 0, 2, 0, 1, 0, 1, 2, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 2, 0, 0, 2, 0, 2, 0, 1, 2, 2, 0, 2, 2, 0, 0, 2, 1, 2, 2, 2, 0, 2, 2, 1, 0, 2, 2, 1, 1, 0, 1, 0, 0, 0, 0, 2, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0, 1, 0, 0, 0];
    const NTT_COEFFICIENTS: [i32; N] = [1324937, -473917, -2439054, -10231118, 7883660, 836686, 4848677, -2272999, 1963327, -1510497, -2897401, -1684829, -4536221, -6097265, -10275862, -15835020, -239926, 2029004, -2959764, 1915074, -12255912, -5176090, -4573092, -406606, -3755953, 2767767, 2519330, 8192708, -6832530, -699898, 5801576, 935400, 7472675, 10990953, 5523952, 1939072, -3057062, 130916, 6683053, 2353425, 4530773, 4871577, -1807019, -1408155, 7951805, 7565533, 7194716, 3843994, 1904348, -3849058, 6460844, 4873286, -327803, 5844995, 5025194, 763154, -1592490, -1683540, -10476880, -5809438, -172396, 4126920, 1208165, 7551435, 7621375, 11556285, 3226148, 7012680, 6982375, -838523, 3318031, -1684739, -3564720, 3917146, 6799143, 4455491, 3096178, -4016104, -2036888, -272646, -1235904, 4231312, -2400736, 737704, 6057633, 10256975, 2220081, 6181087, 7235493, 13856177, 10964618, 16686392, 7405180, 6788062, 12129779, 7356547, -5472207, -2686819, -6160316, -3591874, 592587, -3416893, 2279304, 7354778, -1196158, -3517162, -3513080, -10027440, 3816012, 10256176, 183087, -1974251, -292218, -1612134, 8267103, 1931533, 8230937, 6731531, 513025, 7431967, -7588736, -4109412, -5498606, 1714950, -9678835, -7542319, -159506, -1593504, -116176, -546524, 8162666, 4127334, 3021272, -4749084, -7061437, -5126963, 5185687, 617469, 8925260, 4691928, 7896362, 2647094, 2891498, 8597454, -3944949, 3792265, -7713318, 320670, 1648974, -1649926, 5740220, 5376896, -1587757, -9000271, 2606666, -1780722, -6501489, -7571601, -2435001, 3863983, -8224044, -8302902, -11292931, -5844527, -10065113, -14814547, -6456804, -8118580, -9741109, -1690563, 2747567, -3929995, -4336170, -11593456, -10092051, -10468567, 7894864, 822938, 975578, -3828584, -2120166, 3236766, -2849792, -5517388, -10711699, -2989097, -11108684, -4594932, -2489041, 2203609, 2127561, -2869773, -7327707, -10202521, -15024533, -13025039, -6015338, -12143668, -7404088, -252250, 392371, 175705, 2694655, 4981345, -6787192, 412084, -6735759, -9727809, -5975696, -13542336, -4325145, -2095579, 5933577, 5227575, 936622, -5843130, 9338548, 6500968, -977158, 6608042, 3423279, -3826809, 432508, 5547534, -2160431, -2581961, 8067361, 2263787, 10987804, 3468808, 8007044, 12241772, 9168100, 6903642, 1032720, 5793854, 4911550, 11107972, 7420186, 5541760, -3048105, 458339, -9674009, -2131069, -1338307, 297905, -5995055, -3701955, -82750, -482096, -1374443, 5753413, -1996480, -9061298, -6695919, 1033669];

    const POLY_COEFFICIENTS_MONTGOMERY: [i32; N] = [0, 7167, -4186625, 0, -4186625, 7167, -4186625, 7167, -4186625, 7167, -4186625, -4186625, 0, 7167, 0, 7167, 0, 0, 7167, 7167, -4186625, 0, 7167, 0, 7167, 7167, -4186625, 0, 7167, 0, 0, 0, 7167, -4186625, 7167, -4186625, -4186625, -4186625, 0, -4186625, 0, -4186625, -4186625, -4186625, 0, 7167, 0, 7167, 0, -4186625, -4186625, -4186625, -4186625, 7167, 0, -4186625, 7167, 0, 7167, 7167, -4186625, -4186625, 0, 7167, -4186625, -4186625, -4186625, 0, 0, 7167, 0, -4186625, 0, -4186625, 7167, 7167, -4186625, 0, 0, -4186625, 0, 0, -4186625, 7167, -4186625, 7167, 0, 0, 7167, 0, -4186625, -4186625, 7167, -4186625, -4186625, 7167, -4186625, 0, 0, 0, -4186625, 7167, 0, 7167, -4186625, 7167, -4186625, 0, -4186625, 7167, 0, 7167, 7167, -4186625, -4186625, 7167, -4186625, -4186625, -4186625, -4186625, 7167, -4186625, -4186625, -4186625, 7167, 0, 7167, -4186625, 7167, 7167, -4186625, 0, 7167, -4186625, -4186625, -4186625, 7167, 7167, -4186625, 0, -4186625, -4186625, 0, 7167, 0, 7167, -4186625, -4186625, -4186625, 0, 0, 7167, -4186625, -4186625, -4186625, -4186625, -4186625, -4186625, -4186625, 0, -4186625, -4186625, 7167, 7167, 7167, -4186625, -4186625, 7167, 0, 0, 0, 7167, 7167, -4186625, 0, -4186625, 7167, -4186625, 7167, 0, 0, 7167, 7167, 7167, -4186625, 0, -4186625, 0, 7167, 0, -4186625, 0, -4186625, 7167, 0, -4186625, -4186625, 0, 0, -4186625, -4186625, -4186625, 0, -4186625, 0, 0, 7167, 0, 0, 7167, 0, 7167, 0, -4186625, 7167, 7167, 0, 7167, 7167, 0, 0, 7167, -4186625, 7167, 7167, 7167, 0, 7167, 7167, -4186625, 0, 7167, 7167, -4186625, -4186625, 0, -4186625, 0, 0, 0, 0, 7167, 0, -4186625, 7167, 0, 0, 7167, 0, 7167, 7167, 0, -4186625, 0, 0, 0];

    #[test]
    fn test_should_convert_polynomial_coefficients_into_ntt_domain() {
        let mut poly = Polynomial{coefficients:  POLY_COEFFICIENTS};
        poly_to_ntt(&mut poly);
        assert_eq!(poly.coefficients, NTT_COEFFICIENTS);
    }


    #[test]
    fn test_should_convert_polynomial_coefficients_from_ntt_domain_to_montgomery_domain() {
        let mut poly = Polynomial { coefficients: NTT_COEFFICIENTS };
        poly_from_ntt(&mut poly);
        assert_eq!(poly.coefficients, POLY_COEFFICIENTS_MONTGOMERY)
    }

}