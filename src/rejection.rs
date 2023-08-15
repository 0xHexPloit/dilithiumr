use crate::algebra::poly::Polynomial;
use crate::constants::Q;

pub fn rejection_uniform(
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
        current_index += 1;
        coefficient_value &= 0x7FFFFF; // Logical and between b2 and 2**128-1

        if coefficient_value < Q as i32 {
            poly[offset + counter] = coefficient_value;
            counter += 1;
        }
    }

    counter
}

pub fn rejection_eta<const ETA: usize>(
    poly: &mut Polynomial,
    remaining_coefficients: usize,
    offset: usize,
    buffer: &[u8],
    buffer_len: usize
) -> usize {
    let mut counter = 0;
    let mut current_index = 0;

    while counter < remaining_coefficients && current_index < buffer_len {
        let byte_value = buffer[current_index];
        let mut t_zero = (byte_value & 0x0F) as u32;
        let mut t_one = (byte_value >> 4) as u32;
        current_index += 1;

        if ETA == 2 {
            if t_zero < 15 {
                t_zero = t_zero - ((205 * t_zero) >> 10) * 5;
                poly[offset + counter] = 2 - t_zero as i32;
                counter += 1;
            }

            if t_one < 15 && counter < remaining_coefficients {
                t_one = t_one - ((205 * t_one) >> 10) * 5;
                poly[offset + counter] = 2 - t_one as i32;
                counter += 1;
            }
        } else {
            if t_zero < 9 {
                poly[offset + counter] = 4 - t_zero as i32;
                counter += 1;
            }

            if t_one < 9 && counter < remaining_coefficients {
                poly[offset + counter] = 4 - t_one as i32;
                counter += 1;
            }
        }
    }

    counter
}
