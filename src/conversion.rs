pub fn u16_to_bytes_le(val: u16) -> [u8; 2] {
    let mut output = [0u8; 2];
    output[0] = (val & 0xFF) as u8;
    output[1] = ((val >> 8) & 0xFF) as u8;
    output
}
