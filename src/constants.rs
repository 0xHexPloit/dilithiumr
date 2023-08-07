const Q: usize = 8380417; // 2**23 - 2**13 + 1;
const D: u8 = 13; // Dropped bits from t
const N: u16 = 256; // Ring degree

#[derive(Debug)]
pub struct DilithiumParams {
    pub k: usize,
    pub l: usize,

}

pub mod dilithium2 {
    use super::DilithiumParams;

    const K: usize = 4;
    const L: usize = 4;

    pub const PARAMS: DilithiumParams = DilithiumParams {
        k: K,
        l: L
    };
}

pub mod dilithium3 {
    use super::DilithiumParams;

    const K: usize = 6;
    const L: usize = 5;

    pub const PARAMS: DilithiumParams = DilithiumParams {
        k: K,
        l: L
    };
}


pub mod dilithium5 {
    use super::DilithiumParams;

    const K: usize = 8;
    const L: usize = 7;

    pub const PARAMS: DilithiumParams = DilithiumParams {
        k: K,
        l: L
    };
}
