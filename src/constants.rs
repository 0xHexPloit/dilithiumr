pub const Q: usize = 8380417; // 2**23 - 2**13 + 1;
pub const D: u8 = 13; // Dropped bits from t
pub const N: usize = 256; // Ring degree

const SHAKE128_RATE: usize = 168;
pub const STREAM128_BLOCKBYTES: usize = SHAKE128_RATE;

const REJ_UNIFORM_NBLOCKS: usize  = (768 + STREAM128_BLOCKBYTES - 1) / STREAM128_BLOCKBYTES;
pub const REJ_UNIFORM_BUFLEN: usize = REJ_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES;



