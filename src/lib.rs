use algorithm::Dilithium;

mod constants;
mod algorithm;
mod helper;
mod algebra;


pub const DILITHIUM2: Dilithium<
    4,
    4,
    2420,
    1312,
    false
> = Dilithium;
pub const DILITHIUM2_AES: Dilithium<
    4,
    4,
    2420,
    1312,
    true
> = Dilithium;

pub const DILITHIUM3: Dilithium<
    6,
    5,
    3293,
    1952,
    false
> = Dilithium;
pub const DILITHIUM3_AES: Dilithium<
    6,
    5,
    3293,
    1952,
    true
> = Dilithium;
pub const DILITHIUM5: Dilithium<
    8,
    7,
    4595,
    2592,
    false
> = Dilithium;
pub const DILITHIUM5_AES: Dilithium<
    8,
    7,
    4595,
    2592,
    true
> = Dilithium;

