use algorithm::Dilithium;

mod constants;
mod algorithm;
mod helper;
mod algebra;
mod reduce;
mod round;
mod packing;
mod rejection;

pub const DILITHIUM2: Dilithium<
    4,
    4,
    2,
    2528,
    1312,
> = Dilithium;

pub const DILITHIUM3: Dilithium<
    6,
    5,
    4,
    4000,
    1952,
> = Dilithium;
pub const DILITHIUM5: Dilithium<
    8,
    7,
    2,
    4864,
    2592,
> = Dilithium;

