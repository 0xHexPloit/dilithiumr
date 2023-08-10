use algorithm::Dilithium;

mod constants;
mod algorithm;
mod helper;
mod algebra;
mod poly;
mod reduce;

pub const DILITHIUM2: Dilithium<
    4,
    4,
    2,
    2420,
    1312,
> = Dilithium;

pub const DILITHIUM3: Dilithium<
    6,
    5,
    4,
    3293,
    1952,
> = Dilithium;
pub const DILITHIUM5: Dilithium<
    8,
    7,
    2,
    4595,
    2592,
> = Dilithium;

