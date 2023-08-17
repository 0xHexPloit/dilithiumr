use crate::constants::Q;
use algorithm::Dilithium;

mod algebra;
mod algorithm;
mod constants;
mod conversion;
mod helper;
mod packing;
mod reduce;
mod rejection;
mod rounding;

pub const DILITHIUM2: Dilithium<
    4,
    4,
    2,
    { 1 << 17 },
    { (Q as usize - 1) / 88 },
    39,
    78,
    80,
    576,
    768,
    2528,
    1312,
    2420,
> = Dilithium;

pub const DILITHIUM3: Dilithium<
    6,
    5,
    4,
    { 1 << 19 },
    { (Q as usize - 1) / 32 },
    49,
    196,
    55,
    640,
    768,
    4000,
    1952,
    3293,
> = Dilithium;
pub const DILITHIUM5: Dilithium<
    8,
    7,
    2,
    { 1 << 19 },
    { (Q as usize - 1) / 32 },
    60,
    120,
    75,
    640,
    1024,
    4864,
    2592,
    4595,
> = Dilithium;
