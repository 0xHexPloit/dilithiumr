use algorithm::Dilithium;

mod constants;
mod algorithm;


pub const DILITHIUM2: Dilithium = Dilithium {params: constants::dilithium2::PARAMS};
pub const DILITHIUM3: Dilithium = Dilithium {params: constants::dilithium3::PARAMS};
pub const DILITHIUM5: Dilithium = Dilithium {params: constants::dilithium5::PARAMS};