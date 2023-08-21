use crate::constants::{Q, ZETA_BYTES};
use algorithm::Dilithium;
#[cfg(feature = "js")]
use wasm_bindgen::prelude::*;

#[cfg(feature = "python")]
use pyo3::prelude::*;

mod algebra;
mod algorithm;
mod constants;
mod conversion;
mod helper;
mod packing;
mod reduce;
mod rejection;
mod rounding;

const DILITHIUM2: Dilithium<
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

const DILITHIUM3: Dilithium<
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
const DILITHIUM5: Dilithium<
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

fn check_seed(seed: Vec<u8>) -> [u8; ZETA_BYTES] {
    if seed.len() != ZETA_BYTES {
        panic!("Invalid seed length");
    }

    let mut custom_seed = [0u8; ZETA_BYTES];
    custom_seed.copy_from_slice(&seed);
    custom_seed
}

#[cfg(feature = "js")]
#[wasm_bindgen]
pub fn dilithium_mode2_keygen(seed: Option<Vec<u8>>) -> Vec<u8> {
    let mut output = [0u8; 3840];

    if let Some(seed) = seed {
        let custom_seed = check_seed(seed);
        let (sk, pk) = DILITHIUM2.keygen(Some(&custom_seed));
        output[..2528].copy_from_slice(&sk);
        output[2528..].copy_from_slice(&pk);
        return output.to_vec()
    } else {
        let (sk, pk) = DILITHIUM2.keygen(None);
        output[..2528].copy_from_slice(&sk);
        output[2528..].copy_from_slice(&pk);
        return output.to_vec()
    }
}

#[cfg(feature = "python")]
#[pyfunction]
pub fn dilithium_mode2_keygen(seed: Option<Vec<u8>>) -> (Vec<u8>, Vec<u8>) {
    if let Some(seed) = seed {
        let custom_seed = check_seed(seed);
        let (sk, pk) = DILITHIUM2.keygen(Some(&custom_seed));
        return (sk.to_vec(), pk.to_vec())
    } else {
        let (sk, pk) = DILITHIUM2.keygen(None);
        return (sk.to_vec(), pk.to_vec())
    }
}

#[cfg(feature = "js")]
#[wasm_bindgen]
pub fn dilithium_mode3_keygen(seed: Option<Vec<u8>>) -> Vec<u8> {
    let mut output = [0u8; 5952];

    if let Some(seed) = seed {
        if seed.len() != ZETA_BYTES {
            panic!("Invalid seed length");
        }
        let custom_seed = check_seed(seed);
        let (sk, pk) = DILITHIUM3.keygen(Some(&custom_seed));
        output[..4000].copy_from_slice(&sk);
        output[4000..].copy_from_slice(&pk);
        return output.to_vec()

    } else {
        let (sk, pk) = DILITHIUM3.keygen(None);
        output[..4000].copy_from_slice(&sk);
        output[4000..].copy_from_slice(&pk);
        return output.to_vec()
    }
}

#[cfg(feature = "python")]
#[pyfunction]
pub fn dilithium_mode3_keygen(seed: Option<Vec<u8>>) -> (Vec<u8>, Vec<u8>) {
    if let Some(seed) = seed {
        if seed.len() != ZETA_BYTES {
            panic!("Invalid seed length");
        }
        let custom_seed = check_seed(seed);
        let (sk, pk) = DILITHIUM3.keygen(Some(&custom_seed));
        return (sk.to_vec(), pk.to_vec())
    } else {
        let (sk, pk) = DILITHIUM3.keygen(None);
        return (sk.to_vec(), pk.to_vec())
    }
}


#[cfg(feature = "js")]
#[wasm_bindgen]
pub fn dilithium_mode5_keygen(seed: Option<Vec<u8>>) -> Vec<u8> {
    let mut output = [0u8; 7456];

    if let Some(seed) = seed {
        if seed.len() != ZETA_BYTES {
            panic!("Invalid seed length");
        }
        let custom_seed = check_seed(seed);
        let (sk, pk) = DILITHIUM5.keygen(Some(&custom_seed));
        output[..4864].copy_from_slice(&sk);
        output[4864..].copy_from_slice(&pk);
        return output.to_vec()
    } else {
        let (sk, pk) = DILITHIUM5.keygen(None);
        output[..4864].copy_from_slice(&sk);
        output[4864..].copy_from_slice(&pk);
        return output.to_vec()
    }
}


#[cfg(feature = "python")]
#[pyfunction]
pub fn dilithium_mode5_keygen(seed: Option<Vec<u8>>) -> (Vec<u8>, Vec<u8>) {
    if let Some(seed) = seed {
        if seed.len() != ZETA_BYTES {
            panic!("Invalid seed length");
        }
        let custom_seed = check_seed(seed);
        let (sk, pk) = DILITHIUM5.keygen(Some(&custom_seed));
        return (sk.to_vec(), pk.to_vec())
    } else {
        let (sk, pk) = DILITHIUM5.keygen(None);
        return (sk.to_vec(), pk.to_vec())
    }
}

#[cfg(feature = "js")]
#[wasm_bindgen]
pub fn dilithium_mode2_sign(key: Vec<u8>, message: Vec<u8>, rsigning: bool) -> Vec<u8> {
    DILITHIUM2.sign(&key, &message, rsigning).to_vec()
}

#[cfg(feature = "python")]
#[pyfunction]
pub fn dilithium_mode2_sign(key: Vec<u8>, message: Vec<u8>, rsigning: bool) -> Vec<u8> {
    DILITHIUM2.sign(&key, &message, rsigning).to_vec()
}

#[cfg(feature = "js")]
#[wasm_bindgen]
pub fn dilithium_mode3_sign(key: Vec<u8>, message: Vec<u8>, rsigning: bool) -> Vec<u8> {
    DILITHIUM3.sign(&key, &message, rsigning).to_vec()
}

#[cfg(feature = "python")]
#[pyfunction]
pub fn dilithium_mode3_sign(key: Vec<u8>, message: Vec<u8>, rsigning: bool) -> Vec<u8> {
    DILITHIUM3.sign(&key, &message, rsigning).to_vec()
}

#[cfg(feature = "js")]
#[wasm_bindgen]
pub fn dilithium_mode5_sign(key: Vec<u8>, message: Vec<u8>, rsigning: bool) -> Vec<u8> {
    DILITHIUM5.sign(&key, &message, rsigning).to_vec()
}

#[cfg(feature = "python")]
#[pyfunction]
pub fn dilithium_mode5_sign(key: Vec<u8>, message: Vec<u8>, rsigning: bool) -> Vec<u8> {
    DILITHIUM5.sign(&key, &message, rsigning).to_vec()
}

#[cfg(feature = "js")]
#[wasm_bindgen]
pub fn dilithium_mode2_verify(key: Vec<u8>, message: Vec<u8>, signature: Vec<u8>) -> bool {
    DILITHIUM2.verify(&key, &message, &signature)
}

#[cfg(feature = "python")]
#[pyfunction]
pub fn dilithium_mode2_verify(key: Vec<u8>, message: Vec<u8>, signature: Vec<u8>) -> bool {
    DILITHIUM2.verify(&key, &message, &signature)
}

#[cfg(feature = "js")]
#[wasm_bindgen]
pub fn dilithium_mode3_verify(key: Vec<u8>, message: Vec<u8>, signature: Vec<u8>) -> bool {
    DILITHIUM3.verify(&key, &message, &signature)
}

#[cfg(feature = "python")]
#[pyfunction]
pub fn dilithium_mode3_verify(key: Vec<u8>, message: Vec<u8>, signature: Vec<u8>) -> bool {
    DILITHIUM3.verify(&key, &message, &signature)
}

#[cfg(feature = "js")]
#[wasm_bindgen]
pub fn dilithium_mode5_verify(key: Vec<u8>, message: Vec<u8>, signature: Vec<u8>) -> bool {
    DILITHIUM5.verify(&key, &message, &signature)
}

#[cfg(feature = "python")]
#[pyfunction]
pub fn dilithium_mode5_verify(key: Vec<u8>, message: Vec<u8>, signature: Vec<u8>) -> bool {
    DILITHIUM5.verify(&key, &message, &signature)
}

#[cfg(feature = "python")]
#[pymodule]
fn dilithiumr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(dilithium_mode2_keygen, m)?)?;
    m.add_function(wrap_pyfunction!(dilithium_mode3_keygen, m)?)?;
    m.add_function(wrap_pyfunction!(dilithium_mode5_keygen, m)?)?;

    m.add_function(wrap_pyfunction!(dilithium_mode2_sign, m)?)?;
    m.add_function(wrap_pyfunction!(dilithium_mode3_sign, m)?)?;
    m.add_function(wrap_pyfunction!(dilithium_mode5_sign, m)?)?;

    m.add_function(wrap_pyfunction!(dilithium_mode2_verify, m)?)?;
    m.add_function(wrap_pyfunction!(dilithium_mode3_verify, m)?)?;
    m.add_function(wrap_pyfunction!(dilithium_mode5_verify, m)?)?;

    Ok(())
}
