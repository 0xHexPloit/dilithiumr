use rand::{thread_rng, Rng};
use sha3::digest::core_api::XofReaderCoreWrapper;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128, Shake128ReaderCore, Shake256, Shake256ReaderCore,
};

pub fn random_bytes<const N: usize>() -> [u8; N] {
    let mut rng = thread_rng();
    let mut output = [0u8; N];
    output.iter_mut().for_each(|val| *val = rng.gen::<u8>());
    output
}

pub fn shake_128_reader(data: &[&[u8]]) -> XofReaderCoreWrapper<Shake128ReaderCore> {
    let mut hasher = Shake128::default();
    for item in data {
        hasher.update(*item)
    }
    hasher.finalize_xof()
}

pub fn shake_256_reader(data: &[&[u8]]) -> XofReaderCoreWrapper<Shake256ReaderCore> {
    let mut hasher = Shake256::default();

    for item in data {
        hasher.update(*item)
    }

    hasher.finalize_xof()
}

pub fn shake_256<const N: usize>(data: &[&[u8]]) -> [u8; N] {
    let mut output = [0u8; N];
    let mut reader = shake_256_reader(data);
    reader.read(&mut output);
    output
}
