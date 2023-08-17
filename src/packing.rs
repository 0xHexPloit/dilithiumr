use crate::algebra::poly::Polynomial;
use crate::algebra::vec::PolyVec;
use crate::constants::{
    get_eta_poly_eta_bytes, get_poly_w1_packed_bytes, BIG_K_BYTES, C_TILDE_BYTES, N,
    POLY_T0_PACKED_BYTES, POLY_T1_PACKED_BYTES, RHO_BYTES, TR_BYTES,
};

pub fn pack_public_key<const N: usize, const K: usize>(rho: &[u8], t_one: &PolyVec<K>) -> [u8; N] {
    let mut public_key = [0u8; N];
    public_key[..RHO_BYTES].copy_from_slice(rho);

    for i in 0..K {
        let offset = RHO_BYTES + i * POLY_T1_PACKED_BYTES;
        Polynomial::pack_t1(
            &mut public_key[offset..offset + POLY_T1_PACKED_BYTES],
            &t_one[i],
        );
    }

    public_key
}

pub fn pack_signing_key<const N: usize, const K: usize, const L: usize, const ETA: usize>(
    rho: &[u8],
    big_k: &[u8],
    tr: &[u8],
    s_one: &PolyVec<L>,
    s_two: &PolyVec<K>,
    t_zero: &PolyVec<K>,
) -> [u8; N] {
    let mut private_key = [0u8; N];
    let mut offset = 0;
    let poly_eta_bytes = get_eta_poly_eta_bytes(ETA);

    private_key[offset..RHO_BYTES].copy_from_slice(rho);

    offset += RHO_BYTES;
    private_key[offset..offset + BIG_K_BYTES].copy_from_slice(big_k);

    offset += BIG_K_BYTES;
    private_key[offset..offset + TR_BYTES].copy_from_slice(tr);

    offset += TR_BYTES;
    for i in 0..L {
        Polynomial::pack_eta::<ETA>(&mut private_key[offset + i * poly_eta_bytes..], &s_one[i]);
    }
    offset += L * poly_eta_bytes;

    for i in 0..K {
        Polynomial::pack_eta::<ETA>(&mut private_key[offset + i * poly_eta_bytes..], &s_two[i]);
    }
    offset += K * poly_eta_bytes;

    for i in 0..K {
        Polynomial::pack_t0(
            &mut private_key[offset + i * POLY_T0_PACKED_BYTES..],
            &t_zero[i],
        );
    }

    private_key
}

pub fn unpack_signing_key<'a, const K: usize, const L: usize, const ETA: usize>(
    signing_key: &'a [u8],
    rho: &mut &'a [u8],
    big_k: &mut &'a [u8],
    tr: &mut &'a [u8],
    s_one: &mut PolyVec<L>,
    s_two: &mut PolyVec<K>,
    t_zero: &mut PolyVec<K>,
) {
    let mut offset = 0;
    let poly_eta_bytes = get_eta_poly_eta_bytes(ETA);
    *rho = signing_key[..RHO_BYTES].as_ref();
    offset += RHO_BYTES;
    *big_k = signing_key[offset..offset + BIG_K_BYTES].as_ref();
    offset += BIG_K_BYTES;
    *tr = signing_key[offset..offset + TR_BYTES].as_ref();
    offset += TR_BYTES;

    for i in 0..L {
        Polynomial::unpack_eta::<ETA>(&mut s_one[i], &signing_key[offset + i * poly_eta_bytes..])
    }
    offset += L * poly_eta_bytes;

    for i in 0..K {
        Polynomial::unpack_eta::<ETA>(&mut s_two[i], &signing_key[offset + i * poly_eta_bytes..])
    }
    offset += K * poly_eta_bytes;

    for i in 0..K {
        Polynomial::unpack_t0(&mut t_zero[i], &signing_key[offset + i * poly_eta_bytes..])
    }
}

pub fn pack_w1<const N: usize, const K: usize, const GAMMA_TWO: usize>(
    w_one: &PolyVec<K>,
) -> [u8; N] {
    let mut buffer = [0u8; N];
    let packed_bytes = get_poly_w1_packed_bytes(GAMMA_TWO);

    for i in 0..K {
        Polynomial::pack_w1::<GAMMA_TWO>(&mut buffer[i * packed_bytes..], &w_one[i]);
    }

    buffer
}

fn encode_h_in_signature<const OMEGA: usize, const K: usize>(signature: &mut [u8], h: &PolyVec<K>) {
    let mut k = 0;

    for i in 0..K {
        for j in 0..N {
            if h[i][j] != 0 {
                signature[k] = j as u8;
                k += 1;
            }
        }
        signature[OMEGA + i] = k as u8;
    }
}

pub fn pack_signature<
    const N: usize,
    const K: usize,
    const L: usize,
    const POLY_Z_PACKED_BYTES: usize,
    const OMEGA: usize,
    const GAMMA_ONE: usize,
>(
    c_tile: &[u8],
    z: &PolyVec<L>,
    h: &PolyVec<K>,
) -> [u8; N] {
    let mut offset = 0;
    let mut signature = [0u8; N];
    signature[..C_TILDE_BYTES].copy_from_slice(c_tile);
    offset += C_TILDE_BYTES;

    for i in 0..L {
        Polynomial::pack_z::<GAMMA_ONE>(&mut signature[offset + i * POLY_Z_PACKED_BYTES..], &z[i])
    }
    offset += L * POLY_Z_PACKED_BYTES;

    encode_h_in_signature::<OMEGA, K>(&mut signature[offset..], h);

    signature
}
