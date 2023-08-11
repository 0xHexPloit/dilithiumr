use crate::algebra::poly::Polynomial;
use crate::algebra::vec::PolyVec;
use crate::constants::{BIG_K_BYTES, POLY_T0_PACKED_BYTES, POLY_T1_PACKED_BYTES, RHO_BYTES, TR_BYTES};

pub fn pack_public_key<const N: usize, const K: usize>(rho: &[u8], t_one: &PolyVec<K>) -> [u8; N]  {
    let mut public_key = [0u8; N];
    public_key[..RHO_BYTES].copy_from_slice(rho);

    for i in 0..K {
        let offset = RHO_BYTES + i * POLY_T1_PACKED_BYTES;
        Polynomial::pack_t1(
            &mut public_key[offset..offset+POLY_T1_PACKED_BYTES],
            &t_one[i]
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
    t_zero: &PolyVec<K>
) -> [u8; N] {
    let mut private_key = [0u8; N];
    let mut offset = 0;
    let poly_eta_bytes:usize = if ETA == 2 { 96 } else {128};

    private_key[offset..RHO_BYTES].copy_from_slice(rho);

    offset += RHO_BYTES;
    private_key[offset..offset+BIG_K_BYTES].copy_from_slice(big_k);

    offset += BIG_K_BYTES;
    private_key[offset..offset + TR_BYTES].copy_from_slice(tr);

    offset += TR_BYTES;
    for i in 0..L {
        Polynomial::pack_eta::<ETA>(
            &mut private_key[offset + i * poly_eta_bytes..],
            &s_one[i]
        );
    }
    offset += L * poly_eta_bytes;

    for i in 0..K {
        Polynomial::pack_eta::<ETA>(
            &mut private_key[offset + i * poly_eta_bytes..],
            &s_two[i]
        );
    }
    offset += K * poly_eta_bytes;

    for i in 0..K {
        Polynomial::pack_t0(
            &mut private_key[offset + i * POLY_T0_PACKED_BYTES..],
            &t_zero[i]
        );
    }

    private_key
}