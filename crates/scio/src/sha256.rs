//! SHA-256 (FIPS 180-4), standard library only: the digest the reference scripts record with
//! `digest::digest(algo = "sha256")`, e.g. of R's final `.Random.seed` (`rng_end_sha256` in
//! `summary.json`) and of large outputs, and the checksums of the output manifest (M13).
//!
//! Where the CPU has the SHA extensions (Intel since Ice Lake server, so the Sapphire Rapids nodes;
//! not Cascade Lake), the compression function runs on them ([`shani`]): the same digest, several
//! times faster, which matters for the manifest's checksums of multi-gigabyte outputs. The
//! choice is made at run time and changes no digest.

use std::io::{self, Read};
use std::path::Path;

const K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

/// The SHA-256 of `data`, as 64 lowercase hex digits.
pub fn hex(data: &[u8]) -> String {
    let mut h: [u32; 8] = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19];
    let bits = (data.len() as u64).wrapping_mul(8);
    // the whole 64-byte blocks of the message, then the padded tail (one or two blocks)
    let whole = data.len() / 64 * 64;
    compress_blocks(&mut h, &data[..whole]);
    let mut tail = data[whole..].to_vec();
    tail.push(0x80);
    while tail.len() % 64 != 56 {
        tail.push(0);
    }
    tail.extend_from_slice(&bits.to_be_bytes());
    for chunk in tail.chunks_exact(64) {
        compress(&mut h, chunk);
    }
    h.iter().map(|x| format!("{x:08x}")).collect()
}

/// SHA-256 over data given in pieces (large outputs digested as they are produced, without
/// holding them); [`Sha256::hex`] gives what [`hex`] gives on the pieces concatenated.
#[derive(Clone)]
pub struct Sha256 {
    h: [u32; 8],
    /// bytes not yet compressed (fewer than 64)
    pending: Vec<u8>,
    len: u64,
}

impl Default for Sha256 {
    fn default() -> Self {
        Sha256::new()
    }
}

impl Sha256 {
    pub fn new() -> Sha256 {
        Sha256 { h: [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19], pending: Vec::with_capacity(64), len: 0 }
    }

    pub fn update(&mut self, mut data: &[u8]) {
        self.len = self.len.wrapping_add(data.len() as u64);
        if !self.pending.is_empty() {
            let take = (64 - self.pending.len()).min(data.len());
            self.pending.extend_from_slice(&data[..take]);
            data = &data[take..];
            if self.pending.len() < 64 {
                return;
            }
            let block = std::mem::take(&mut self.pending);
            compress(&mut self.h, &block);
            self.pending = block;
            self.pending.clear();
        }
        let whole = data.len() / 64 * 64;
        compress_blocks(&mut self.h, &data[..whole]);
        self.pending.extend_from_slice(&data[whole..]);
    }

    /// The digest, as 64 lowercase hex digits.
    pub fn hex(mut self) -> String {
        let bits = self.len.wrapping_mul(8);
        let mut tail = std::mem::take(&mut self.pending);
        tail.push(0x80);
        while tail.len() % 64 != 56 {
            tail.push(0);
        }
        tail.extend_from_slice(&bits.to_be_bytes());
        for chunk in tail.chunks_exact(64) {
            compress(&mut self.h, chunk);
        }
        self.h.iter().map(|x| format!("{x:08x}")).collect()
    }
}

/// The bytes and SHA-256 of the file at `path`, read in pieces of 4 MiB.
pub fn file(path: &Path) -> io::Result<(u64, String)> {
    let mut f = std::fs::File::open(path).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 4 << 20];
    let mut n = 0u64;
    loop {
        let k = f.read(&mut buf).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
        if k == 0 {
            break;
        }
        n += k as u64;
        h.update(&buf[..k]);
    }
    Ok((n, h.hex()))
}

/// Whole 64-byte blocks, on the SHA extensions where the CPU has them.
fn compress_blocks(h: &mut [u32; 8], blocks: &[u8]) {
    debug_assert_eq!(blocks.len() % 64, 0);
    if blocks.is_empty() {
        return;
    }
    #[cfg(target_arch = "x86_64")]
    {
        if shani::available() {
            // SAFETY: the CPU has the instructions the function is compiled for (checked above).
            unsafe { shani::compress_blocks(h, blocks) };
            return;
        }
    }
    for chunk in blocks.chunks_exact(64) {
        compress(h, chunk);
    }
}

/// The compression function on the x86 SHA extensions (`sha256rnds2`, `sha256msg1`,
/// `sha256msg2`), in the sequence of Intel's reference code: four rounds per instruction pair,
/// the message schedule four words at a time.
#[cfg(target_arch = "x86_64")]
pub(crate) mod shani {
    use super::K;
    use std::arch::x86_64::*;

    pub fn available() -> bool {
        std::arch::is_x86_feature_detected!("sha") && std::arch::is_x86_feature_detected!("sse4.1") && std::arch::is_x86_feature_detected!("ssse3")
    }

    #[inline(always)]
    unsafe fn k4(i: usize) -> __m128i {
        _mm_set_epi32(K[4 * i + 3] as i32, K[4 * i + 2] as i32, K[4 * i + 1] as i32, K[4 * i] as i32)
    }

    /// The next four words of the message schedule from the previous sixteen.
    #[inline(always)]
    unsafe fn schedule(w0: __m128i, w1: __m128i, w2: __m128i, w3: __m128i) -> __m128i {
        let t1 = _mm_sha256msg1_epu32(w0, w1);
        let t2 = _mm_alignr_epi8(w3, w2, 4);
        _mm_sha256msg2_epu32(_mm_add_epi32(t1, t2), w3)
    }

    /// # Safety
    /// The CPU must have the SHA, SSE4.1 and SSSE3 extensions ([`available`]); `blocks` holds
    /// whole 64-byte blocks.
    #[target_feature(enable = "sha,sse2,ssse3,sse4.1")]
    pub unsafe fn compress_blocks(h: &mut [u32; 8], blocks: &[u8]) {
        // big-endian words: reverse the bytes of each 32-bit lane
        let mask = _mm_set_epi64x(0x0c0d_0e0f_0809_0a0bu64 as i64, 0x0405_0607_0001_0203u64 as i64);
        let p = h.as_ptr() as *const __m128i;
        let dcba = _mm_loadu_si128(p);
        let efgh = _mm_loadu_si128(p.add(1));
        let cdab = _mm_shuffle_epi32(dcba, 0xb1);
        let efgh = _mm_shuffle_epi32(efgh, 0x1b);
        let mut abef = _mm_alignr_epi8(cdab, efgh, 8);
        let mut cdgh = _mm_blend_epi16(efgh, cdab, 0xf0);
        macro_rules! rounds4 {
            ($w:expr, $i:expr) => {{
                let t1 = _mm_add_epi32($w, k4($i));
                cdgh = _mm_sha256rnds2_epu32(cdgh, abef, t1);
                let t2 = _mm_shuffle_epi32(t1, 0x0e);
                abef = _mm_sha256rnds2_epu32(abef, cdgh, t2);
            }};
        }
        for block in blocks.chunks_exact(64) {
            let (abef0, cdgh0) = (abef, cdgh);
            let d = block.as_ptr() as *const __m128i;
            let mut w0 = _mm_shuffle_epi8(_mm_loadu_si128(d), mask);
            let mut w1 = _mm_shuffle_epi8(_mm_loadu_si128(d.add(1)), mask);
            let mut w2 = _mm_shuffle_epi8(_mm_loadu_si128(d.add(2)), mask);
            let mut w3 = _mm_shuffle_epi8(_mm_loadu_si128(d.add(3)), mask);
            rounds4!(w0, 0);
            rounds4!(w1, 1);
            rounds4!(w2, 2);
            rounds4!(w3, 3);
            let mut i = 4;
            while i < 16 {
                let w4 = schedule(w0, w1, w2, w3);
                rounds4!(w4, i);
                (w0, w1, w2, w3) = (w1, w2, w3, w4);
                i += 1;
            }
            abef = _mm_add_epi32(abef, abef0);
            cdgh = _mm_add_epi32(cdgh, cdgh0);
        }
        let feba = _mm_shuffle_epi32(abef, 0x1b);
        let dchg = _mm_shuffle_epi32(cdgh, 0xb1);
        let dcba = _mm_blend_epi16(feba, dchg, 0xf0);
        let hgef = _mm_alignr_epi8(dchg, feba, 8);
        let q = h.as_mut_ptr() as *mut __m128i;
        _mm_storeu_si128(q, dcba);
        _mm_storeu_si128(q.add(1), hgef);
    }
}

fn compress(h: &mut [u32; 8], chunk: &[u8]) {
    let mut w = [0u32; 64];
    for i in 0..16 {
        w[i] = u32::from_be_bytes([chunk[4 * i], chunk[4 * i + 1], chunk[4 * i + 2], chunk[4 * i + 3]]);
    }
    for i in 16..64 {
        let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
        let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16].wrapping_add(s0).wrapping_add(w[i - 7]).wrapping_add(s1);
    }
    let mut v = *h;
    for i in 0..64 {
        let s1 = v[4].rotate_right(6) ^ v[4].rotate_right(11) ^ v[4].rotate_right(25);
        let ch = (v[4] & v[5]) ^ (!v[4] & v[6]);
        let t1 = v[7].wrapping_add(s1).wrapping_add(ch).wrapping_add(K[i]).wrapping_add(w[i]);
        let s0 = v[0].rotate_right(2) ^ v[0].rotate_right(13) ^ v[0].rotate_right(22);
        let maj = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
        let t2 = s0.wrapping_add(maj);
        v = [t1.wrapping_add(t2), v[0], v[1], v[2], v[3].wrapping_add(t1), v[4], v[5], v[6]];
    }
    for (a, b) in h.iter_mut().zip(v) {
        *a = a.wrapping_add(b);
    }
}

#[cfg(test)]
mod tests {
    use super::{hex, Sha256};

    #[test]
    fn streaming_equals_one_shot() {
        let data: Vec<u8> = (0..1000u32).map(|i| (i * 7 + 3) as u8).collect();
        for pieces in [1usize, 3, 63, 64, 65, 200] {
            let mut h = Sha256::new();
            for c in data.chunks(pieces) {
                h.update(c);
            }
            assert_eq!(h.hex(), hex(&data), "pieces of {pieces}");
        }
        assert_eq!(Sha256::new().hex(), hex(b""));
    }

    #[test]
    fn the_sha_extensions_give_the_scalar_digest() {
        // every length from 0 to 300 bytes and a few large ones, through the run-time choice
        // (the extensions where the CPU has them) against the scalar compression function
        let mut x = 0x2545f4914f6cdd1du64;
        let data: Vec<u8> = (0..(1 << 20) + 77)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                x as u8
            })
            .collect();
        let scalar = |d: &[u8]| -> String {
            let mut h: [u32; 8] = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19];
            let whole = d.len() / 64 * 64;
            for c in d[..whole].chunks_exact(64) {
                super::compress(&mut h, c);
            }
            let mut tail = d[whole..].to_vec();
            tail.push(0x80);
            while tail.len() % 64 != 56 {
                tail.push(0);
            }
            tail.extend_from_slice(&((d.len() as u64) * 8).to_be_bytes());
            for c in tail.chunks_exact(64) {
                super::compress(&mut h, c);
            }
            h.iter().map(|v| format!("{v:08x}")).collect()
        };
        for n in (0..300).chain([4096, 65_537, data.len()]) {
            assert_eq!(hex(&data[..n]), scalar(&data[..n]), "{n} bytes");
        }
        #[cfg(target_arch = "x86_64")]
        if super::shani::available() {
            let mut a = [1u32, 2, 3, 4, 5, 6, 7, 8];
            let mut b = a;
            unsafe { super::shani::compress_blocks(&mut a, &data[..64 * 1000]) };
            for c in data[..64 * 1000].chunks_exact(64) {
                super::compress(&mut b, c);
            }
            assert_eq!(a, b);
        }
    }

    #[test]
    fn a_file_s_digest() {
        let dir = crate::test_dir("sha256_file");
        let p = dir.join("f.bin");
        let data: Vec<u8> = (0..(5u32 << 20) + 3).map(|i| (i * 31 + 7) as u8).collect();
        std::fs::write(&p, &data).unwrap();
        assert_eq!(super::file(&p).unwrap(), (data.len() as u64, hex(&data)));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn fips_180_examples() {
        assert_eq!(hex(b""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
        assert_eq!(hex(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
        assert_eq!(
            hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
        // a million 'a': many whole blocks
        assert_eq!(hex(&vec![b'a'; 1_000_000]), "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
        // tails of 55, 56 and 64 bytes (one or two padding blocks)
        for n in [55usize, 56, 63, 64, 65] {
            assert_eq!(hex(&vec![0u8; n]).len(), 64);
        }
    }
}
