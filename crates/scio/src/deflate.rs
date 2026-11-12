//! gzip compression (RFC 1952 members holding RFC 1951 DEFLATE streams), in the standard library
//! only, for the 10x Matrix Market directories the tools write ([`crate::tenx::write_columns`]).
//! [`gzip`] compresses a whole buffer. [`Chunk`] and [`GzWriter`] compress a stream in
//! independent pieces, on several threads.
//!
//! - **Matching:** LZ77 over a 32 KiB window. Hash chains on 3-byte prefixes, at most
//!   [`MAX_CHAIN`] candidates per position, and one step of lazy evaluation, as zlib does.
//! - **Coding:** each block of at most [`BLOCK_SYMBOLS`] symbols gets its own dynamic Huffman
//!   codes. The package-merge algorithm limits them to 15 bits (7 for the code-length code).
//! - **Pieces:** [`Chunk::compress`] codes one piece of the stream without references into the
//!   previous piece. It ends the piece on a byte boundary with an empty stored block, as zlib's
//!   `Z_SYNC_FLUSH` does. [`GzWriter`] concatenates the pieces in order, combines their CRC-32s,
//!   and closes the stream with an empty final block. The bytes depend on the data and on where
//!   the pieces are cut, never on the number of threads.
//!
//! The result is one gzip member that every gzip reader decodes (`gzip -d`, zlib, R's `gzfile`,
//! Python's `gzip`, [`crate::gz::gunzip`]). It is not byte-identical to zlib's output, and need
//! not be.

use crate::gz::crc32;
use std::io::{self, Write};
use std::sync::OnceLock;

const WINDOW: usize = 1 << 15;
/// zlib's MAX_DIST: matches stay within the window minus the look-ahead, so that a chain never
/// reaches a ring slot that the current position has overwritten.
const MAX_DIST: usize = WINDOW - (MAX_MATCH + MIN_MATCH + 1);
const MIN_MATCH: usize = 3;
const MAX_MATCH: usize = 258;
const HASH_BITS: u32 = 15;
const NONE: u32 = u32::MAX;
/// Candidates examined per position.
pub const MAX_CHAIN: usize = 16;
/// A match at least this long ends the search.
const NICE_MATCH: usize = 128;
/// A match at least this long is taken without looking for a longer one at the next position.
const MAX_LAZY: usize = 16;
/// A length-3 match farther than this costs more than three literals (zlib's TOO_FAR).
const TOO_FAR: usize = 4096;
/// Symbols per block; each block has its own codes.
pub const BLOCK_SYMBOLS: usize = 1 << 15;

/// A token: a literal byte, or `MATCH | length << 16 | distance`.
const MATCH: u32 = 1 << 31;

const LBASE: [u16; 29] = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
const LEXT: [u8; 29] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
const DBASE: [u16; 30] = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577];
const DEXT: [u8; 30] = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];
/// The order in which the code-length code's lengths are sent.
const CL_ORDER: [usize; 19] = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];

/// The gzip header: deflate, no flags, no time stamp (so the output depends on the data only),
/// OS "unknown".
const HEADER: [u8; 10] = [0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 255];
/// An empty final block with fixed codes: BFINAL = 1, BTYPE = 01, end of block (seven 0 bits).
const FINAL_BLOCK: [u8; 2] = [0x03, 0x00];

/// Length code (0..=28) of every match length 3..=258; 258 has its own code, 285.
fn length_code(len: usize) -> usize {
    static TABLE: OnceLock<[u8; MAX_MATCH + 1]> = OnceLock::new();
    let t = TABLE.get_or_init(|| {
        let mut t = [0u8; MAX_MATCH + 1];
        for (l, e) in t.iter_mut().enumerate().skip(MIN_MATCH) {
            *e = LBASE.iter().rposition(|&b| b as usize <= l).unwrap() as u8;
        }
        t
    });
    t[len] as usize
}

/// Distance code (0..=29) of a distance 1..=32768.
fn dist_code(dist: usize) -> usize {
    let d = (dist - 1) as u32;
    if d < 4 {
        d as usize
    } else {
        let n = 31 - d.leading_zeros();
        (2 * n + ((d >> (n - 1)) & 1)) as usize
    }
}

/// Bits, least significant first, as DEFLATE packs them.
struct Bits {
    out: Vec<u8>,
    acc: u64,
    n: u32,
}

impl Bits {
    fn new(capacity: usize) -> Bits {
        Bits { out: Vec::with_capacity(capacity), acc: 0, n: 0 }
    }

    /// The `len` low bits of `v` (len <= 32).
    #[inline]
    fn put(&mut self, v: u32, len: u32) {
        debug_assert!(len <= 32 && (len == 32 || v >> len == 0) && self.n < 32);
        self.acc |= (v as u64) << self.n;
        self.n += len;
        if self.n >= 32 {
            self.out.extend_from_slice(&(self.acc as u32).to_le_bytes());
            self.acc >>= 32;
            self.n -= 32;
        }
    }

    /// Pads with 0 bits to a byte boundary.
    fn align(&mut self) {
        while self.n > 0 {
            self.out.push(self.acc as u8);
            self.acc >>= 8;
            self.n = self.n.saturating_sub(8);
        }
        self.acc = 0;
    }
}

/// Code lengths for `freq` (0 for unused symbols), at most `limit` bits, by package-merge
/// (Larmore and Hirschberg): optimal among codes of bounded length. With at least two used
/// symbols the code is complete (Kraft sum exactly 1); a single used symbol gets length 1.
fn code_lengths(freq: &[u32], limit: usize) -> Vec<u8> {
    let mut lens = vec![0u8; freq.len()];
    let mut syms: Vec<usize> = (0..freq.len()).filter(|&s| freq[s] > 0).collect();
    let n = syms.len();
    if n == 0 {
        return lens;
    }
    if n == 1 {
        lens[syms[0]] = 1;
        return lens;
    }
    assert!(n <= 1 << limit, "{n} symbols cannot have codes of at most {limit} bits");
    syms.sort_by_key(|&s| (freq[s], s));
    let leaves: Vec<u64> = syms.iter().map(|&s| freq[s] as u64).collect();
    // levels[0] holds the leaves alone (the deepest level); each next level merges the leaves
    // with the packages (sums of consecutive pairs) of the level below it, a leaf first on ties.
    let mut levels: Vec<Vec<(u64, bool)>> = Vec::with_capacity(limit);
    levels.push(leaves.iter().map(|&w| (w, true)).collect());
    for _ in 1..limit {
        let below = levels.last().unwrap();
        let packages: Vec<u64> = below.as_chunks::<2>().0.iter().map(|[a, b]| a.0 + b.0).collect();
        let mut merged = Vec::with_capacity(n + packages.len());
        let (mut i, mut j) = (0, 0);
        while i < n || j < packages.len() {
            if j == packages.len() || (i < n && leaves[i] <= packages[j]) {
                merged.push((leaves[i], true));
                i += 1;
            } else {
                merged.push((packages[j], false));
                j += 1;
            }
        }
        levels.push(merged);
    }
    // the 2n - 2 lightest items of the top level; each package selected at one level selects
    // its two items at the level below. A symbol's length is the number of levels at which its
    // leaf is selected; the selected leaves of a level are the lightest ones.
    let mut take = 2 * n - 2;
    for level in levels.iter().rev() {
        if take == 0 {
            break;
        }
        let selected = &level[..take];
        let n_leaves = selected.iter().filter(|e| e.1).count();
        for &s in &syms[..n_leaves] {
            lens[s] += 1;
        }
        take = 2 * (take - n_leaves);
    }
    lens
}

/// Canonical codes (RFC 1951 3.2.2) for `lens`, bit-reversed so that [`Bits::put`] sends their
/// most significant bit first.
fn canonical_codes(lens: &[u8]) -> Vec<u32> {
    let mut count = [0u32; 16];
    for &l in lens {
        if l > 0 {
            count[l as usize] += 1;
        }
    }
    let mut next = [0u32; 16];
    let mut code = 0u32;
    for bits in 1..16 {
        code = (code + count[bits - 1]) << 1;
        next[bits] = code;
    }
    lens.iter()
        .map(|&l| {
            if l == 0 {
                return 0;
            }
            let c = next[l as usize];
            next[l as usize] += 1;
            c.reverse_bits() >> (32 - l as u32)
        })
        .collect()
}

/// At least two used symbols, as zlib ensures: a lone code would be incomplete, which inflaters
/// reject for the code-length code. Unused symbols with the lowest numbers get frequency 1.
fn at_least_two(freq: &mut [u32]) {
    let mut used = freq.iter().filter(|&&f| f > 0).count();
    for f in freq.iter_mut() {
        if used >= 2 {
            break;
        }
        if *f == 0 {
            *f = 1;
            used += 1;
        }
    }
}

/// Run-length codes for a sequence of code lengths: (symbol 0..=18, extra bits value).
fn rle_lengths(seq: &[u8]) -> Vec<(u8, u8)> {
    let mut out = Vec::with_capacity(seq.len());
    let mut i = 0;
    while i < seq.len() {
        let l = seq[i];
        let mut run = 1;
        while i + run < seq.len() && seq[i + run] == l {
            run += 1;
        }
        let mut r = run;
        if l == 0 {
            while r >= 11 {
                let k = r.min(138);
                out.push((18, (k - 11) as u8));
                r -= k;
            }
            if r >= 3 {
                out.push((17, (r - 3) as u8));
                r = 0;
            }
        } else {
            out.push((l, 0));
            r -= 1;
            while r >= 3 {
                let k = r.min(6);
                out.push((16, (k - 3) as u8));
                r -= k;
            }
        }
        for _ in 0..r {
            out.push((l, 0));
        }
        i += run;
    }
    out
}

/// One block with dynamic Huffman codes (BTYPE = 10), never the final one.
fn write_block(bits: &mut Bits, tokens: &[u32]) {
    let mut lfreq = [0u32; 286];
    let mut dfreq = [0u32; 30];
    for &t in tokens {
        if t & MATCH != 0 {
            lfreq[257 + length_code(((t >> 16) & 0x1ff) as usize)] += 1;
            dfreq[dist_code((t & 0xffff) as usize)] += 1;
        } else {
            lfreq[t as usize] += 1;
        }
    }
    lfreq[256] = 1;
    at_least_two(&mut lfreq);
    at_least_two(&mut dfreq);
    let llen = code_lengths(&lfreq, 15);
    let dlen = code_lengths(&dfreq, 15);
    let hlit = (llen.iter().rposition(|&l| l > 0).unwrap() + 1).max(257);
    let hdist = (dlen.iter().rposition(|&l| l > 0).unwrap() + 1).max(1);
    let seq: Vec<u8> = llen[..hlit].iter().chain(&dlen[..hdist]).copied().collect();
    let rle = rle_lengths(&seq);
    let mut cfreq = [0u32; 19];
    for &(s, _) in &rle {
        cfreq[s as usize] += 1;
    }
    at_least_two(&mut cfreq);
    let clen = code_lengths(&cfreq, 7);
    let hclen = (CL_ORDER.iter().rposition(|&s| clen[s] > 0).unwrap() + 1).max(4);

    bits.put(0, 1); // BFINAL
    bits.put(2, 2); // BTYPE = 10
    bits.put((hlit - 257) as u32, 5);
    bits.put((hdist - 1) as u32, 5);
    bits.put((hclen - 4) as u32, 4);
    for &s in &CL_ORDER[..hclen] {
        bits.put(clen[s] as u32, 3);
    }
    let ccode = canonical_codes(&clen);
    for &(s, e) in &rle {
        bits.put(ccode[s as usize], clen[s as usize] as u32);
        match s {
            16 => bits.put(e as u32, 2),
            17 => bits.put(e as u32, 3),
            18 => bits.put(e as u32, 7),
            _ => {}
        }
    }
    let lcode = canonical_codes(&llen);
    let dcode = canonical_codes(&dlen);
    for &t in tokens {
        if t & MATCH != 0 {
            let len = ((t >> 16) & 0x1ff) as usize;
            let dist = (t & 0xffff) as usize;
            let lc = length_code(len);
            bits.put(lcode[257 + lc], llen[257 + lc] as u32);
            bits.put((len - LBASE[lc] as usize) as u32, LEXT[lc] as u32);
            let dc = dist_code(dist);
            bits.put(dcode[dc], dlen[dc] as u32);
            bits.put((dist - DBASE[dc] as usize) as u32, DEXT[dc] as u32);
        } else {
            bits.put(lcode[t as usize], llen[t as usize] as u32);
        }
    }
    bits.put(lcode[256], llen[256] as u32);
}

/// LZ77 state for one piece: hash heads and a ring of chain links over the window.
struct Matcher<'a> {
    data: &'a [u8],
    head: Vec<u32>,
    prev: Vec<u32>,
}

impl<'a> Matcher<'a> {
    fn new(data: &'a [u8]) -> Matcher<'a> {
        Matcher { data, head: vec![NONE; 1 << HASH_BITS], prev: vec![NONE; WINDOW] }
    }

    /// Inserts position `pos` into its hash chain; returns the chain's previous head (NONE when
    /// fewer than 3 bytes remain).
    #[inline]
    fn insert(&mut self, pos: usize) -> u32 {
        let d = self.data;
        if pos + MIN_MATCH > d.len() {
            return NONE;
        }
        let v = d[pos] as u32 | (d[pos + 1] as u32) << 8 | (d[pos + 2] as u32) << 16;
        let h = (v.wrapping_mul(0x9e37_79b1) >> (32 - HASH_BITS)) as usize;
        let old = self.head[h];
        self.prev[pos & (WINDOW - 1)] = old;
        self.head[h] = pos as u32;
        old
    }

    /// The longest match for `pos` along the chain from `cand` (the nearest one among equals):
    /// (length, distance), length 0 if none reaches 3.
    #[inline]
    fn longest(&self, pos: usize, mut cand: u32) -> (usize, usize) {
        let d = self.data;
        let max_len = MAX_MATCH.min(d.len() - pos);
        if max_len < MIN_MATCH {
            return (0, 0);
        }
        let min_pos = pos.saturating_sub(MAX_DIST);
        let (mut best, mut dist) = (MIN_MATCH - 1, 0);
        let mut chain = MAX_CHAIN;
        while cand != NONE && cand as usize >= min_pos && chain > 0 {
            let c = cand as usize;
            if d[c + best] == d[pos + best] && d[c] == d[pos] && d[c + 1] == d[pos + 1] {
                let len = match_len(d, c, pos, max_len);
                if len > best {
                    best = len;
                    dist = pos - c;
                    if len >= max_len.min(NICE_MATCH) {
                        break;
                    }
                }
            }
            cand = self.prev[c & (WINDOW - 1)];
            chain -= 1;
        }
        if best < MIN_MATCH || (best == MIN_MATCH && dist > TOO_FAR) {
            (0, 0)
        } else {
            (best, dist)
        }
    }
}

/// Length of the common prefix of `d[a..]` and `d[b..]`, at most `max` (a < b, b + max <= len).
#[inline]
fn match_len(d: &[u8], a: usize, b: usize, max: usize) -> usize {
    let mut n = 0;
    while n + 8 <= max {
        let x = u64::from_le_bytes(d[a + n..a + n + 8].try_into().unwrap());
        let y = u64::from_le_bytes(d[b + n..b + n + 8].try_into().unwrap());
        let diff = x ^ y;
        if diff != 0 {
            return n + (diff.trailing_zeros() / 8) as usize;
        }
        n += 8;
    }
    while n < max && d[a + n] == d[b + n] {
        n += 1;
    }
    n
}

/// One piece of a gzip stream: DEFLATE blocks ending on a byte boundary (never the final
/// block), the piece's CRC-32 and its uncompressed length.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chunk {
    pub deflate: Vec<u8>,
    pub crc: u32,
    pub len: u64,
}

impl Chunk {
    /// Compresses `data` on its own (no references outside it). Empty data gives an empty
    /// piece.
    pub fn compress(data: &[u8]) -> Chunk {
        let crc = crc32(data);
        if data.is_empty() {
            return Chunk { deflate: Vec::new(), crc, len: 0 };
        }
        fn push(tokens: &mut Vec<u32>, bits: &mut Bits, t: u32) {
            tokens.push(t);
            if tokens.len() == BLOCK_SYMBOLS {
                write_block(bits, tokens);
                tokens.clear();
            }
        }
        let mut bits = Bits::new(data.len() / 3 + 64);
        let mut tokens: Vec<u32> = Vec::with_capacity(BLOCK_SYMBOLS);
        let mut m = Matcher::new(data);
        let n = data.len();
        // zlib's deflate_slow: the match found at pos - 1 is kept pending until pos shows
        // whether a longer one starts there
        let (mut pending, mut prev_len, mut prev_dist) = (false, 0usize, 0usize);
        let mut pos = 0;
        while pos < n {
            let head = m.insert(pos);
            let (len, dist) = if head != NONE && prev_len < MAX_LAZY { m.longest(pos, head) } else { (0, 0) };
            if pending && prev_len >= MIN_MATCH && len <= prev_len {
                push(&mut tokens, &mut bits, MATCH | (prev_len as u32) << 16 | prev_dist as u32);
                let end = pos - 1 + prev_len;
                for k in pos + 1..end {
                    m.insert(k);
                }
                pos = end;
                pending = false;
                prev_len = 0;
                continue;
            }
            if pending {
                push(&mut tokens, &mut bits, data[pos - 1] as u32);
            }
            pending = true;
            prev_len = len;
            prev_dist = dist;
            pos += 1;
        }
        if pending {
            // a match cannot start at the last byte
            push(&mut tokens, &mut bits, data[n - 1] as u32);
        }
        if !tokens.is_empty() {
            write_block(&mut bits, &tokens);
        }
        // Z_SYNC_FLUSH: an empty stored block (BFINAL = 0, BTYPE = 00), then LEN = 0, NLEN = !0
        bits.put(0, 3);
        bits.align();
        bits.out.extend_from_slice(&[0, 0, 0xff, 0xff]);
        Chunk { deflate: bits.out, crc, len: n as u64 }
    }
}

/// A gzip member written piece by piece: the header first, then each [`Chunk`] in order, then
/// [`finish`](Self::finish) for the final block and the trailer (CRC-32 and length).
pub struct GzWriter<W: Write> {
    inner: W,
    crc: u32,
    len: u64,
}

impl<W: Write> GzWriter<W> {
    pub fn new(mut inner: W) -> io::Result<GzWriter<W>> {
        inner.write_all(&HEADER)?;
        Ok(GzWriter { inner, crc: 0, len: 0 })
    }

    pub fn write_chunk(&mut self, chunk: &Chunk) -> io::Result<()> {
        self.inner.write_all(&chunk.deflate)?;
        self.crc = crc32_combine(self.crc, chunk.crc, chunk.len);
        self.len += chunk.len;
        Ok(())
    }

    /// Ends the stream and returns the inner writer, flushed.
    pub fn finish(mut self) -> io::Result<W> {
        self.inner.write_all(&FINAL_BLOCK)?;
        self.inner.write_all(&self.crc.to_le_bytes())?;
        self.inner.write_all(&(self.len as u32).to_le_bytes())?;
        self.inner.flush()?;
        Ok(self.inner)
    }
}

/// `data` as one gzip member.
pub fn gzip(data: &[u8]) -> Vec<u8> {
    let mut w = GzWriter::new(Vec::with_capacity(data.len() / 3 + 32)).unwrap();
    w.write_chunk(&Chunk::compress(data)).unwrap();
    w.finish().unwrap()
}

const POLY: u32 = 0xedb8_8320;

/// `a * b mod p(x)`, polynomials over GF(2) in the reflected bit order of the CRC (zlib's
/// `multmodp`).
fn multmodp(a: u32, mut b: u32) -> u32 {
    let mut p = 0u32;
    let mut m = 1u32 << 31;
    while m != 0 {
        if a & m != 0 {
            p ^= b;
        }
        b = if b & 1 != 0 { (b >> 1) ^ POLY } else { b >> 1 };
        m >>= 1;
    }
    p
}

/// `x^(n * 2^k) mod p(x)` (zlib's `x2nmodp`).
fn x2nmodp(mut n: u64, mut k: u32) -> u32 {
    static TABLE: OnceLock<[u32; 32]> = OnceLock::new();
    let t = TABLE.get_or_init(|| {
        let mut t = [0u32; 32];
        let mut p = 1u32 << 30; // x^1
        t[0] = p;
        for e in t.iter_mut().skip(1) {
            p = multmodp(p, p);
            *e = p;
        }
        t
    });
    let mut p = 1u32 << 31; // x^0
    while n != 0 {
        if n & 1 != 0 {
            p = multmodp(t[(k & 31) as usize], p);
        }
        n >>= 1;
        k += 1;
    }
    p
}

/// The CRC-32 of `A` followed by `B`, from `crc32(A)`, `crc32(B)` and B's length (zlib's
/// `crc32_combine`).
pub fn crc32_combine(crc_a: u32, crc_b: u32, len_b: u64) -> u32 {
    multmodp(x2nmodp(len_b, 3), crc_a) ^ crc_b
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gz::gunzip;

    /// A small deterministic generator (xorshift64*), so the tests need no crate.
    struct Gen(u64);
    impl Gen {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// Matrix Market text like the counts': increasing rows within columns, small values.
    fn mtx_like(entries: usize, g: &mut Gen) -> Vec<u8> {
        let mut s = String::from("%%MatrixMarket matrix coordinate integer general\n2000 300 0\n");
        let (mut col, mut row) = (1u64, 0u64);
        for _ in 0..entries {
            row += 1 + g.below(40);
            if row > 2000 {
                col += 1;
                row = 1 + g.below(20);
            }
            let v = 1 + g.below(4) * g.below(4) * g.below(20);
            s.push_str(&format!("{row} {col} {v}\n"));
        }
        s.into_bytes()
    }

    fn round_trip(data: &[u8]) -> Vec<u8> {
        let z = gzip(data);
        assert_eq!(gunzip(&z).unwrap(), data, "{} bytes", data.len());
        z
    }

    #[test]
    fn round_trips() {
        let mut g = Gen(0x9e37_79b9_7f4a_7c15);
        round_trip(b"");
        round_trip(b"a");
        round_trip(b"ab");
        round_trip(b"abc");
        round_trip(b"abcabcabcabcabcabc");
        round_trip(&[0u8; 100_000]); // runs up to the longest match, distance 1
        round_trip(&(0..=255u8).cycle().take(70_000).collect::<Vec<u8>>());
        let random: Vec<u8> = (0..200_000).map(|_| g.next() as u8).collect();
        round_trip(&random);
        // repeats beyond the window, near its edge, and length-3 matches far away
        let mut far = random[..40_000].to_vec();
        far.extend_from_slice(&random[..40_000]);
        far.extend_from_slice(&random[7_000..7_100]);
        round_trip(&far);
        let text = mtx_like(150_000, &mut g);
        let z = round_trip(&text);
        assert!(z.len() * 2 < text.len(), "{} -> {}", text.len(), z.len());
    }

    #[test]
    fn pieces_make_one_stream() {
        let mut g = Gen(42);
        let text = mtx_like(60_000, &mut g);
        let cuts = [0, 1, 10, 4_096, 33_000, 200_000, 512_345, text.len()];
        let mut w = GzWriter::new(Vec::new()).unwrap();
        for c in cuts.windows(2) {
            let (a, b) = (c[0].min(text.len()), c[1].min(text.len()));
            w.write_chunk(&Chunk::compress(&text[a..b])).unwrap();
        }
        let z = w.finish().unwrap();
        assert_eq!(gunzip(&z).unwrap(), text);
        // the same cuts give the same bytes
        let mut w = GzWriter::new(Vec::new()).unwrap();
        for c in cuts.windows(2) {
            let (a, b) = (c[0].min(text.len()), c[1].min(text.len()));
            w.write_chunk(&Chunk::compress(&text[a..b])).unwrap();
        }
        assert_eq!(w.finish().unwrap(), z);
    }

    #[test]
    fn crc_combines() {
        let mut g = Gen(7);
        let data: Vec<u8> = (0..10_000).map(|_| g.next() as u8).collect();
        for cut in [0, 1, 2, 3, 100, 4_999, 9_999, 10_000] {
            let (a, b) = data.split_at(cut);
            assert_eq!(crc32_combine(crc32(a), crc32(b), b.len() as u64), crc32(&data), "cut {cut}");
        }
    }

    #[test]
    fn lengths_are_optimal_and_bounded() {
        // Huffman's lengths when no limit binds
        assert_eq!(code_lengths(&[1, 1, 2, 4], 15), vec![3, 3, 2, 1]);
        assert_eq!(code_lengths(&[0, 5, 0], 15), vec![0, 1, 0]);
        // Fibonacci weights: Huffman would need 30 bits
        let fib: Vec<u32> = (0..31).scan((1u32, 1u32), |s, _| {
            let v = s.0;
            *s = (s.1, s.0 + s.1);
            Some(v)
        }).collect();
        for limit in [7, 15] {
            let freq: Vec<u32> = fib[..19.min(fib.len())].to_vec();
            let lens = code_lengths(&freq, limit);
            assert!(lens.iter().all(|&l| l >= 1 && l as usize <= limit), "{lens:?}");
            let kraft: f64 = lens.iter().map(|&l| 0.5f64.powi(l as i32)).sum();
            assert_eq!(kraft, 1.0, "{lens:?}");
        }
        let lens = code_lengths(&fib, 15);
        assert_eq!(lens.iter().max(), Some(&15));
        assert_eq!(lens.iter().map(|&l| 0.5f64.powi(l as i32)).sum::<f64>(), 1.0);
    }

    #[test]
    fn length_and_distance_codes() {
        assert_eq!((length_code(3), length_code(10), length_code(11), length_code(257), length_code(258)), (0, 7, 8, 27, 28));
        for d in 1..=32768usize {
            let c = dist_code(d);
            assert!(DBASE[c] as usize <= d && d < DBASE[c] as usize + (1 << DEXT[c]), "{d}");
        }
    }

    #[test]
    fn rle_of_code_lengths() {
        let seq: Vec<u8> = [vec![0; 150], vec![8; 9], vec![0; 4], vec![5], vec![0; 2]].concat();
        let rle = rle_lengths(&seq);
        let mut back = Vec::new();
        for &(s, e) in &rle {
            match s {
                16 => {
                    let l = *back.last().unwrap();
                    back.extend(std::iter::repeat_n(l, 3 + e as usize));
                }
                17 => back.extend(std::iter::repeat_n(0, 3 + e as usize)),
                18 => back.extend(std::iter::repeat_n(0, 11 + e as usize)),
                l => back.push(l),
            }
        }
        assert_eq!(back, seq);
    }

    /// The system's gzip, where there is one, reads what this module writes.
    #[test]
    fn system_gzip_reads_it() {
        let mut g = Gen(3);
        let text = mtx_like(20_000, &mut g);
        let mut w = GzWriter::new(Vec::new()).unwrap();
        for piece in text.chunks(50_000) {
            w.write_chunk(&Chunk::compress(piece)).unwrap();
        }
        let z = w.finish().unwrap();
        let path = crate::test_dir("deflate").join("mtx_like.gz");
        std::fs::write(&path, &z).unwrap();
        let Ok(out) = std::process::Command::new("gzip").arg("-dc").arg(&path).output() else {
            return;
        };
        assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
        assert_eq!(out.stdout, text);
        let t = std::process::Command::new("gzip").arg("-t").arg(&path).status().unwrap();
        assert!(t.success());
    }
}
