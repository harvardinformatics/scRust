//! gzip decompression (RFC 1952 members holding RFC 1951 DEFLATE streams), in the standard
//! library only, for the 10x Matrix Market files (`matrix.mtx.gz`, `barcodes.tsv.gz`,
//! `features.tsv.gz`). Concatenated members are decoded in turn; each member's CRC-32 and length
//! are checked.
//!
//! Huffman codes are decoded through a table indexed by the next 10 bits, falling back to a
//! canonical bit-by-bit decode for longer codes.

use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GzError(pub String);

impl fmt::Display for GzError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "gzip: {}", self.0)
    }
}

impl std::error::Error for GzError {}

fn err<T>(m: &str) -> Result<T, GzError> {
    Err(GzError(m.to_string()))
}

/// True if `bytes` starts with the gzip magic number.
pub fn is_gzip(bytes: &[u8]) -> bool {
    bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b
}

/// Decompresses every gzip member in `bytes`.
pub fn gunzip(bytes: &[u8]) -> Result<Vec<u8>, GzError> {
    let mut out = Vec::with_capacity(bytes.len() * 4);
    let mut pos = 0;
    while pos < bytes.len() {
        // trailing zero padding after the last member is tolerated, as gzip does
        if bytes[pos..].iter().all(|&b| b == 0) && pos > 0 {
            break;
        }
        pos = member(bytes, pos, &mut out)?;
    }
    Ok(out)
}

fn member(b: &[u8], mut pos: usize, out: &mut Vec<u8>) -> Result<usize, GzError> {
    if b.len() < pos + 18 || !is_gzip(&b[pos..]) {
        return err("not a gzip member");
    }
    if b[pos + 2] != 8 {
        return err("compression method is not deflate");
    }
    let flags = b[pos + 3];
    pos += 10;
    if flags & 0x04 != 0 {
        // FEXTRA
        let xlen = u16::from_le_bytes([b[pos], b[pos + 1]]) as usize;
        pos += 2 + xlen;
    }
    for bit in [0x08u8, 0x10] {
        // FNAME, FCOMMENT: zero-terminated
        if flags & bit != 0 {
            while *b.get(pos).ok_or(GzError("truncated header".into()))? != 0 {
                pos += 1;
            }
            pos += 1;
        }
    }
    if flags & 0x02 != 0 {
        pos += 2; // FHCRC
    }
    let start = out.len();
    let mut r = BitReader { b, pos, bit: 0, nbits: 0 };
    inflate(&mut r, out)?;
    let mut pos = r.byte_pos();
    let trailer = b.get(pos..pos + 8).ok_or(GzError("truncated trailer".into()))?;
    let crc = u32::from_le_bytes(trailer[..4].try_into().unwrap());
    let isize = u32::from_le_bytes(trailer[4..].try_into().unwrap());
    if crc32(&out[start..]) != crc {
        return err("CRC-32 mismatch");
    }
    if (out.len() - start) as u32 != isize {
        return err("length mismatch");
    }
    pos += 8;
    Ok(pos)
}

struct BitReader<'a> {
    b: &'a [u8],
    pos: usize,
    bit: u64,
    nbits: u32,
}

impl BitReader<'_> {
    fn refill(&mut self) {
        while self.nbits <= 56 {
            let byte = if self.pos < self.b.len() { self.b[self.pos] } else { 0 };
            self.pos += 1;
            self.bit |= (byte as u64) << self.nbits;
            self.nbits += 8;
        }
    }

    fn bits(&mut self, n: u32) -> Result<u32, GzError> {
        if n == 0 {
            return Ok(0);
        }
        if self.nbits < n {
            self.refill();
        }
        if self.pos > self.b.len() + 8 {
            return err("unexpected end of data");
        }
        let v = (self.bit & ((1u64 << n) - 1)) as u32;
        self.bit >>= n;
        self.nbits -= n;
        Ok(v)
    }

    fn peek(&mut self, n: u32) -> u32 {
        if self.nbits < n {
            self.refill();
        }
        (self.bit & ((1u64 << n) - 1)) as u32
    }

    fn consume(&mut self, n: u32) {
        self.bit >>= n;
        self.nbits -= n;
    }

    /// Drops the bits up to the next byte boundary.
    fn align(&mut self) {
        let r = self.nbits % 8;
        self.consume(r);
    }

    /// The position of the next unread byte (after `align`).
    fn byte_pos(&self) -> usize {
        self.pos - (self.nbits / 8) as usize
    }
}

const FAST_BITS: u32 = 10;

/// A canonical Huffman code: counts per length, symbols in code order, and a lookup table for
/// codes of up to `FAST_BITS` bits (entry: symbol << 4 | length, or 0 when longer).
struct Huffman {
    count: [u16; 16],
    symbol: Vec<u16>,
    fast: Vec<u16>,
}

impl Huffman {
    #[allow(clippy::needless_range_loop)]
    fn new(lengths: &[u8]) -> Result<Huffman, GzError> {
        let mut count = [0u16; 16];
        for &l in lengths {
            count[l as usize] += 1;
        }
        count[0] = 0;
        let mut offs = [0u16; 16];
        for len in 1..16 {
            offs[len] = offs[len - 1] + count[len - 1];
        }
        let mut symbol = vec![0u16; lengths.len()];
        for (s, &l) in lengths.iter().enumerate() {
            if l != 0 {
                symbol[offs[l as usize] as usize] = s as u16;
                offs[l as usize] += 1;
            }
        }
        // the fast table: walk codes in canonical order
        let mut fast = vec![0u16; 1 << FAST_BITS];
        let (mut first, mut index) = (0u32, 0usize);
        for len in 1..=FAST_BITS as usize {
            let c = count[len] as u32;
            for k in 0..c {
                let codeval = first + k; // canonical code, MSB-first
                let sym = symbol[index + k as usize];
                // DEFLATE sends codes MSB-first within the LSB-first bit stream: reverse
                let mut rev = 0u32;
                for b in 0..len {
                    rev |= ((codeval >> b) & 1) << (len - 1 - b);
                }
                let step = 1u32 << len;
                let mut f = rev;
                while f < (1 << FAST_BITS) {
                    fast[f as usize] = (sym << 4) | len as u16;
                    f += step;
                }
            }
            index += c as usize;
            first = (first + c) << 1;
        }
        Ok(Huffman { count, symbol, fast })
    }

    fn decode(&self, r: &mut BitReader) -> Result<u16, GzError> {
        let e = self.fast[r.peek(FAST_BITS) as usize];
        if e != 0 {
            r.consume((e & 15) as u32);
            return Ok(e >> 4);
        }
        // canonical decode, bit by bit (puff.c's method)
        let (mut code, mut first, mut index) = (0i32, 0i32, 0i32);
        for len in 1..16 {
            code |= r.bits(1)? as i32;
            let count = self.count[len] as i32;
            if code - count < first {
                return Ok(self.symbol[(index + (code - first)) as usize]);
            }
            index += count;
            first += count;
            first <<= 1;
            code <<= 1;
        }
        err("invalid Huffman code")
    }
}

const LBASE: [u16; 29] = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
const LEXT: [u8; 29] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
const DBASE: [u16; 30] = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577];
const DEXT: [u8; 30] = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];

fn inflate(r: &mut BitReader, out: &mut Vec<u8>) -> Result<(), GzError> {
    let start = out.len();
    loop {
        let last = r.bits(1)?;
        match r.bits(2)? {
            0 => {
                r.align();
                let len = r.bits(16)?;
                let nlen = r.bits(16)?;
                if len != !nlen & 0xffff {
                    return err("stored block length check failed");
                }
                for _ in 0..len {
                    out.push(r.bits(8)? as u8);
                }
            }
            1 => {
                let mut l = [0u8; 288];
                l[..144].fill(8);
                l[144..256].fill(9);
                l[256..280].fill(7);
                l[280..].fill(8);
                let lit = Huffman::new(&l)?;
                let dist = Huffman::new(&[5u8; 30])?;
                codes(r, out, start, &lit, &dist)?;
            }
            2 => {
                let nlen = r.bits(5)? as usize + 257;
                let ndist = r.bits(5)? as usize + 1;
                let ncode = r.bits(4)? as usize + 4;
                const ORDER: [usize; 19] = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
                let mut cl = [0u8; 19];
                for &o in ORDER.iter().take(ncode) {
                    cl[o] = r.bits(3)? as u8;
                }
                let clh = Huffman::new(&cl)?;
                let mut lengths = vec![0u8; nlen + ndist];
                let mut i = 0;
                while i < nlen + ndist {
                    let sym = clh.decode(r)?;
                    match sym {
                        0..=15 => {
                            lengths[i] = sym as u8;
                            i += 1;
                        }
                        16 => {
                            if i == 0 {
                                return err("repeat with no previous length");
                            }
                            let prev = lengths[i - 1];
                            for _ in 0..3 + r.bits(2)? {
                                *lengths.get_mut(i).ok_or(GzError("too many lengths".into()))? = prev;
                                i += 1;
                            }
                        }
                        17 | 18 => {
                            let n = if sym == 17 { 3 + r.bits(3)? } else { 11 + r.bits(7)? };
                            for _ in 0..n {
                                *lengths.get_mut(i).ok_or(GzError("too many lengths".into()))? = 0;
                                i += 1;
                            }
                        }
                        _ => return err("bad code length symbol"),
                    }
                }
                let lit = Huffman::new(&lengths[..nlen])?;
                let dist = Huffman::new(&lengths[nlen..])?;
                codes(r, out, start, &lit, &dist)?;
            }
            _ => return err("invalid block type"),
        }
        if last == 1 {
            r.align();
            return Ok(());
        }
    }
}

fn codes(r: &mut BitReader, out: &mut Vec<u8>, start: usize, lit: &Huffman, dist: &Huffman) -> Result<(), GzError> {
    loop {
        let sym = lit.decode(r)? as usize;
        if sym < 256 {
            out.push(sym as u8);
        } else if sym == 256 {
            return Ok(());
        } else {
            let s = sym - 257;
            if s >= 29 {
                return err("bad length symbol");
            }
            let len = LBASE[s] as usize + r.bits(LEXT[s] as u32)? as usize;
            let ds = dist.decode(r)? as usize;
            if ds >= 30 {
                return err("bad distance symbol");
            }
            let d = DBASE[ds] as usize + r.bits(DEXT[ds] as u32)? as usize;
            if d > out.len() - start {
                return err("distance too far back");
            }
            let from = out.len() - d;
            if d >= len {
                out.extend_from_within(from..from + len);
            } else {
                for k in 0..len {
                    let c = out[from + k];
                    out.push(c);
                }
            }
        }
    }
}

/// CRC-32 (IEEE 802.3, reflected, as gzip uses it).
pub fn crc32(data: &[u8]) -> u32 {
    static TABLE: std::sync::OnceLock<[u32; 256]> = std::sync::OnceLock::new();
    let t = TABLE.get_or_init(|| {
        let mut t = [0u32; 256];
        for (i, e) in t.iter_mut().enumerate() {
            let mut c = i as u32;
            for _ in 0..8 {
                c = if c & 1 != 0 { 0xedb8_8320 ^ (c >> 1) } else { c >> 1 };
            }
            *e = c;
        }
        t
    });
    let mut c = !0u32;
    for &b in data {
        c = t[((c ^ b as u32) & 0xff) as usize] ^ (c >> 8);
    }
    !c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crc_of_known_text() {
        assert_eq!(crc32(b"123456789"), 0xcbf4_3926);
    }

    /// Decompresses gzip files of the committed fixtures and compares with the system's gzip.
    #[test]
    fn matches_system_gzip() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/doubletfinder");
        let mut checked = 0;
        for obj in ["soupx_tenx_test", "cellbender_fromraw_test"] {
            for f in ["matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"] {
                let path = root.join(obj).join("input").join(f);
                let Ok(bytes) = std::fs::read(&path) else { continue };
                let got = gunzip(&bytes).unwrap();
                let Ok(sys) = std::process::Command::new("gzip").arg("-dc").arg(&path).output() else {
                    continue;
                };
                assert!(sys.status.success());
                assert_eq!(got, sys.stdout, "{}", path.display());
                checked += 1;
            }
        }
        assert!(checked > 0 || std::process::Command::new("gzip").arg("--version").output().is_err());
    }
}
