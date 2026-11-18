//! Copying files with their checksums (docs/PLAN.md C3f, milestone M13): the export's `counts/`
//! is the input counts copied byte for byte, and the manifest records each file's SHA-256, which
//! must equal the input's.

use crate::sha256::Sha256;
use std::fs;
use std::io::{self, Read, Write};
use std::path::Path;

/// A file copied: its name in the target directory, its bytes and the SHA-256 of what was read
/// (which is what was written).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Copied {
    pub name: String,
    pub bytes: u64,
    pub sha256: String,
}

fn ctx(e: io::Error, p: &Path) -> io::Error {
    io::Error::new(e.kind(), format!("{}: {e}", p.display()))
}

/// Copies the file `from` to `to`, byte for byte (a symbolic link is followed: its target's bytes
/// are copied), digesting it as it goes. `to` must not exist.
pub fn copy_file(from: &Path, to: &Path) -> io::Result<(u64, String)> {
    let mut src = fs::File::open(from).map_err(|e| ctx(e, from))?;
    let mut dst = fs::OpenOptions::new().write(true).create_new(true).open(to).map_err(|e| ctx(e, to))?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 4 << 20];
    let mut n = 0u64;
    loop {
        let k = src.read(&mut buf).map_err(|e| ctx(e, from))?;
        if k == 0 {
            break;
        }
        h.update(&buf[..k]);
        dst.write_all(&buf[..k]).map_err(|e| ctx(e, to))?;
        n += k as u64;
    }
    dst.flush().map_err(|e| ctx(e, to))?;
    dst.sync_all().map_err(|e| ctx(e, to))?;
    Ok((n, h.hex()))
}

/// Copies `files` (paths of existing files) into the directory `to`, created if needed, each
/// under its own file name, byte for byte, with its checksum. The copy is read back and digested
/// again, and must give the same SHA-256 and size as the source.
pub fn copy_files(files: &[&Path], to: &Path) -> io::Result<Vec<Copied>> {
    fs::create_dir_all(to).map_err(|e| ctx(e, to))?;
    let mut out = Vec::with_capacity(files.len());
    for &from in files {
        let name = from
            .file_name()
            .and_then(|n| n.to_str())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, format!("{}: no file name", from.display())))?
            .to_string();
        let target = to.join(&name);
        let (bytes, sha256) = copy_file(from, &target)?;
        let (b2, s2) = crate::sha256::file(&target)?;
        if b2 != bytes || s2 != sha256 {
            return Err(io::Error::other(format!("{}: the copy of {} differs from it ({b2} bytes, sha256 {s2}; the source {bytes} bytes, {sha256})", target.display(), from.display())));
        }
        out.push(Copied { name, bytes, sha256 });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn copies_byte_for_byte_through_a_link() {
        let dir = crate::test_dir("copy");
        let src = dir.join("src");
        let dst = dir.join("dst");
        let _ = fs::remove_dir_all(&src);
        let _ = fs::remove_dir_all(&dst);
        fs::create_dir_all(&src).unwrap();
        let data: Vec<u8> = (0..(9u32 << 20) + 11).map(|i| (i % 251) as u8).collect();
        fs::write(src.join("matrix.mtx.gz"), &data).unwrap();
        fs::write(src.join("real.tsv.gz"), b"AAAC-1\n").unwrap();
        std::os::unix::fs::symlink("real.tsv.gz", src.join("barcodes.tsv.gz")).unwrap();
        let got = copy_files(&[&src.join("matrix.mtx.gz"), &src.join("barcodes.tsv.gz")], &dst).unwrap();
        assert_eq!(got[0], Copied { name: "matrix.mtx.gz".into(), bytes: data.len() as u64, sha256: crate::sha256::hex(&data) });
        assert_eq!(got[1].name, "barcodes.tsv.gz");
        assert_eq!(fs::read(dst.join("barcodes.tsv.gz")).unwrap(), b"AAAC-1\n");
        assert!(!fs::symlink_metadata(dst.join("barcodes.tsv.gz")).unwrap().file_type().is_symlink());
        assert_eq!(fs::read(dst.join("matrix.mtx.gz")).unwrap(), data);
        // an existing target is never overwritten
        assert!(copy_files(&[&src.join("matrix.mtx.gz")], &dst).is_err());
        fs::remove_dir_all(&dir).unwrap();
    }
}
