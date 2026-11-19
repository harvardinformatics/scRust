//! The allocator at in-process hand-offs (A3's in-memory rule, M15).
//!
//! A chain of tools in one process peaks at its largest step's peak only if the memory a step
//! frees goes back to the system before the next step. With glibc's malloc it does not, by
//! default, after a step that freed large buffers. glibc serves a block above its mmap threshold
//! (128 KiB at first) with `mmap`, which returns to the system when freed, and a smaller one from
//! its heaps, which keep freed pages for reuse. Each time an mmapped block is freed, glibc raises
//! the threshold to that block's size, up to 32 MiB. So after a step has freed its large
//! buffers, the next step's medium-sized blocks (128 KiB to 32 MiB) come from the heaps, and their
//! pages stay resident between uses, on top of the next step's own peak. In M15, clustering on the
//! fly before DoubletFinder raised DoubletFinder's peak by up to 16% above M7's runs this way
//! (docs/validation/clustering.md, "Memory and time").
//!
//! [`settle`] puts the allocator back: called between two steps of one process, at the hand-off,
//! it fixes the mmap threshold at its default, 128 KiB (which also ends the adjustment), and
//! returns the heaps' free pages to the system (`malloc_trim(0)`). Call it at in-process hand-off
//! boundaries only, after a step that freed large buffers: never in a command's own path, whose
//! memory profile it would change. It changes no result.

#[cfg(all(target_os = "linux", target_env = "gnu"))]
mod glibc {
    use std::os::raw::c_int;

    extern "C" {
        /// glibc's `mallopt(3)`.
        pub fn mallopt(param: c_int, value: c_int) -> c_int;
        /// glibc's `malloc_trim(3)`: returns the free pages of every arena to the system.
        pub fn malloc_trim(pad: usize) -> c_int;
    }

    /// `M_MMAP_THRESHOLD` (glibc's malloc.h).
    pub const M_MMAP_THRESHOLD: c_int = -3;
}

/// glibc's default mmap threshold, at which [`settle`] fixes it.
pub const MMAP_THRESHOLD: usize = 128 * 1024;

/// At an in-process hand-off, between one step and the next: fixes glibc's mmap threshold at its
/// default, 128 KiB, and returns the heaps' free pages to the system, so that the next step's
/// peak is its own (the module's documentation). On Linux with glibc; elsewhere it does nothing.
/// Returns whether the allocator was settled.
pub fn settle() -> bool {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    {
        // SAFETY: mallopt and malloc_trim only change the allocator's policy and release free
        // pages; both may be called at any time, from any thread
        unsafe {
            let fixed = glibc::mallopt(glibc::M_MMAP_THRESHOLD, MMAP_THRESHOLD as std::os::raw::c_int) == 1;
            glibc::malloc_trim(0);
            fixed
        }
    }
    #[cfg(not(all(target_os = "linux", target_env = "gnu")))]
    {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Settling is allowed at any time, again and again, and leaves the heap usable: what was
    /// allocated before keeps its values, and blocks of every size come and go after it.
    #[test]
    fn settling_keeps_the_heap_usable() {
        let kept: Vec<u64> = (0..100_000).collect();
        {
            // a step that frees a large block (which would raise glibc's threshold)
            let big = vec![1u8; 64 << 20];
            assert_eq!(big[big.len() - 1], 1);
        }
        let settled = settle();
        assert_eq!(settled, cfg!(all(target_os = "linux", target_env = "gnu")));
        for size in [16usize, 4096, MMAP_THRESHOLD - 1, MMAP_THRESHOLD, MMAP_THRESHOLD + 1, 1 << 20, 40 << 20] {
            let v = vec![7u8; size];
            assert!(v.iter().all(|&b| b == 7), "{size}");
        }
        assert_eq!(settle(), settled, "again");
        assert_eq!(kept.iter().sum::<u64>(), 100_000 * 99_999 / 2);
    }
}
