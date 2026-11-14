//! The libstdc++ algorithms Annoy's search relies on, statement by statement (GCC's
//! `bits/stl_heap.h` and `bits/stl_algo.h`, unchanged from GCC 14, which built RcppAnnoy, to
//! GCC 16, whose headers the R env carries): `std::priority_queue` (`push_heap`, `pop_heap`) and
//! `std::partial_sort` (`__heap_select`, `__sort_heap`), with `std::less<std::pair<float, int>>`.
//!
//! For pairs of finite floats and distinct ids the order is total, and any correct priority
//! queue or partial sort gives the same result. The copies here make the result identical in
//! every case, including equal pairs and NaN priorities, where the outcome depends on the
//! algorithm's exact moves.

/// `std::pair<float, int>`'s `operator<` (C++17, as RcppAnnoy is compiled: gnu++17):
/// `a.first < b.first || (!(b.first < a.first) && a.second < b.second)`.
#[inline]
pub fn pair_less(a: &(f32, i32), b: &(f32, i32)) -> bool {
    a.0 < b.0 || (!(b.0 < a.0) && a.1 < b.1)
}

/// `std::__push_heap`: moves `value` up from `hole` towards `top`.
#[inline]
fn push_heap_hole<T: Copy>(c: &mut [T], mut hole: usize, top: usize, value: T, less: &impl Fn(&T, &T) -> bool) {
    while hole > top {
        let parent = (hole - 1) / 2;
        if !less(&c[parent], &value) {
            break;
        }
        c[hole] = c[parent];
        hole = parent;
    }
    c[hole] = value;
}

/// `std::__adjust_heap` on `c[..len]`: sifts the hole at `hole` down to a leaf, always taking
/// the larger child (the left one on ties), then pushes `value` back up.
#[inline]
fn adjust_heap<T: Copy>(c: &mut [T], mut hole: usize, len: usize, value: T, less: &impl Fn(&T, &T) -> bool) {
    let top = hole;
    let mut second = hole;
    // the C++ uses signed ptrdiff_t: (len - 1) / 2 and (len - 2) / 2 truncate towards zero
    let len_i = len as isize;
    while (second as isize) < (len_i - 1) / 2 {
        second = 2 * (second + 1);
        if less(&c[second], &c[second - 1]) {
            second -= 1;
        }
        c[hole] = c[second];
        hole = second;
    }
    if len & 1 == 0 && (second as isize) == (len_i - 2) / 2 {
        second = 2 * (second + 1);
        c[hole] = c[second - 1];
        hole = second - 1;
    }
    push_heap_hole(c, hole, top, value, less);
}

/// `std::__make_heap` on `c[..len]`.
fn make_heap<T: Copy>(c: &mut [T], len: usize, less: &impl Fn(&T, &T) -> bool) {
    if len < 2 {
        return;
    }
    let mut parent = (len - 2) / 2;
    loop {
        let value = c[parent];
        adjust_heap(c, parent, len, value, less);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

/// `std::priority_queue<T>` over a `std::vector<T>` with libstdc++'s heap: `top` is the element
/// that is not less than any other.
#[derive(Debug, Clone)]
pub struct PriorityQueue<T: Copy, F: Fn(&T, &T) -> bool> {
    c: Vec<T>,
    less: F,
}

impl<T: Copy, F: Fn(&T, &T) -> bool> PriorityQueue<T, F> {
    pub fn new(less: F) -> Self {
        PriorityQueue { c: Vec::new(), less }
    }

    /// Empties the queue, keeping its allocation.
    pub fn clear(&mut self) {
        self.c.clear();
    }

    pub fn is_empty(&self) -> bool {
        self.c.is_empty()
    }

    pub fn len(&self) -> usize {
        self.c.len()
    }

    /// `push`: `c.push_back(x); std::push_heap(c.begin(), c.end(), comp)`.
    #[inline]
    pub fn push(&mut self, x: T) {
        self.c.push(x);
        let hole = self.c.len() - 1;
        push_heap_hole(&mut self.c, hole, 0, x, &self.less);
    }

    /// `top()` (panics when empty, where the C++ is undefined).
    #[inline]
    pub fn top(&self) -> T {
        self.c[0]
    }

    /// `pop()`: `std::pop_heap(c.begin(), c.end(), comp); c.pop_back()`.
    #[inline]
    pub fn pop(&mut self) {
        let len = self.c.len();
        if len > 1 {
            // __pop_heap(first, last - 1, last - 1): the last element is re-inserted from the top
            let last = len - 1;
            let value = self.c[last];
            self.c[last] = self.c[0];
            adjust_heap(&mut self.c, 0, last, value, &self.less);
        }
        self.c.pop();
    }
}

/// `std::partial_sort(c.begin(), c.begin() + middle, c.end(), comp)`: `__heap_select`, then
/// `__sort_heap`. Afterwards `c[..middle]` holds the `middle` smallest elements in order.
pub fn partial_sort<T: Copy>(c: &mut [T], middle: usize, less: &impl Fn(&T, &T) -> bool) {
    debug_assert!(middle <= c.len());
    // __heap_select(first, middle, last)
    make_heap(c, middle, less);
    for i in middle..c.len() {
        if less(&c[i], &c[0]) {
            // __pop_heap(first, middle, i)
            let value = c[i];
            c[i] = c[0];
            adjust_heap(c, 0, middle, value, less);
        }
    }
    // __sort_heap(first, middle)
    let mut last = middle;
    while last > 1 {
        last -= 1;
        let value = c[last];
        c[last] = c[0];
        adjust_heap(c, 0, last, value, less);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            self.0 >> 33
        }
    }

    #[test]
    fn queue_pops_in_pair_order() {
        let mut g = Lcg(3);
        let mut q = PriorityQueue::new(pair_less);
        let mut all = Vec::new();
        for _ in 0..500 {
            let x = ((g.next() % 7) as f32 - 3.0, (g.next() % 50) as i32);
            q.push(x);
            all.push(x);
        }
        // max first; equal priorities by the larger id
        all.sort_by(|a, b| if pair_less(a, b) { std::cmp::Ordering::Greater } else if pair_less(b, a) { std::cmp::Ordering::Less } else { std::cmp::Ordering::Equal });
        let mut got = Vec::new();
        while !q.is_empty() {
            got.push(q.top());
            q.pop();
        }
        assert_eq!(got, all);
    }

    #[test]
    fn signed_zeros_compare_equal() {
        assert!(!pair_less(&(-0.0, 1), &(0.0, 1)));
        assert!(!pair_less(&(0.0, 1), &(-0.0, 1)));
        assert!(pair_less(&(0.0, 1), &(-0.0, 2)));
    }

    #[test]
    fn partial_sort_selects_the_smallest_in_order() {
        let mut g = Lcg(11);
        for (n, m) in [(1usize, 1usize), (2, 1), (10, 3), (100, 20), (100, 100), (1000, 20), (37, 0)] {
            let v: Vec<(f32, i32)> = (0..n).map(|i| ((g.next() % 13) as f32, i as i32)).collect();
            let mut a = v.clone();
            partial_sort(&mut a, m, &pair_less);
            let mut b = v.clone();
            b.sort_by(|x, y| x.partial_cmp(y).unwrap());
            assert_eq!(&a[..m], &b[..m], "n={n} m={m}");
            // the rest is a permutation of the remaining elements
            let mut rest: Vec<_> = a[m..].to_vec();
            rest.sort_by(|x, y| x.partial_cmp(y).unwrap());
            assert_eq!(rest, b[m..].to_vec());
        }
    }

    /// The libstdc++ moves themselves, on a case with NaN priorities, where the result depends
    /// on the algorithm: worked out by hand from stl_heap.h.
    #[test]
    fn heap_moves_follow_libstdcxx() {
        // push 3 then 1 then 2 (as (p, id)): heap [3, 1, 2]; pop -> [2, 1]
        let mut q = PriorityQueue::new(pair_less);
        for p in [3.0, 1.0, 2.0] {
            q.push((p, 0));
        }
        assert_eq!(q.c, vec![(3.0, 0), (1.0, 0), (2.0, 0)]);
        q.pop();
        assert_eq!(q.c, vec![(2.0, 0), (1.0, 0)]);
        // a NaN priority compares by id only: pushed after (5, 1), (NaN, 0) stays a child
        let mut q = PriorityQueue::new(pair_less);
        q.push((5.0, 1));
        q.push((f32::NAN, 0));
        assert_eq!(q.top(), (5.0, 1));
        let mut q = PriorityQueue::new(pair_less);
        q.push((5.0, 0));
        q.push((f32::NAN, 1));
        assert!(q.top().0.is_nan());
    }
}
