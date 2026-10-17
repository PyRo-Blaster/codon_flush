//! Small deterministic RNG (SplitMix64), so a seed reproduces a run exactly.

#[derive(Debug, Clone)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed.wrapping_add(0x9E37_79B9_7F4A_7C15))
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `0..n`. Panics if `n == 0`.
    pub fn below(&mut self, n: usize) -> usize {
        assert!(n > 0, "below(0) has no value to return");
        (self.next_u64() % n as u64) as usize
    }

    /// Uniform in `[lo, hi]`.
    pub fn between(&mut self, lo: usize, hi: usize) -> usize {
        lo + self.below(hi - lo + 1)
    }

    pub fn chance(&mut self, p: f64) -> bool {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64 > 1.0 - p
    }

    pub fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            items.swap(i, self.below(i + 1));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_gives_the_same_stream() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        let mut c = Rng::new(43);
        let xs: Vec<u64> = (0..16).map(|_| a.next_u64()).collect();
        let ys: Vec<u64> = (0..16).map(|_| b.next_u64()).collect();
        let zs: Vec<u64> = (0..16).map(|_| c.next_u64()).collect();
        assert_eq!(xs, ys);
        assert_ne!(xs, zs);
    }

    #[test]
    fn below_stays_in_range_and_covers_it() {
        let mut rng = Rng::new(7);
        let mut seen = [false; 5];
        for _ in 0..500 {
            let v = rng.below(5);
            assert!(v < 5);
            seen[v] = true;
        }
        assert!(seen.iter().all(|&s| s), "every value should appear");
    }

    #[test]
    fn chance_is_roughly_calibrated() {
        let mut rng = Rng::new(11);
        let hits = (0..10_000).filter(|_| rng.chance(0.25)).count();
        assert!((2000..3000).contains(&hits), "got {hits} hits out of 10000");
        let mut rng = Rng::new(11);
        assert_eq!((0..100).filter(|_| rng.chance(0.0)).count(), 0);
        let mut rng = Rng::new(11);
        assert_eq!((0..100).filter(|_| rng.chance(1.0)).count(), 100);
    }

    #[test]
    fn shuffle_keeps_every_item() {
        let mut rng = Rng::new(3);
        let mut items: Vec<u32> = (0..50).collect();
        rng.shuffle(&mut items);
        assert_ne!(items, (0..50).collect::<Vec<u32>>());
        items.sort_unstable();
        assert_eq!(items, (0..50).collect::<Vec<u32>>());
    }
}
