//! A fixed-capacity vector, for the hot paths where a `Vec` would mean an allocation
//! per candidate move. A move touches at most two codons, so every such buffer has a
//! small, known upper bound.

#[derive(Debug, Clone, Copy)]
pub struct Small<T: Copy + Default, const N: usize> {
    items: [T; N],
    len: usize,
}

impl<T: Copy + Default, const N: usize> Default for Small<T, N> {
    fn default() -> Self {
        Self { items: [T::default(); N], len: 0 }
    }
}

impl<T: Copy + Default, const N: usize> Small<T, N> {
    pub fn new() -> Self {
        Self::default()
    }

    /// Pushes if there is room. The capacities here are derived from how far a move can
    /// reach, so a full buffer means an assumption changed.
    pub fn push(&mut self, item: T) {
        debug_assert!(self.len < N, "fixed buffer of {N} is full");
        if self.len < N {
            self.items[self.len] = item;
            self.len += 1;
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn as_slice(&self) -> &[T] {
        &self.items[..self.len]
    }

    pub fn last_mut(&mut self) -> Option<&mut T> {
        if self.len == 0 {
            None
        } else {
            Some(&mut self.items[self.len - 1])
        }
    }

    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.as_slice().iter()
    }
}

impl<'a, T: Copy + Default, const N: usize> IntoIterator for &'a Small<T, N> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn holds_items_in_order() {
        let mut small: Small<(usize, usize), 4> = Small::new();
        assert!(small.is_empty());
        small.push((1, 2));
        small.push((3, 4));
        assert_eq!(small.len(), 2);
        assert_eq!(small.as_slice(), &[(1, 2), (3, 4)]);
        small.last_mut().unwrap().1 = 9;
        assert_eq!(small.as_slice(), &[(1, 2), (3, 9)]);
        assert_eq!(small.iter().count(), 2);
    }
}
