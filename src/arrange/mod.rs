//! Arrangement: where each codon goes.
//!
//! The composition is fixed before this module runs, so every move here either keeps
//! the codon multiset exactly (a swap) or spends an explicit slice of the deviation
//! budget (a substitution).

pub mod constraints;
pub mod search;
pub mod state;
