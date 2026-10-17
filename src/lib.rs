//! Codon optimization engine.
//!
//! The problem splits in two. *Composition* decides which codons are used and how many
//! of each, which fixes the deviation from the host profile and the codon adaptation
//! index. *Arrangement* decides where each one goes, which fixes every sequence
//! constraint. Moves that swap two positions holding the same amino acid cannot change
//! the multiset, so the search repairs constraints without ever disturbing the profile.

pub mod arrange;
pub mod code;
pub mod composition;
pub mod fasta;
pub mod host;
pub mod metrics;
pub mod small;
pub mod rng;
pub mod testdata;
