//! Restriction enzyme recognition sites.

use std::ops::Range;

use crate::arrange::state::{SeqState, View};
use crate::code::{base_index, Base};

use super::{window_delta, Constraint, Kind, Violation};

/// The enzymes the tool has always offered, with their recognition sequences.
pub const DEFAULT_ENZYMES: [(&str, &str); 7] = [
    ("NdeI", "CATATG"),
    ("XhoI", "CTCGAG"),
    ("HpaI", "GTTAAC"),
    ("PstI", "CTGCAG"),
    ("EcoRV", "GATATC"),
    ("NcoI", "CCATGG"),
    ("BamHI", "GGATCC"),
];

#[derive(Debug)]
pub struct RestrictionSites {
    sites: Vec<(String, Vec<Base>)>,
    longest: usize,
    total: u32,
}

impl RestrictionSites {
    pub fn new(sites: &[(String, String)]) -> Self {
        let sites: Vec<(String, Vec<Base>)> = sites
            .iter()
            .map(|(name, seq)| {
                let bases = seq.bytes().filter_map(base_index).collect::<Vec<Base>>();
                (name.clone(), bases)
            })
            .filter(|(_, bases)| !bases.is_empty())
            .collect();
        let longest = sites.iter().map(|(_, b)| b.len()).max().unwrap_or(0);
        Self { sites, longest, total: 0 }
    }

    pub fn from_names(names: &[String]) -> Result<Self, String> {
        let mut chosen = Vec::new();
        for name in names {
            match DEFAULT_ENZYMES.iter().find(|(n, _)| n.eq_ignore_ascii_case(name)) {
                Some((n, seq)) => chosen.push((n.to_string(), seq.to_string())),
                None => {
                    return Err(format!(
                        "unknown restriction enzyme {name:?}; known ones are {}",
                        DEFAULT_ENZYMES.iter().map(|(n, _)| *n).collect::<Vec<_>>().join(", ")
                    ))
                }
            }
        }
        Ok(Self::new(&chosen))
    }

    fn count_in(&self, view: &View<'_>, range: Range<usize>) -> u32 {
        let n = view.n_bases();
        let mut total = 0;
        for start in range.start..range.end.min(n) {
            for (_, site) in &self.sites {
                if start + site.len() <= n
                    && (0..site.len()).all(|k| view.base(start + k) == site[k])
                {
                    total += 1;
                }
            }
        }
        total
    }
}

impl Constraint for RestrictionSites {
    fn kind(&self) -> Kind {
        Kind::Restriction
    }

    fn rebuild(&mut self, state: &SeqState) {
        self.total = self.count_in(&state.view(), 0..state.n_bases());
    }

    fn total(&self) -> u32 {
        self.total
    }

    fn delta(&self, before: &View<'_>, after: &View<'_>) -> i32 {
        window_delta(before, after, self.longest, |v, r| self.count_in(v, r))
    }

    fn apply(&mut self, before: &View<'_>, after: &View<'_>) {
        self.total = (self.total as i32 + self.delta(before, after)) as u32;
    }

    fn report(&self, state: &SeqState) -> Vec<Violation> {
        let view = state.view();
        let n = state.n_bases();
        let mut out = Vec::new();
        for start in 0..n {
            for (name, site) in &self.sites {
                if start + site.len() <= n
                    && (0..site.len()).all(|k| view.base(start + k) == site[k])
                {
                    out.push(Violation {
                        kind: Kind::Restriction,
                        at: start as u32,
                        len: site.len() as u32,
                        detail: format!("{name} site at base {}", start + 1),
                    });
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrange::constraints::harness;

    fn make() -> RestrictionSites {
        RestrictionSites::new(
            &DEFAULT_ENZYMES
                .iter()
                .map(|(n, s)| (n.to_string(), s.to_string()))
                .collect::<Vec<_>>(),
        )
    }

    fn oracle(state: &SeqState) -> u32 {
        let dna = state.dna();
        let mut total = 0;
        for (_, site) in DEFAULT_ENZYMES.iter() {
            let mut from = 0;
            while let Some(at) = dna[from..].find(site) {
                total += 1;
                from += at + 1;
            }
        }
        total
    }

    #[test]
    fn matches_the_oracle() {
        harness::assert_matches_oracle(make, oracle, 200);
    }

    #[test]
    fn incremental_matches_rebuild() {
        harness::assert_incremental_matches_rebuild(make, 20_000);
    }

    #[test]
    fn unknown_enzyme_is_rejected() {
        assert!(RestrictionSites::from_names(&["NdeI".into()]).is_ok());
        let err = RestrictionSites::from_names(&["Nonsense".into()]).unwrap_err();
        assert!(err.contains("Nonsense") && err.contains("NdeI"));
    }
}
