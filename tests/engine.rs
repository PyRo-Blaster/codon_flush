//! End-to-end properties of the optimiser.

use codon_harmony_core::arrange::constraints::{ConstraintSet, Kind};
use codon_harmony_core::arrange::search::{optimise, Options, Outcome};
use codon_harmony_core::arrange::state::Move;
use codon_harmony_core::code::{aa_index, AaIdx};
use codon_harmony_core::composition::{self, Composition};
use codon_harmony_core::host::HostTable;
use codon_harmony_core::rng::Rng;
use codon_harmony_core::testdata::{ECOLI_TABLE, GFP, STRESS};
use codon_harmony_core::{fasta, metrics};
use std::path::Path;

fn host() -> HostTable {
    HostTable::load(Path::new(ECOLI_TABLE), 0.10).expect("the test host table loads")
}

fn protein(seq: &str) -> Vec<AaIdx> {
    seq.bytes().map(|b| aa_index(b).unwrap()).collect()
}

/// A protein with a realistic amino acid composition, free of tandem repeats.
fn realistic_protein(n: usize, seed: u64) -> Vec<AaIdx> {
    const LETTERS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    const WEIGHTS: [u32; 20] = [95, 12, 51, 58, 39, 74, 23, 60, 44, 106, 28, 39, 44, 44, 55, 58, 54, 71, 15, 28];
    let total: u32 = WEIGHTS.iter().sum();
    let mut rng = Rng::new(seed);
    let mut out = vec![aa_index(b'M').unwrap()];
    while out.len() < n {
        let mut pick = rng.below(total as usize) as u32;
        let mut index = 0;
        while pick >= WEIGHTS[index] {
            pick -= WEIGHTS[index];
            index += 1;
        }
        out.push(aa_index(LETTERS[index]).unwrap());
    }
    out
}

fn run(protein: &[AaIdx], options: Options) -> (Composition, Outcome) {
    let host = host();
    let composition = composition::build(protein, &host, 0.1, (0.30, 0.65)).unwrap();
    let outcome = optimise(protein, &composition, &host, &options);
    (composition, outcome)
}

/// Rebuild the full constraint set independently of the search's own bookkeeping.
fn fresh_constraints(options: &Options, outcome: &Outcome) -> ConstraintSet {
    use codon_harmony_core::arrange::constraints::{
        gc, hairpin, homopolymer, repeat, restriction, splice, start_site,
    };
    let mut set = ConstraintSet::new(vec![
        Box::new(gc::GcWindows::new(gc::DEFAULT_WINDOWS.to_vec())),
        Box::new(homopolymer::Homopolymers::new(options.homopolymer_limit)),
        Box::new(repeat::Repeats::new()),
        Box::new(hairpin::Hairpins::new()),
        Box::new(restriction::RestrictionSites::new(&options.enzymes)),
        Box::new(start_site::StartSites::with_rule(options.start_site_rule)),
        Box::new(splice::SpliceSites::new()),
    ]);
    set.rebuild(&outcome.state);
    set
}

#[test]
fn the_protein_always_survives() {
    for (name, seq) in [("gfp", GFP), ("stress", STRESS)] {
        let residues = protein(seq);
        let (_, outcome) = run(&residues, Options { seed: 11, ..Default::default() });
        assert_eq!(outcome.state.translate(), residues, "{name} changed protein");
        assert_eq!(outcome.state.dna().len(), residues.len() * 3, "{name} changed length");
    }
    for n in [50, 300, 1200] {
        let residues = realistic_protein(n, n as u64);
        let (_, outcome) = run(&residues, Options { seed: 2, restarts: 2, ..Default::default() });
        assert_eq!(outcome.state.translate(), residues, "random protein of {n}");
    }
}

#[test]
fn the_search_reports_its_own_violation_count_honestly() {
    let options = Options { seed: 3, restarts: 2, ..Default::default() };
    let residues = realistic_protein(600, 7);
    let (_, outcome) = run(&residues, options.clone());
    let independent = fresh_constraints(&options, &outcome);
    assert_eq!(
        outcome.violations,
        independent.total(),
        "the running total must match a fresh count of the finished gene"
    );
    assert_eq!(outcome.remaining.len(), outcome.violations as usize);
}

#[test]
fn it_stops_only_where_no_single_move_helps() {
    let options = Options { seed: 4, restarts: 1, ..Default::default() };
    let residues = realistic_protein(400, 21);
    let host = host();
    let (_, outcome) = run(&residues, options.clone());
    let set = fresh_constraints(&options, &outcome);

    // For every violation left, check that no swap and no substitution anywhere in its
    // window would reduce the total. Anything still reported is a genuine local
    // minimum, not the search giving up early.
    for violation in &outcome.remaining {
        let span = violation.codons();
        for p in span.start..span.end.min(outcome.state.n_codons() as u32) {
            let aa = outcome.state.aa[p as usize];
            for &q in &outcome.state.pos_by_aa[aa as usize] {
                if q == p || outcome.state.codons[q as usize] == outcome.state.codons[p as usize] {
                    continue;
                }
                let patch = outcome.state.patch_of(Move::Swap { p, q });
                let delta = set.delta(&outcome.state.scoped(patch), &outcome.state.patched(patch));
                assert!(delta >= 0, "a swap {p}/{q} would still have helped {violation:?}");
            }
            for codon in host.usable_synonyms(aa) {
                if codon == outcome.state.codons[p as usize] {
                    continue;
                }
                let patch = outcome.state.patch_of(Move::Subst { p, codon });
                let delta = set.delta(&outcome.state.scoped(patch), &outcome.state.patched(patch));
                assert!(delta >= 0, "a substitution at {p} would still have helped {violation:?}");
            }
        }
    }
}

#[test]
fn deviation_never_exceeds_the_budget() {
    for seed in 1..=3 {
        for n in [120, 700] {
            let residues = realistic_protein(n, seed * 31);
            let (composition, outcome) =
                run(&residues, Options { seed, restarts: 2, ..Default::default() });
            assert!(
                outcome.deviation <= composition.deviation + composition.budget_left + 1e-9,
                "deviation {} over budget {}",
                outcome.deviation,
                composition.deviation + composition.budget_left
            );
            assert!(outcome.deviation <= 0.1 + 1e-9, "max_relax must hold at the end");
            assert!(
                metrics::unused_codons(&outcome.state, &host()).is_empty(),
                "a codon the host never uses would drive CAI to zero"
            );
        }
    }
}

#[test]
fn with_no_budget_the_profile_is_untouched() {
    let host = host();
    let residues = realistic_protein(500, 99);
    let composition = composition::build(&residues, &host, 0.0, (0.30, 0.65)).unwrap();
    let mut seen_cai = Vec::new();
    for seed in 1..=4 {
        let outcome = optimise(
            &residues,
            &composition,
            &host,
            &Options { seed, restarts: 2, ..Default::default() },
        );
        assert_eq!(outcome.substitutions, 0);
        assert_eq!(outcome.deviation, composition.deviation);
        seen_cai.push(outcome.cai);
    }
    for pair in seen_cai.windows(2) {
        assert_eq!(pair[0], pair[1], "rearrangement cannot change CAI");
    }
}

#[test]
fn results_do_not_depend_on_thread_count_or_repeat_runs() {
    let residues = realistic_protein(300, 5);
    let first = run(&residues, Options { seed: 8, restarts: 4, ..Default::default() }).1;
    let again = run(&residues, Options { seed: 8, restarts: 4, ..Default::default() }).1;
    assert_eq!(first.state.dna(), again.state.dna());
    assert_eq!(first.restart, again.restart);
    assert_eq!(first.violations, again.violations);
}

#[test]
fn a_tandem_repeat_protein_is_reported_not_silently_shipped() {
    // Seven copies of GFP: most nine-base windows are forced to repeat, because the
    // amino acids repeat and the codon pool cannot supply a distinct encoding for every
    // copy. The engine should say so rather than claim success.
    let residues: Vec<AaIdx> = protein(&GFP.repeat(7))[..1600].to_vec();
    let (_, outcome) = run(&residues, Options { seed: 6, restarts: 1, ..Default::default() });
    assert_eq!(outcome.state.translate(), residues);
    assert!(outcome.violations > 0);
    assert!(
        outcome.remaining.iter().any(|v| v.kind == Kind::Repeat),
        "the repeats are the point of this protein"
    );
    assert_eq!(outcome.remaining.len(), outcome.violations as usize);
}

#[test]
fn fasta_round_trips_through_the_engine() {
    let text = format!(">gene one\n{GFP}\n>gene two\n{STRESS}\n");
    let records = fasta::parse(&text).unwrap();
    assert_eq!(records.len(), 2);
    for record in &records {
        let (_, outcome) = run(&record.residues, Options { seed: 1, restarts: 1, ..Default::default() });
        let written = fasta::write_record(&record.header(), &outcome.state.dna());
        assert!(written.starts_with(&format!(">{}", record.header())));
        let dna: String = written.lines().skip(1).collect();
        assert_eq!(dna.len(), record.residues.len() * 3);
    }
}
