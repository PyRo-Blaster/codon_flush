use codon_harmony_core::{arrange, code, composition, host::HostTable};
use arrange::constraints::{Constraint, Kind, Violation, repeat::Repeats, hairpin::Hairpins};
use arrange::state::{SeqState, Move};
use std::{hint::black_box, time::Instant, path::Path};
#[path = "probe_support/fast_repeat.rs"] mod fast_repeat;
fn state(dna: &str) -> SeqState {
    let codons: Vec<_> = dna.as_bytes().chunks_exact(3).map(|c| code::parse_codon(std::str::from_utf8(c).unwrap()).unwrap()).collect();
    SeqState::new(codons.iter().map(|&c| code::codon_aa(c)).collect(), codons)
}
fn hairpin_oracle(s: &SeqState) -> u32 {
    (0..s.n_bases().saturating_sub(9)).step_by(10).map(|p| {
        (0..s.n_bases().saturating_sub(9)).filter(|&q| (0..10).all(|j| {
            let partner = code::complement(s.nt[p + 9 - j]);
            s.nt[q+j] == partner || partner == 1 && s.nt[q+j] == 3
        })).count() as u32
    }).sum()
}
fn micros() {
    println!("repeat_delta,aa,baseline_ns,prototype_ns,speedup");
    for n in [100, 400, 1600, 6400] {
        let s = state(&"GCT".repeat(n));
        let patch = s.patch_of(Move::Subst {p: (n/2) as u32, codon: code::parse_codon("GCC").unwrap()});
        let before = s.scoped(patch); let after = s.patched(patch);
        let mut a = Repeats::new(); a.rebuild(&s);
        let mut b = fast_repeat::Repeats::new(); b.rebuild(&s);
        assert_eq!(a.delta(&before, &after), b.delta(&before, &after));
        let repetitions = 10000;
        let start = Instant::now(); for _ in 0..repetitions { black_box(a.delta(black_box(&before), black_box(&after))); }
        let old = start.elapsed().as_nanos() as f64 / repetitions as f64;
        let start = Instant::now(); for _ in 0..repetitions { black_box(b.delta(black_box(&before), black_box(&after))); }
        let new = start.elapsed().as_nanos() as f64 / repetitions as f64;
        println!("repeat_delta,{n},{old:.0},{new:.0},{:.2}", old/new);
    }
    // Dynamic comparison, with repeated and varied codons, including overlapping swaps.
    let mut rng = codon_harmony_core::rng::Rng::new(1234);
    let mut s = state(&"GCTGCCGCAGCG".repeat(100));
    let mut a = Repeats::new(); a.rebuild(&s);
    let mut b = fast_repeat::Repeats::new(); b.rebuild(&s);
    for _ in 0..20000 {
        let p = rng.below(s.n_codons()) as u32;
        let mv = if rng.below(2) == 0 { Move::Swap {p, q: rng.below(s.n_codons()) as u32} }
        else { Move::Subst {p, codon: code::parse_codon(["GCT","GCC","GCA","GCG"][rng.below(4)]).unwrap()} };
        let patch = s.patch_of(mv);
        let before=s.scoped(patch); let after=s.patched(patch);
        assert_eq!(a.delta(&before,&after),b.delta(&before,&after));
        a.apply(&before,&after); b.apply(&before,&after); s.apply(mv);
        assert_eq!(a.total(),b.total());
    }
    println!("differential_repeat_moves,20000,passed");
    let s = state("GGGGGGGAAATTTCTCCCCCAC");
    let mut h = Hairpins::new(); h.rebuild(&s);
    println!("hairpin_wobble,reported={},oracle={}",h.total(),hairpin_oracle(&s));
    for n in [100,200,400] {
        let s = state(&"AT".repeat(3*n/2));
        let mut h = Hairpins::new(); h.rebuild(&s);
        let start = Instant::now(); let report = h.report(&s);
        println!("hairpin_report,aa={n},entries={},micros={}",report.len(),start.elapsed().as_micros());
        let mut r = Repeats::new(); r.rebuild(&s);
        let start = Instant::now(); let report = r.report(&s);
        println!("repeat_report,aa={n},detail_bytes={},micros={}",report.iter().map(|v|v.detail.len()).sum::<usize>(),start.elapsed().as_micros());
    }
    let host = HostTable::load(Path::new("tests/data/e_coli_316407.json"),0.1).unwrap();
    let aa: Vec<_> = codon_harmony_core::testdata::GFP.bytes().map(|b| code::aa_index(b).unwrap()).collect();
    let c = composition::build(&aa,&host,0.001,(0.3,0.65)).unwrap();
    println!("max_relax,requested=0.001,composition_deviation={},budget_left={}",c.deviation,c.budget_left);
}
fn engine() {
    use arrange::search::{optimise,Options};
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global().unwrap();
    let host = HostTable::load(Path::new("tests/data/e_coli_316407.json"),0.1).unwrap();
    println!("case,seed,millis,violations,deviation,moves,dna_hash");
    for (name,n) in [("GFP",238),("random1600",1600),("random6400",6400)] {
        let mut rng = codon_harmony_core::rng::Rng::new(2026);
        let aa: Vec<_> = if name == "GFP" { codon_harmony_core::testdata::GFP.bytes().map(|b|code::aa_index(b).unwrap()).collect() }
        else { (0..n).map(|_| rng.below(20) as u8).collect() };
        let c = composition::build(&aa,&host,0.01,(0.3,0.65)).unwrap();
        for seed in [1,2,3] {
            let options = Options {seed, restarts:1,..Default::default()};
            let start=Instant::now(); let o=optimise(&aa,&c,&host,&options); let elapsed=start.elapsed().as_secs_f64()*1000.;
            let hash=o.state.nt.iter().fold(0u64, |h,&b|h.wrapping_mul(131).wrapping_add(b as u64));
            println!("{name},{seed},{elapsed:.3},{},{:.8},{},{hash}",o.violations,o.deviation,o.moves);
        }
    }
}
fn main() { if std::env::args().nth(1).as_deref()==Some("engine") {engine()} else {micros()} }
