//! End-to-end runs of the binary: the candidate FASTA, the JSON report, the exit
//! status, and the settings the CLI refuses.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use codon_harmony_core::testdata::{ECOLI_TABLE, GFP};

/// A directory of this test's own, under the target directory so `cargo clean`
/// removes it.
fn workspace(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn write(path: &Path, text: &str) {
    std::fs::write(path, text).unwrap();
}

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_codon_harmony"))
        .args(args)
        .output()
        .expect("the binary should run")
}

/// Every record in a FASTA, as (header, sequence).
fn records(text: &str) -> Vec<(String, String)> {
    let mut out: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        match line.strip_prefix('>') {
            Some(header) => out.push((header.to_string(), String::new())),
            None => out.last_mut().expect("a header first").1.push_str(line.trim()),
        }
    }
    out
}

fn report_of(path: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn translate(dna: &str) -> String {
    use codon_harmony_core::code::{aa_letter, codon_aa, parse_codon};
    dna.as_bytes()
        .chunks(3)
        .map(|c| {
            let codon = parse_codon(std::str::from_utf8(c).unwrap()).unwrap();
            aa_letter(codon_aa(codon)) as char
        })
        .collect()
}

#[test]
fn five_ranked_candidates_per_protein_with_a_matching_report() {
    let dir = workspace("five");
    let input = dir.join("proteins.fasta");
    let output = dir.join("genes.fasta");
    let report_path = dir.join("report.json");
    write(&input, &format!(">gfp green fluorescent protein\n{GFP}\n>half\n{}\n", &GFP[..90]));

    let out = run(&[
        "--input", input.to_str().unwrap(),
        "--output", output.to_str().unwrap(),
        "--host-table", ECOLI_TABLE,
        "--report-json", report_path.to_str().unwrap(),
        "--quiet",
    ]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));

    let written = records(&std::fs::read_to_string(&output).unwrap());
    assert_eq!(written.len(), 10, "five candidates for each of two proteins");

    let report = report_of(&report_path);
    assert_eq!(report["candidates_requested"], 5);
    assert_eq!(report["cai_threshold"], 0.8);
    assert_eq!(report["summary"]["exit_code"], 0);
    assert!(report["cai_convention"].as_str().unwrap().contains("Met and Trp"));

    let mut seen_ids = std::collections::HashSet::new();
    for (index, record) in report["records"].as_array().unwrap().iter().enumerate() {
        let candidates = record["candidates"].as_array().unwrap();
        assert_eq!(candidates.len(), 5, "{}", record["id"]);
        assert_eq!(record["candidate_count"], 5);
        let mut last = f64::INFINITY;
        let mut genes = std::collections::HashSet::new();
        for (rank, candidate) in candidates.iter().enumerate() {
            let cai = candidate["cai"].as_f64().unwrap();
            assert!(cai <= last + 1e-12, "candidates must descend by CAI");
            last = cai;
            assert_eq!(candidate["rank"], rank as u64 + 1);
            let dna = candidate["dna"].as_str().unwrap();
            assert!(genes.insert(dna.to_string()), "candidates must be distinct");
            assert!(seen_ids.insert(candidate["output_id"].as_str().unwrap().to_string()));

            // the gene in the report is the gene in the FASTA, and it translates back
            let (header, sequence) = &written[index * 5 + rank];
            assert!(header.starts_with(candidate["output_id"].as_str().unwrap()));
            assert!(header.contains(&format!("source={}", record["id"].as_str().unwrap())));
            assert!(header.contains(&format!("rank={}/5", rank + 1)));
            assert_eq!(sequence, dna);
            assert_eq!(translate(dna).len(), record["residues"].as_u64().unwrap() as usize);
        }
        assert_eq!(record["status"], "ok");
        assert_eq!(record["low_cai"], false);
        assert!(candidates[0]["cai"].as_f64().unwrap() > 0.8);
    }
}

#[test]
fn the_retry_happens_once_and_is_recorded() {
    let dir = workspace("retry");
    let input = dir.join("p.fasta");
    let output = dir.join("g.fasta");
    let report_path = dir.join("r.json");
    write(&input, &format!(">gfp\n{GFP}\n"));

    // A threshold nothing can reach forces the second run and the low-CAI fallback.
    let out = run(&[
        "--input", input.to_str().unwrap(),
        "--output", output.to_str().unwrap(),
        "--host-table", ECOLI_TABLE,
        "--report-json", report_path.to_str().unwrap(),
        "--cai-threshold", "0.99",
        "--quiet",
    ]);
    assert!(out.status.success(), "a complete low-CAI set is not a failure");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("neither run produced CAI > 0.99"), "{stderr}");

    let report = report_of(&report_path);
    let record = &report["records"][0];
    assert_eq!(record["status"], "low_cai");
    assert_eq!(record["low_cai"], true);
    assert_eq!(record["candidate_count"], 5);
    let attempts = record["attempts"].as_array().unwrap();
    assert_eq!(attempts.len(), 2, "at most two runs per protein");
    assert_eq!(attempts[0]["max_relax"], 0.1);
    assert_eq!(attempts[1]["max_relax"], 0.2);
    assert_eq!(attempts[1]["restarts"].as_u64().unwrap(), 2 * attempts[0]["restarts"].as_u64().unwrap());
    assert_eq!(
        attempts[1]["expansion_budget"].as_u64().unwrap(),
        2 * attempts[0]["expansion_budget"].as_u64().unwrap()
    );
    assert!(record["retry_reason"].as_str().unwrap().contains("threshold"));
    assert_eq!(report["summary"]["low_cai"], 1);
    assert_eq!(report["summary"]["exit_code"], 0);
}

#[test]
fn a_protein_with_one_encoding_reports_the_shortage_and_exits_nonzero() {
    let dir = workspace("shortage");
    let input = dir.join("p.fasta");
    let output = dir.join("g.fasta");
    let report_path = dir.join("r.json");
    write(&input, ">mw only methionine and tryptophan\nMWMWM\n");

    let out = run(&[
        "--input", input.to_str().unwrap(),
        "--output", output.to_str().unwrap(),
        "--host-table", ECOLI_TABLE,
        "--report-json", report_path.to_str().unwrap(),
        "--quiet",
    ]);
    assert_eq!(out.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("insufficient_candidates"), "{stderr}");

    let written = records(&std::fs::read_to_string(&output).unwrap());
    assert_eq!(written.len(), 1, "the one sequence that exists, never padded");
    assert_eq!(written[0].1, "ATGTGGATGTGGATG");
    let report = report_of(&report_path);
    assert_eq!(report["records"][0]["status"], "insufficient_candidates");
    assert_eq!(report["records"][0]["candidate_count"], 1);
    assert_eq!(report["summary"]["exit_code"], 1);
}

#[test]
fn a_failing_protein_does_not_stop_its_peers() {
    let dir = workspace("mixed");
    let input = dir.join("p.fasta");
    let output = dir.join("g.fasta");
    let report_path = dir.join("r.json");
    write(&input, &format!(">mw\nMWM\n>gfp\n{GFP}\n"));

    let out = run(&[
        "--input", input.to_str().unwrap(),
        "--output", output.to_str().unwrap(),
        "--host-table", ECOLI_TABLE,
        "--report-json", report_path.to_str().unwrap(),
        "--quiet",
    ]);
    assert_eq!(out.status.code(), Some(1), "the batch fails because one record did");
    let written = records(&std::fs::read_to_string(&output).unwrap());
    assert_eq!(written.len(), 6, "one for MWM and five for GFP");
    assert!(written[0].0.starts_with("record_000001_candidate_01"));
    assert!(written[1].0.starts_with("record_000002_candidate_01"));
    let report = report_of(&report_path);
    assert_eq!(report["records"][1]["status"], "ok", "the peer still succeeded");
}

#[test]
fn the_same_seed_gives_the_same_output_whatever_the_thread_count() {
    let dir = workspace("threads");
    let input = dir.join("p.fasta");
    write(&input, &format!(">gfp\n{GFP}\n>half\n{}\n", &GFP[..90]));
    let mut outputs = Vec::new();
    for threads in ["1", "4"] {
        let output = dir.join(format!("g{threads}.fasta"));
        let out = run(&[
            "--input", input.to_str().unwrap(),
            "--output", output.to_str().unwrap(),
            "--host-table", ECOLI_TABLE,
            "--threads", threads,
            "--quiet",
        ]);
        assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
        outputs.push(std::fs::read_to_string(&output).unwrap());
    }
    assert_eq!(outputs[0], outputs[1]);
}

#[test]
fn impossible_settings_are_refused_before_anything_is_optimised() {
    let dir = workspace("settings");
    let input = dir.join("p.fasta");
    let output = dir.join("g.fasta");
    write(&input, &format!(">gfp\n{}\n", &GFP[..30]));
    let base = [
        "--input", input.to_str().unwrap(),
        "--output", output.to_str().unwrap(),
        "--host-table", ECOLI_TABLE,
    ];
    let refuse = |extra: &[&str], expect: &str| {
        let mut args = base.to_vec();
        args.extend_from_slice(extra);
        let out = run(&args);
        assert!(!out.status.success(), "{extra:?} should be refused");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(stderr.contains(expect), "{extra:?} said {stderr}");
    };
    refuse(&["--cai-threshold", "1.0"], "below 1");
    refuse(&["--max-relax", "1.5"], "between 0 and 1");
    refuse(&["--candidates", "0"], "at least 1");
    refuse(&["--retry-max-relax", "0.05"], "at least --max-relax");
    refuse(&["--gc-low", "0.7"], "--gc-low must be below --gc-high");

    // and a report written over the input would destroy it
    let mut args = base.to_vec();
    args.extend_from_slice(&["--report-json", input.to_str().unwrap()]);
    let out = run(&args);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("must be different files"));
    assert!(std::fs::read_to_string(&input).unwrap().starts_with('>'), "input untouched");
}

// --- regressions for the review on 29bfb4d -------------------------------------

#[test]
fn a_positive_deviation_limit_below_the_rounding_floor_is_refused() {
    // Was: five candidates at deviation 0.0083 under --max-relax 0.001, exit 0. The
    // limit the caller set is a bound, not a suggestion; the rounding floor cannot
    // quietly replace it.
    let dir = workspace("floor");
    let input = dir.join("p.fasta");
    let output = dir.join("g.fasta");
    write(&input, &format!(">gfp\n{GFP}\n"));
    let out = run(&[
        "--input", input.to_str().unwrap(),
        "--output", output.to_str().unwrap(),
        "--host-table", ECOLI_TABLE,
        "--max-relax", "0.001",
        "--retry-max-relax", "0.001",
        "--quiet",
    ]);
    assert_eq!(out.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("smallest deviation"), "{stderr}");
    assert!(stderr.contains("0.008333"), "it should say what the floor is: {stderr}");
    assert_eq!(
        std::fs::read_to_string(&output).unwrap(),
        "",
        "nothing may be written for a protein that cannot meet the request"
    );
}

#[test]
fn every_candidate_stays_inside_the_requested_gc_band() {
    // Was: five candidates at GC 0.458 under --gc-high 0.40, status low_cai, exit 0.
    // The fallback relaxes the CAI expectation and nothing else.
    let dir = workspace("band");
    let input = dir.join("p.fasta");
    let output = dir.join("g.fasta");
    write(&input, &format!(">gfp\n{GFP}\n"));
    let out = run(&[
        "--input", input.to_str().unwrap(),
        "--output", output.to_str().unwrap(),
        "--host-table", ECOLI_TABLE,
        "--max-relax", "0.01",
        "--retry-max-relax", "0.01",
        "--gc-high", "0.40",
        "--quiet",
    ]);
    assert_eq!(out.status.code(), Some(1), "an unreachable band is a failure");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("GC cannot reach"), "{stderr}");
    assert_eq!(std::fs::read_to_string(&output).unwrap(), "");

    // and with a band it can reach, every sequence written is actually inside it
    let output = dir.join("ok.fasta");
    let out = run(&[
        "--input", input.to_str().unwrap(),
        "--output", output.to_str().unwrap(),
        "--host-table", ECOLI_TABLE,
        "--gc-low", "0.40",
        "--gc-high", "0.50",
        "--quiet",
    ]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    for (header, dna) in records(&std::fs::read_to_string(&output).unwrap()) {
        let gc = dna.bytes().filter(|b| *b == b'G' || *b == b'C').count() as f64
            / dna.len() as f64;
        assert!((0.40..=0.50).contains(&gc), "{header} has GC {gc}");
    }
}

#[test]
fn a_frozen_composition_still_fills_the_candidate_set() {
    // Was: three candidates from four arrangements and exit 1, because one frozen
    // multiset with four restarts cannot produce five distinct sequences.
    let dir = workspace("frozen");
    let input = dir.join("p.fasta");
    let output = dir.join("g.fasta");
    let report_path = dir.join("r.json");
    write(&input, &format!(">gfp\n{GFP}\n"));
    let out = run(&[
        "--input", input.to_str().unwrap(),
        "--output", output.to_str().unwrap(),
        "--host-table", ECOLI_TABLE,
        "--report-json", report_path.to_str().unwrap(),
        "--max-relax", "0",
        "--cai-threshold", "0.7",
        "--quiet",
    ]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let written = records(&std::fs::read_to_string(&output).unwrap());
    assert_eq!(written.len(), 5, "five were asked for and five are obtainable");
    let unique: std::collections::HashSet<&String> = written.iter().map(|(_, d)| d).collect();
    assert_eq!(unique.len(), 5, "and they are distinct, not padded");

    let report = report_of(&report_path);
    let attempt = &report["records"][0]["attempts"][0];
    assert_eq!(report["records"][0]["status"], "ok");
    assert_eq!(attempt["compositions"], 1, "the multiset stayed frozen");
    assert!(
        attempt["rounds"].as_u64().unwrap() > 1,
        "restarts alone could not fill the set, so more rounds were needed"
    );
    assert_eq!(report["records"][0]["attempts"].as_array().unwrap().len(), 1,
               "filling the set is not a retry");
    // frozen means frozen: one deviation across the whole set
    let deviations: std::collections::HashSet<String> = report["records"][0]["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| format!("{:.12}", c["deviation"].as_f64().unwrap()))
        .collect();
    assert_eq!(deviations.len(), 1);
}

#[cfg(unix)]
#[test]
fn two_spellings_of_one_output_path_are_refused_before_anything_is_written() {
    // Was: a symlinked parent made two names for the same not-yet-existing file look
    // different, and the report overwrote the FASTA.
    let dir = workspace("alias");
    let real = dir.join("out");
    std::fs::create_dir(&real).unwrap();
    let link = dir.join("link");
    std::os::unix::fs::symlink(&real, &link).unwrap();
    let input = dir.join("p.fasta");
    write(&input, &format!(">gfp\n{}\n", &GFP[..30]));

    let out = run(&[
        "--input", input.to_str().unwrap(),
        "--output", real.join("genes.fasta").to_str().unwrap(),
        "--host-table", ECOLI_TABLE,
        "--report-json", link.join("genes.fasta").to_str().unwrap(),
        "--quiet",
    ]);
    assert!(!out.status.success(), "the collision must be refused");
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("must be different files"),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        !real.join("genes.fasta").exists(),
        "nothing may be created once the destinations collide"
    );
}
