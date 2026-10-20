//! Command line interface.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::Parser;
use rayon::prelude::*;

use codon_harmony_core::arrange::constraints::restriction::DEFAULT_ENZYMES;
use codon_harmony_core::arrange::search::Options;
use codon_harmony_core::batch::{self, BatchOptions, RecordReport, RetryPolicy, Status};
use codon_harmony_core::composition::SearchLimits;
use codon_harmony_core::fasta::{self, Record};
use codon_harmony_core::host::HostTable;
use codon_harmony_core::code;

#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
enum SpliceRule {
    NoMotif,
    NoPair,
}

impl From<SpliceRule> for codon_harmony_core::arrange::constraints::splice::Rule {
    fn from(value: SpliceRule) -> Self {
        match value {
            SpliceRule::NoMotif => Self::NoMotif,
            SpliceRule::NoPair => Self::NoPair,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
enum StartSiteRule {
    ShineDalgarno,
    PurineRun,
}

impl From<StartSiteRule> for codon_harmony_core::arrange::constraints::start_site::Rule {
    fn from(value: StartSiteRule) -> Self {
        match value {
            StartSiteRule::PurineRun => Self::PurineRun,
            StartSiteRule::ShineDalgarno => Self::ShineDalgarno,
        }
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "codon_harmony",
    version,
    about = "Reverse translate a protein and optimise the gene for a host's codon usage",
    after_help = "Each protein yields several distinct genes, ranked by codon adaptation \
                  index. If none of them exceeds --cai-threshold the whole search runs a \
                  second time with more deviation allowed and more effort, and the best \
                  candidates of both runs are returned."
)]
struct Cli {
    /// Input FASTA file of amino acid sequences
    #[arg(long)]
    input: PathBuf,

    /// Output FASTA file for the DNA sequences
    #[arg(long)]
    output: PathBuf,

    /// Host codon usage table, as JSON
    #[arg(long, value_name = "FILE")]
    host_table: PathBuf,

    /// Lowest codon fraction per amino acid the host is allowed to use
    #[arg(long, default_value_t = 0.10)]
    host_threshold: f64,

    /// Distinct genes to return per protein, best CAI first
    #[arg(long, default_value_t = 5)]
    candidates: usize,

    /// CAI a candidate must strictly exceed, or the search is run a second time.
    /// Equality does not pass.
    #[arg(long, default_value_t = 0.80)]
    cai_threshold: f64,

    /// Deviation from the host profile the search may spend to clear a constraint
    #[arg(long, default_value_t = 0.10)]
    max_relax: f64,

    /// Deviation allowed on the retry. Defaults to min(1, --max-relax + 0.10) for a
    /// positive limit, and to 0 for a frozen composition. Must be at least --max-relax.
    #[arg(long)]
    retry_max_relax: Option<f64>,

    /// Arrangements per composition on the retry. Defaults to twice --restarts.
    #[arg(long)]
    retry_restarts: Option<usize>,

    /// States the composition search may expand. The retry gets twice as many.
    #[arg(long, default_value_t = 2000)]
    composition_budget: usize,

    /// Compositions carried from one round of the composition search to the next
    #[arg(long, default_value_t = 12)]
    composition_beam: usize,

    /// Distinct codon multisets handed on to the arrangement search
    #[arg(long, default_value_t = 8)]
    max_compositions: usize,

    /// Write a machine-readable report of every candidate and both runs
    #[arg(long, value_name = "FILE")]
    report_json: Option<PathBuf>,

    /// Longest run of identical bases allowed
    #[arg(long, default_value_t = 4)]
    local_homopolymer_threshold: usize,

    /// Restriction sites to remove; pass with no values to keep them all
    #[arg(long, num_args = 0.., default_values_t = DEFAULT_ENZYMES.iter().map(|(n, _)| n.to_string()).collect::<Vec<_>>())]
    restriction_enzymes: Vec<String>,

    /// Remove RNA splice sites. Use for mammalian hosts.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    remove_splice_sites: bool,

    /// What counts as a splice problem. "no-motif", the default, rejects every
    /// donor and every acceptor. "no-pair" rejects only a donor and an acceptor
    /// together, on the grounds that splicing needs both: a lone half site cannot
    /// define an intron within the coding sequence. Some proteins leave no choice —
    /// 30 tripeptides, all of the form [MW][CWY][ADEGV], force a donor in every
    /// synonymous encoding — and "no-pair" is what lets those be built at all.
    /// CAVEAT: a lone donor is inert only relative to the coding sequence. A
    /// transgene is expressed from a vector, and a donor inside the insert can still
    /// pair with an acceptor downstream in the vector's expression cassette, which
    /// is the usual mammalian arrangement. Every retained site is reported with its
    /// position so the gene can be checked against that vector.
    #[arg(long, value_enum, default_value_t = SpliceRule::NoMotif)]
    splice_rule: SpliceRule,

    /// Remove alternate start sites. Use for bacterial hosts.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    remove_start_sites: bool,

    /// How to recognise an internal start site. "shine-dalgarno" asks whether the
    /// sequence could actually pair with the 16S rRNA, and looks at ATG, GTG and TTG.
    /// "purine-run" is the original rule: five purines before TTG, CTG or ATG.
    #[arg(long, value_enum, default_value_t = StartSiteRule::ShineDalgarno)]
    start_site_rule: StartSiteRule,

    /// Whole-gene GC bounds
    #[arg(long, default_value_t = 0.30)]
    gc_low: f64,
    #[arg(long, default_value_t = 0.65)]
    gc_high: f64,

    /// Random seed; the same seed reproduces the same genes
    #[arg(long, default_value_t = 1)]
    seed: u64,

    /// Independent arrangements to try per composition
    #[arg(long, default_value_t = 4)]
    restarts: usize,

    /// Worker threads; 0 means one per core
    #[arg(long, default_value_t = 0)]
    threads: usize,

    /// Accepted for compatibility with the old sampling loop, and ignored
    #[arg(long)]
    cycles: Option<usize>,
    #[arg(long)]
    inner_cycles: Option<usize>,

    /// Print a line per sequence rather than a block
    #[arg(long)]
    quiet: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    for (flag, value) in [("--cycles", cli.cycles), ("--inner-cycles", cli.inner_cycles)] {
        if value.is_some() {
            eprintln!(
                "note: {flag} describes the old sampling loop and no longer does anything; \
                 use --restarts to try more arrangements"
            );
        }
    }
    let options = check(&cli)?;
    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .context("setting up the thread pool")?;
    }

    let host = HostTable::load(&cli.host_table, cli.host_threshold)?;
    let dropped = host.dropped();
    if !dropped.is_empty() && !cli.quiet {
        println!(
            "host threshold {:.2} drops {} codons: {}",
            cli.host_threshold,
            dropped.len(),
            dropped.iter().map(|&c| code::codon_string(c)).collect::<Vec<_>>().join(" ")
        );
    }

    let text = std::fs::read_to_string(&cli.input)
        .with_context(|| format!("reading {}", cli.input.display()))?;
    let records = fasta::parse(&text).with_context(|| format!("in {}", cli.input.display()))?;

    if !cli.quiet {
        announce(&cli, &options);
    }

    let started = Instant::now();
    let reports: BTreeMap<usize, RecordReport> = records
        .par_iter()
        .enumerate()
        .map(|(index, record)| {
            (index, batch::run_record(index, &record.residues, &host, &options))
        })
        .collect();

    let mut output = String::new();
    let mut failures = 0usize;
    for (index, report) in &reports {
        let record = &records[*index];
        for (rank, candidate) in report.candidates.iter().enumerate() {
            output.push_str(&fasta::write_record(
                &header_for(record, *index, rank + 1, report, candidate),
                &candidate.dna,
            ));
        }
        if report.status.is_failure() {
            failures += 1;
        }
        tell(record, report, &cli);
    }

    write_atomically(&cli.output, output.as_bytes())
        .with_context(|| format!("writing {}", cli.output.display()))?;
    if let Some(path) = &cli.report_json {
        let json = report_json(&cli, &options, &records, &reports, failures);
        write_atomically(path, format!("{json:#}\n").as_bytes())
            .with_context(|| format!("writing {}", path.display()))?;
    }

    if !cli.quiet {
        let written: usize = reports.values().map(|r| r.candidates.len()).sum();
        println!(
            "\n{} candidates for {} of {} proteins written to {} in {:.3} s",
            written,
            reports.values().filter(|r| !r.candidates.is_empty()).count(),
            records.len(),
            cli.output.display(),
            started.elapsed().as_secs_f64()
        );
    }

    if failures > 0 {
        std::process::exit(1);
    }
    Ok(())
}

/// Validate the settings and turn them into the batch options. Everything numeric is
/// checked here so no run starts on a limit that cannot mean anything.
fn check(cli: &Cli) -> Result<BatchOptions> {
    for (flag, value) in [
        ("--host-threshold", cli.host_threshold),
        ("--cai-threshold", cli.cai_threshold),
        ("--max-relax", cli.max_relax),
        ("--gc-low", cli.gc_low),
        ("--gc-high", cli.gc_high),
    ] {
        if !value.is_finite() {
            bail!("{flag} must be a finite number");
        }
    }
    if cli.gc_low >= cli.gc_high {
        bail!("--gc-low must be below --gc-high");
    }
    if !(0.0..1.0).contains(&cli.cai_threshold) {
        bail!("--cai-threshold must be at least 0 and below 1");
    }
    if !(0.0..=1.0).contains(&cli.max_relax) {
        bail!("--max-relax must be between 0 and 1");
    }
    if cli.candidates == 0 {
        bail!("--candidates must be at least 1");
    }
    if cli.restarts == 0 {
        bail!("--restarts must be at least 1");
    }
    if cli.composition_budget == 0 || cli.composition_beam == 0 || cli.max_compositions == 0 {
        bail!("--composition-budget, --composition-beam and --max-compositions must be positive");
    }
    let default_retry =
        RetryPolicy::doubling(cli.max_relax, cli.restarts, cli.composition_budget);
    let retry_max_relax = match cli.retry_max_relax {
        None => default_retry.max_relax,
        Some(value) => {
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                bail!("--retry-max-relax must be a finite number between 0 and 1");
            }
            if value < cli.max_relax {
                bail!(
                    "--retry-max-relax ({value}) must be at least --max-relax ({})",
                    cli.max_relax
                );
            }
            value
        }
    };
    let retry_restarts = match cli.retry_restarts {
        None => default_retry.restarts,
        Some(0) => bail!("--retry-restarts must be at least 1"),
        Some(value) => value,
    };
    let Some(retry_budget) = cli.composition_budget.checked_mul(2) else {
        bail!("--composition-budget is too large to double for the retry");
    };

    let mut paths: Vec<(&Path, &str)> =
        vec![(cli.input.as_path(), "--input"), (cli.output.as_path(), "--output")];
    if let Some(report) = &cli.report_json {
        paths.push((report.as_path(), "--report-json"));
    }
    for i in 0..paths.len() {
        for j in i + 1..paths.len() {
            if same_path(paths[i].0, paths[j].0) {
                bail!("{} and {} must be different files", paths[i].1, paths[j].1);
            }
        }
    }

    Ok(BatchOptions {
        wanted: cli.candidates,
        cai_threshold: cli.cai_threshold,
        max_relax: cli.max_relax,
        gc_band: (cli.gc_low, cli.gc_high),
        limits: SearchLimits {
            beam_width: cli.composition_beam,
            max_expansions: cli.composition_budget,
            max_compositions: cli.max_compositions,
        },
        retry: RetryPolicy {
            max_relax: retry_max_relax,
            restarts: retry_restarts,
            expansions: retry_budget,
        },
        arrange: Options {
            seed: cli.seed,
            restarts: cli.restarts,
            gc_band: (cli.gc_low, cli.gc_high),
            homopolymer_limit: cli.local_homopolymer_threshold,
            enzymes: enzymes_from(&cli.restriction_enzymes)?,
            remove_start_sites: cli.remove_start_sites,
            start_site_rule: cli.start_site_rule.into(),
            remove_splice_sites: cli.remove_splice_sites,
            splice_rule: cli.splice_rule.into(),
            ..Default::default()
        },
    })
}

/// Where a path will actually land, as far as can be known before the file exists.
///
/// Canonicalising the whole path fails for a file that is not there yet, which is the
/// normal case for an output. Canonicalising the *parent* succeeds, which is what
/// matters: two spellings of a directory that differ only through a symlink resolve
/// to the same place, so `out/genes.fasta` and `link/genes.fasta` compare equal
/// before either file is created.
fn destination(path: &Path) -> PathBuf {
    let parent = match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => Path::new("."),
    };
    let resolved = parent.canonicalize().unwrap_or_else(|_| parent.to_path_buf());
    match path.file_name() {
        Some(name) => resolved.join(name),
        None => resolved,
    }
}

/// Two settings naming the same file.
fn same_path(a: &Path, b: &Path) -> bool {
    a == b || destination(a) == destination(b)
}

/// Write through a temporary sibling and rename. Each file lands whole or not at all;
/// several files are not one transaction, and a later failure leaves earlier ones in
/// place.
fn write_atomically(path: &Path, bytes: &[u8]) -> Result<()> {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "output".to_string());
    let temp = path.with_file_name(format!(".{name}.{}.tmp", std::process::id()));
    std::fs::write(&temp, bytes).with_context(|| format!("writing {}", temp.display()))?;
    if let Err(error) = std::fs::rename(&temp, path) {
        let _ = std::fs::remove_file(&temp);
        return Err(error)
            .with_context(|| format!("renaming {} to {}", temp.display(), path.display()));
    }
    Ok(())
}

fn announce(cli: &Cli, options: &BatchOptions) {
    println!(
        "asking for {} candidates per protein; a run stands when one exceeds CAI {:.6}",
        options.wanted, options.cai_threshold
    );
    println!(
        "  first run  deviation <= {:.4}  restarts {}  composition budget {}",
        options.max_relax, cli.restarts, options.limits.max_expansions
    );
    println!(
        "  retry      deviation <= {:.4}  restarts {}  composition budget {}{}",
        options.retry.max_relax,
        options.retry.restarts,
        options.retry.expansions,
        if options.max_relax == 0.0 && options.retry.max_relax == 0.0 {
            "  (composition frozen; only arrangement effort grows)"
        } else {
            ""
        }
    );
}

fn enzymes_from(names: &[String]) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    for name in names {
        match DEFAULT_ENZYMES.iter().find(|(n, _)| n.eq_ignore_ascii_case(name)) {
            Some((n, site)) => out.push((n.to_string(), site.to_string())),
            None => bail!(
                "unknown restriction enzyme {name:?}; known ones are {}",
                DEFAULT_ENZYMES.iter().map(|(n, _)| *n).collect::<Vec<_>>().join(", ")
            ),
        }
    }
    Ok(out)
}

/// A collision-free output ID. Two input records may share an ID; the input ordinal
/// never repeats, so it, not the name, identifies the source.
fn output_id(index: usize, rank: usize) -> String {
    format!("record_{:06}_candidate_{:02}", index + 1, rank)
}

fn header_for(
    record: &Record,
    index: usize,
    rank: usize,
    report: &RecordReport,
    candidate: &batch::Candidate,
) -> String {
    let mut header = format!(
        "{} rank={}/{} cai={:.6} deviation={:.6} gc={:.4} run={} composition={} restart={}",
        output_id(index, rank),
        rank,
        report.candidates.len(),
        candidate.cai,
        candidate.deviation,
        candidate.gc,
        candidate.run,
        candidate.composition,
        candidate.restart,
    );
    if report.low_cai {
        header.push_str(" low_cai=true");
    }
    // A gene the pair rule let through keeps its half sites, and the header is
    // where someone pasting this into a vector will see them.
    if !candidate.retained_splice_sites.is_empty() {
        header.push_str(&format!(
            " retained_splice_sites={}",
            describe_sites(&candidate.retained_splice_sites)
        ));
    }
    header.push_str(&format!(" source={}", record.id));
    if !record.description.is_empty() {
        header.push(' ');
        header.push_str(&record.description);
    }
    header
}

/// Say what happened to one protein. Warnings go to stderr so they survive `--quiet`
/// and a redirected stdout.
fn tell(record: &Record, report: &RecordReport, cli: &Cli) {
    let best = report.best_cai();
    if let Some(error) = &report.error {
        eprintln!("{}: {error}", record.id);
    }
    match report.status {
        Status::NoValidCandidates if report.error.is_none() => eprintln!(
            "{}: no valid candidate after {} run(s); nothing written for this protein{}",
            record.id,
            report.attempts.len(),
            why_none(report)
        ),
        Status::InsufficientCandidates => eprintln!(
            "{}: insufficient_candidates: {} distinct valid sequence(s), {} asked for{}",
            record.id,
            report.candidates.len(),
            cli.candidates,
            if report.search_exhausted {
                "; the composition search hit its budget, so other encodings may exist"
            } else {
                ""
            }
        ),
        _ => {}
    }
    if report.low_cai && !report.candidates.is_empty() {
        eprintln!(
            "{}: neither run produced CAI > {:.2}; returning {} ranked candidates (best CAI {:.6})",
            record.id,
            cli.cai_threshold,
            report.candidates.len(),
            best.unwrap_or(0.0)
        );
    }

    if cli.quiet {
        println!(
            "{}\tstatus {}\tcandidates {}\tbest CAI {}\truns {}",
            record.id,
            report.status.label(),
            report.candidates.len(),
            best.map(|c| format!("{c:.6}")).unwrap_or_else(|| "-".to_string()),
            report.attempts.len()
        );
        return;
    }

    println!("\n{} ({} residues)", record.id, record.residues.len());
    println!("  status                       {}", report.status.label());
    if let Some(reason) = &report.retry_reason {
        println!("  retried because              {reason}");
    }
    for attempt in &report.attempts {
        println!(
            "  run {}  deviation <= {:.4}  restarts {}  compositions {}  expansions {}/{}{}",
            attempt.run,
            attempt.max_relax,
            attempt.restarts,
            attempt.compositions,
            attempt.expansions,
            attempt.expansion_budget,
            if attempt.search_exhausted { "  (budget spent)" } else { "" }
        );
        println!(
            "         arrangements {} in {} round(s)  valid {}  best CAI {}",
            attempt.arrangements,
            attempt.rounds,
            attempt.valid,
            attempt.best_cai.map(|c| format!("{c:.6}")).unwrap_or_else(|| "-".to_string())
        );
        let rejected = &attempt.rejected;
        if rejected.unresolved + rejected.over_limit + rejected.unsound > 0 {
            println!(
                "         rejected {} with unresolved constraints{}, {} over a hard \
                 limit, {} unsound",
                rejected.unresolved,
                match (&rejected.fewest_violations, &rejected.violation_detail) {
                    (Some(n), Some(detail)) => format!(" (fewest {n}: {detail})"),
                    _ => String::new(),
                },
                rejected.over_limit,
                rejected.unsound
            );
        }
    }
    for (rank, candidate) in report.candidates.iter().enumerate() {
        println!(
            "  {:>2}. {}  CAI {:.6}  deviation {:.4}  GC {:.1}%  run {}  composition {}  restart {}",
            rank + 1,
            output_id(report.index, rank + 1),
            candidate.cai,
            candidate.deviation,
            candidate.gc * 100.0,
            candidate.run,
            candidate.composition,
            candidate.restart
        );
        if !candidate.retained_splice_sites.is_empty() {
            println!(
                "      keeps {} splice site(s) no encoding can avoid: {}",
                candidate.retained_splice_sites.len(),
                describe_sites(&candidate.retained_splice_sites)
            );
        }
    }
}

/// Retained splice sites as `donor@124,acceptor@871`, in position order.
fn describe_sites(sites: &[codon_harmony_core::arrange::constraints::splice::Site]) -> String {
    sites
        .iter()
        .map(|site| format!("{}@{}", site.kind(), site.at + 1))
        .collect::<Vec<_>>()
        .join(",")
}

/// The shortest honest account of why a protein produced nothing, from the run that
/// came closest.
fn why_none(report: &RecordReport) -> String {
    let best = report
        .attempts
        .iter()
        .filter_map(|a| {
            a.rejected
                .fewest_violations
                .map(|n| (n, a.rejected.violation_detail.clone().unwrap_or_default()))
        })
        .min_by_key(|(n, _)| *n);
    match best {
        Some((n, detail)) if !detail.is_empty() => {
            format!("; the closest arrangement still left {n} violations: {detail}")
        }
        Some((n, _)) => format!("; the closest arrangement still left {n} violations"),
        None => String::new(),
    }
}

/// The CAI definition the engine uses, recorded with every report so a number can
/// always be reproduced.
const CAI_CONVENTION: &str = "w_c = f_c / max(f_synonyms) after thresholding and \
    renormalisation; CAI = exp(sum(n_c * ln w_c) / N) over codons whose amino acid has \
    more than one codon, so Met and Trp are excluded and stops are included; CAI is 1 \
    when no codon counts. Legacy Biopython divides by N - 1 and will differ slightly.";

fn report_json(
    cli: &Cli,
    options: &BatchOptions,
    records: &[Record],
    reports: &BTreeMap<usize, RecordReport>,
    failures: usize,
) -> serde_json::Value {
    use serde_json::json;
    let records_json: Vec<serde_json::Value> = reports
        .iter()
        .map(|(index, report)| {
            let record = &records[*index];
            json!({
                "index": index + 1,
                "id": record.id,
                "description": record.description,
                "header": record.header(),
                "residues": record.residues.len(),
                "status": report.status.label(),
                "low_cai": report.low_cai,
                "search_exhausted": report.search_exhausted,
                "cai_threshold": options.cai_threshold,
                "candidates_requested": options.wanted,
                "candidate_count": report.candidates.len(),
                "retry_reason": report.retry_reason,
                "error": report.error,
                "attempts": report.attempts.iter().map(|attempt| json!({
                    "run": attempt.run,
                    "max_relax": attempt.max_relax,
                    "dev_limit": attempt.dev_limit,
                    "restarts": attempt.restarts,
                    "expansion_budget": attempt.expansion_budget,
                    "compositions": attempt.compositions,
                    "expansions": attempt.expansions,
                    "search_exhausted": attempt.search_exhausted,
                    "arrangements": attempt.arrangements,
                    "rounds": attempt.rounds,
                    "valid": attempt.valid,
                    "best_cai": attempt.best_cai,
                    "rejected": {
                        "unresolved_constraints": attempt.rejected.unresolved,
                        "fewest_violations": attempt.rejected.fewest_violations,
                        "violation_detail": attempt.rejected.violation_detail,
                        "over_hard_limit": attempt.rejected.over_limit,
                        "unsound": attempt.rejected.unsound,
                    },
                })).collect::<Vec<_>>(),
                "candidates": report.candidates.iter().enumerate().map(|(rank, c)| json!({
                    "rank": rank + 1,
                    "output_id": output_id(*index, rank + 1),
                    "dna": c.dna,
                    "cai": c.cai,
                    "deviation": c.deviation,
                    "gc": c.gc,
                    "run": c.run,
                    "composition": c.composition,
                    "restart": c.restart,
                    "seed": c.seed,
                    "max_relax": c.max_relax,
                    "dev_limit": c.dev_limit,
                    "moves": c.moves,
                    "substitutions": c.substitutions,
                    "retained_splice_sites": c.retained_splice_sites.iter().map(|site| json!({
                        "at": site.at + 1,
                        "kind": site.kind(),
                        "sequence": site.text,
                    })).collect::<Vec<_>>(),
                })).collect::<Vec<_>>(),
            })
        })
        .collect();

    json!({
        "version": env!("CARGO_PKG_VERSION"),
        "cai_convention": CAI_CONVENTION,
        "seed": cli.seed,
        "host_table": cli.host_table.display().to_string(),
        "host_threshold": cli.host_threshold,
        "gc_band": [cli.gc_low, cli.gc_high],
        "splice_rule": match cli.splice_rule {
            SpliceRule::NoMotif => "no-motif",
            SpliceRule::NoPair => "no-pair",
        },
        "candidates_requested": options.wanted,
        "cai_threshold": options.cai_threshold,
        "policy": {
            "first": {
                "max_relax": options.max_relax,
                "restarts": options.arrange.restarts,
                "expansion_budget": options.limits.max_expansions,
            },
            "retry": {
                "max_relax": options.retry.max_relax,
                "restarts": options.retry.restarts,
                "expansion_budget": options.retry.expansions,
            },
            "composition_beam": options.limits.beam_width,
            "max_compositions": options.limits.max_compositions,
        },
        "summary": {
            "proteins": records.len(),
            "ok": reports.values().filter(|r| r.status == Status::Ok).count(),
            "low_cai": reports.values().filter(|r| r.status == Status::LowCai).count(),
            "insufficient_candidates":
                reports.values().filter(|r| r.status == Status::InsufficientCandidates).count(),
            "no_valid_candidates":
                reports.values().filter(|r| r.status == Status::NoValidCandidates).count(),
            "errors": reports.values().filter(|r| r.status == Status::Error).count(),
            "failed_records": failures,
            "exit_code": i32::from(failures > 0),
        },
        "records": records_json,
    })
}
