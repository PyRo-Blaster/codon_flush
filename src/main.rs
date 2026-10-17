//! Command line interface.

use std::path::PathBuf;
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::Parser;
use rayon::prelude::*;

use codon_harmony_core::arrange::constraints::restriction::DEFAULT_ENZYMES;
use codon_harmony_core::arrange::search::{optimise, Options, Outcome};
use codon_harmony_core::composition::{self, Composition};
use codon_harmony_core::fasta::{self, Record};
use codon_harmony_core::host::HostTable;
use codon_harmony_core::{code, metrics};

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
    after_help = "The gene's codon composition is chosen exactly, once; the search then only \
                  rearranges it, so the codon profile and CAI cannot drift while constraints \
                  are repaired."
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

    /// Deviation from the host profile the search may spend to clear a constraint
    #[arg(long, default_value_t = 0.10)]
    max_relax: f64,

    /// Longest run of identical bases allowed
    #[arg(long, default_value_t = 4)]
    local_homopolymer_threshold: usize,

    /// Restriction sites to remove; pass with no values to keep them all
    #[arg(long, num_args = 0.., default_values_t = DEFAULT_ENZYMES.iter().map(|(n, _)| n.to_string()).collect::<Vec<_>>())]
    restriction_enzymes: Vec<String>,

    /// Remove RNA splice sites. Use for mammalian hosts.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    remove_splice_sites: bool,

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

    /// Independent arrangements to try per sequence
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
    if cli.gc_low >= cli.gc_high {
        bail!("--gc-low must be below --gc-high");
    }
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

    let options = Options {
        seed: cli.seed,
        restarts: cli.restarts.max(1),
        homopolymer_limit: cli.local_homopolymer_threshold,
        enzymes: enzymes_from(&cli.restriction_enzymes)?,
        remove_start_sites: cli.remove_start_sites,
        start_site_rule: cli.start_site_rule.into(),
        remove_splice_sites: cli.remove_splice_sites,
        ..Default::default()
    };

    let started = Instant::now();
    let results: Vec<(usize, Result<(Composition, Outcome)>)> = records
        .par_iter()
        .enumerate()
        .map(|(index, record)| {
            let options = Options { seed: options.seed ^ (index as u64 + 1), ..options.clone() };
            (index, optimise_record(record, &host, &options, &cli))
        })
        .collect();

    let mut output = String::new();
    let mut failures = 0usize;
    let mut unfixed = 0usize;
    for (index, result) in results.into_iter().collect::<std::collections::BTreeMap<_, _>>() {
        let record = &records[index];
        match result {
            Ok((composition, outcome)) => {
                output.push_str(&fasta::write_record(&record.header(), &outcome.state.dna()));
                unfixed += outcome.violations as usize;
                if !cli.quiet {
                    report(record, &composition, &outcome, &host);
                } else {
                    println!(
                        "{}\tdeviation {:.4}\tCAI {:.3}\tGC {:.3}\tviolations {}",
                        record.id, outcome.deviation, outcome.cai, outcome.gc, outcome.violations
                    );
                }
            }
            Err(error) => {
                failures += 1;
                eprintln!("{}: {error:#}", record.id);
            }
        }
    }

    std::fs::write(&cli.output, &output)
        .with_context(|| format!("writing {}", cli.output.display()))?;
    if !cli.quiet {
        println!(
            "\n{} of {} sequences written to {} in {:.3} s",
            records.len() - failures,
            records.len(),
            cli.output.display(),
            started.elapsed().as_secs_f64()
        );
    }

    if failures > 0 || unfixed > 0 {
        std::process::exit(1);
    }
    Ok(())
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

fn optimise_record(
    record: &Record,
    host: &HostTable,
    options: &Options,
    cli: &Cli,
) -> Result<(Composition, Outcome)> {
    let composition = composition::build(
        &record.residues,
        host,
        cli.max_relax,
        (cli.gc_low, cli.gc_high),
    )?;
    let outcome = optimise(&record.residues, &composition, host, options);
    if outcome.state.translate() != record.residues {
        bail!("internal error: the optimised gene does not translate back to the protein");
    }
    Ok((composition, outcome))
}

fn report(record: &Record, composition: &Composition, outcome: &Outcome, host: &HostTable) {
    println!("\n{} ({} residues)", record.id, record.residues.len());
    println!(
        "  deviation from host profile  {:.4}{}",
        outcome.deviation,
        if composition.gc_moves > 0 {
            format!(" ({} codons spent reaching the GC band)", composition.gc_moves)
        } else {
            String::new()
        }
    );
    println!("  codon adaptation index       {:.3}", outcome.cai);
    println!("  GC content                   {:.1}%", outcome.gc * 100.0);
    println!(
        "  moves                        {} ({} substitutions, restart {})",
        outcome.moves, outcome.substitutions, outcome.restart
    );
    let unused = metrics::unused_codons(&outcome.state, host);
    if !unused.is_empty() {
        println!(
            "  WARNING codons the host never uses: {}",
            unused.iter().map(|&c| code::codon_string(c)).collect::<Vec<_>>().join(" ")
        );
    }
    if outcome.violations == 0 {
        println!("  constraints                  all satisfied");
    } else {
        println!(
            "  constraints                  {} unresolved: {}",
            outcome.violations,
            outcome
                .by_kind
                .iter()
                .map(|(kind, n)| format!("{} x{n}", kind.label()))
                .collect::<Vec<_>>()
                .join(", ")
        );
        for violation in outcome.remaining.iter().take(10) {
            println!("      base {:>5}: {}", violation.at + 1, violation.detail);
        }
        if outcome.remaining.len() > 10 {
            println!("      ... and {} more", outcome.remaining.len() - 10);
        }
    }
}
