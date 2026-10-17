#!/usr/bin/env python
"""Audit the Rust engine's output with the original Python package.

Reads the protein FASTA it was given, the DNA FASTA the Rust tool produced and, when
offered, the JSON report, then re-derives every claim independently: that each gene
translates back to its source protein, how far its codon usage sits from the host
profile (codon_harmony's own `compare_profiles`), its GC, and its codon adaptation
index. It also checks what the Rust tool says about its own candidate sets: that they
are distinct, ranked by descending CAI, no longer than was asked for, and that the
retry policy the report records is the one that was declared.

    python scripts/validate.py proteins.fasta genes.fasta tests/data/e_coli_316407.json
    python scripts/validate.py proteins.fasta genes.fasta table.json report.json

Each output record carries `source=<original id>` in its description, so several
candidates map back to one protein without the auditor trusting the Rust tool's own
grouping. The CAI printed is Biopython's, whose denominator is `len - 1` where the
Rust engine uses `len`; the two are labelled separately and are expected to differ in
the last digits rather than to match.
"""
import json
import sys
from collections import OrderedDict

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

from codon_harmony.util import codon_use, seq_opt


def source_of(record):
    """The original protein ID a candidate names, or its own ID for old output."""
    for field in record.description.split():
        if field.startswith("source="):
            return field.split("=", 1)[1]
    return record.id


def annotation(record, key):
    for field in record.description.split():
        if field.startswith(key + "="):
            return field.split("=", 1)[1]
    return None


def main(protein_path, dna_path, table_path, report_path=None, threshold=0.10):
    cut, host_profile, cra = codon_use.host_codon_usage("316407", threshold, table_path)
    proteins = OrderedDict(
        (r.id, str(r.seq)) for r in SeqIO.parse(protein_path, "fasta", IUPAC.protein)
    )
    genes = list(SeqIO.parse(dna_path, "fasta", IUPAC.unambiguous_dna))
    report = json.load(open(report_path)) if report_path else None

    by_source = OrderedDict((name, []) for name in proteins)
    orphans = []
    for gene in genes:
        name = source_of(gene)
        if name in by_source:
            by_source[name].append(gene)
        else:
            orphans.append(gene)

    dropped = {c for c, f in host_profile.items() if f == 0}
    ok = True
    print("%-32s %7s %10s %8s %7s  %s"
          % ("candidate", "len_aa", "deviation", "CAI(bio)", "GC", "check"))
    for name, protein in proteins.items():
        candidates = by_source[name]
        if not candidates:
            print("%-32s MISSING from the Rust output" % name)
            ok = False
            continue
        seen, last_cai = set(), None
        for gene in candidates:
            dna = Seq(str(gene.seq), IUPAC.unambiguous_dna)
            problems = []
            if str(dna.translate()) != protein:
                problems.append("TRANSLATION CHANGED")
            if str(dna) in seen:
                problems.append("DUPLICATE SEQUENCE")
            seen.add(str(dna))
            _, deviation = seq_opt.compare_profiles(
                codon_use.count_codons(dna), host_profile, 1.0
            )
            codons = [str(dna[i:i + 3]) for i in range(0, len(dna), 3)]
            used_dropped = sorted({c for c in codons if c in dropped})
            if used_dropped:
                problems.append("uses dropped codons: " + " ".join(used_dropped))
            try:
                cai = cra.cai_for_gene(str(dna))
            except Exception:
                cai = 0.0
            gc = sum(1 for b in str(dna) if b in "GC") / len(dna)

            # the engine's own CAI is annotated on the header; ranking is checked
            # against that, since Biopython's definition is a different number
            claimed = annotation(gene, "cai")
            if claimed is not None:
                claimed = float(claimed)
                if last_cai is not None and claimed > last_cai + 1e-9:
                    problems.append("OUT OF RANK: CAI rises from %.6f" % last_cai)
                last_cai = claimed
            print("%-32s %7d %10.4f %8.3f %7.3f  %s"
                  % (gene.id, len(protein), deviation, cai, gc,
                     "; ".join(problems) if problems else "ok"))
            ok &= not problems

    for gene in orphans:
        print("%-32s names a source protein that is not in the input" % gene.id)
        ok = False

    if report is not None:
        ok &= audit_report(report, by_source)
    print("\n%s" % ("all candidates verified" if ok else "PROBLEMS FOUND"))
    return 0 if ok else 1


def audit_report(report, by_source):
    """Check the report against the FASTA and against the policy it declares."""
    ok = True
    wanted = report["candidates_requested"]
    threshold = report["cai_threshold"]
    policy = report["policy"]
    for record in report["records"]:
        name, candidates = record["id"], record["candidates"]
        written = [str(g.seq) for g in by_source.get(name, [])]
        if [c["dna"] for c in candidates] != written:
            print("%-32s report and FASTA disagree about its candidates" % name)
            ok = False
        if len(candidates) > wanted:
            print("%-32s returned more candidates than asked for" % name)
            ok = False
        if len(candidates) != record["candidate_count"]:
            print("%-32s candidate_count does not match the candidates" % name)
            ok = False
        if len(set(c["dna"] for c in candidates)) != len(candidates):
            print("%-32s padded its set with a repeated sequence" % name)
            ok = False
        best = candidates[0]["cai"] if candidates else None
        if record["low_cai"] != (best is None or best <= threshold):
            print("%-32s low_cai does not match its best CAI" % name)
            ok = False
        attempts = record["attempts"]
        if not 1 <= len(attempts) <= 2:
            print("%-32s ran %d times; the policy allows one retry" % (name, len(attempts)))
            ok = False
        if len(attempts) == 2:
            first_best = attempts[0]["best_cai"]
            if first_best is not None and first_best > threshold:
                print("%-32s retried although its first run qualified" % name)
                ok = False
            for key, declared in (("max_relax", policy["retry"]["max_relax"]),
                                  ("restarts", policy["retry"]["restarts"]),
                                  ("expansion_budget", policy["retry"]["expansion_budget"])):
                if attempts[1][key] != declared:
                    print("%-32s retry %s was %r, not the declared %r"
                          % (name, key, attempts[1][key], declared))
                    ok = False
        elif best is not None and best <= threshold:
            print("%-32s did not retry although nothing exceeded the threshold" % name)
            ok = False
        for candidate in candidates:
            if candidate["deviation"] > candidate["dev_limit"] + 1e-9:
                print("%-32s candidate %d exceeded the deviation bound that applied"
                      % (name, candidate["rank"]))
                ok = False
    return ok


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
