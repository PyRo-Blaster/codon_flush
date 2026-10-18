#!/usr/bin/env python
"""Audit the Rust engine's output independently of what it says about itself.

Reads the protein FASTA it was given, the DNA FASTA the Rust tool produced and, when
offered, the JSON report. Every number is re-derived from the sequences themselves:
translation, whole-gene GC, deviation from the host profile (the original Python
package's own `compare_profiles`) and the codon adaptation index, recomputed here from
the host table under the convention the engine documents. The claimed values in the
FASTA descriptions and in the report are then checked *against* those, so a wrong CAI,
a wrong ranking or an out-of-band sequence fails the audit rather than passing it.

    python scripts/validate.py proteins.fasta genes.fasta tests/data/e_coli_316407.json
    python scripts/validate.py proteins.fasta genes.fasta table.json report.json

Each output record carries `source=<original id>` in its description, so several
candidates map back to one protein without the auditor trusting the Rust tool's own
grouping.

Biopython's CAI divides the log sum by `len - 1` where the engine uses `len`; it is
reported in a separate column and is never what a claim is checked against.

The metric and audit helpers deliberately import nothing from Biopython, so they can
be exercised on falsified input by `tests/test_validate_audit.py` without it.
"""
import json
import sys
from collections import OrderedDict

#: Tolerance for comparing a claimed float against a recomputed one.
CLAIM_TOLERANCE = 1e-6
#: Tolerance for a hard bound, matching `composition::LIMIT_TOLERANCE` on the Rust side.
LIMIT_TOLERANCE = 1e-9

BASES = "TCAG"
#: The standard genetic code, in the codon order TTT, TTC, TTA ... GGG.
AMINO_ACIDS = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)
CODON_TO_AA = {
    a + b + c: AMINO_ACIDS[i * 16 + j * 4 + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)
}
#: How many codons the genetic code gives each amino acid. Met and Trp have one, so
#: they carry no information about codon choice and the CAI skips them.
GROUP_SIZE = {}
for _codon, _aa in CODON_TO_AA.items():
    GROUP_SIZE[_aa] = GROUP_SIZE.get(_aa, 0) + 1


def host_weights(table, threshold=0.10):
    """Relative adaptiveness per codon, by the engine's stated rules.

    Accepts the nested amino-acid shape and the flat codon map, normalises RNA
    spellings, renormalises each synonymous group, zeroes codons below the threshold,
    renormalises again, and divides by the group maximum.
    """
    freq = {codon: 0.0 for codon in CODON_TO_AA}
    for key, entry in table.items():
        if isinstance(entry, dict):
            for codon, value in entry.items():
                freq[codon.upper().replace("U", "T")] = float(value)
        else:
            freq[key.upper().replace("U", "T")] = float(entry)

    groups = {}
    for codon, aa in CODON_TO_AA.items():
        groups.setdefault(aa, []).append(codon)

    def renormalise():
        for codons in groups.values():
            total = sum(freq[c] for c in codons)
            if total <= 0:
                raise ValueError("host table has no usable codon for an amino acid")
            for c in codons:
                freq[c] /= total

    renormalise()
    for codon in freq:
        if freq[codon] < threshold:
            freq[codon] = 0.0
    renormalise()

    weights = {}
    for codons in groups.values():
        top = max(freq[c] for c in codons)
        for c in codons:
            weights[c] = freq[c] / top if top > 0 else 0.0
    return freq, weights


def codons_of(dna):
    return [dna[i:i + 3] for i in range(0, len(dna), 3)]


def engine_cai(dna, weights):
    """CAI under the engine's convention: geometric mean over codons whose amino acid
    has more than one codon, so Met and Trp are skipped and stops are included; a
    codon the host never uses collapses it to zero; no counted codon gives 1."""
    import math

    total, counted = 0.0, 0
    for codon in codons_of(dna):
        aa = CODON_TO_AA[codon]
        if GROUP_SIZE[aa] < 2:
            continue
        weight = weights[codon]
        if weight <= 0:
            return 0.0
        total += math.log(weight)
        counted += 1
    if counted == 0:
        return 1.0
    return math.exp(total / counted)


def gc_fraction(dna):
    return sum(1 for b in dna if b in "GC") / len(dna) if dna else 0.0


def source_of_description(description):
    """The original protein ID a candidate names, or None for unannotated output."""
    for field in description.split():
        if field.startswith("source="):
            return field.split("=", 1)[1]
    return None


def annotation(description, key):
    for field in description.split():
        if field.startswith(key + "="):
            return field.split("=", 1)[1]
    return None


def close(claimed, actual, tolerance=CLAIM_TOLERANCE):
    return claimed is not None and abs(claimed - actual) <= tolerance


def check_candidate(dna, description, weights):
    """Recompute a candidate's metrics and check the FASTA's own claims against them.

    Returns (metrics, problems).
    """
    metrics = {"cai": engine_cai(dna, weights), "gc": gc_fraction(dna)}
    problems = []
    claimed_cai = annotation(description, "cai")
    if claimed_cai is not None and not close(float(claimed_cai), metrics["cai"]):
        problems.append(
            "CLAIMED CAI %s, recomputed %.6f" % (claimed_cai, metrics["cai"])
        )
    claimed_gc = annotation(description, "gc")
    if claimed_gc is not None and not close(float(claimed_gc), metrics["gc"], 1e-4):
        problems.append("CLAIMED GC %s, recomputed %.4f" % (claimed_gc, metrics["gc"]))
    return metrics, problems


def audit_report(report, genes_by_source, weights, deviations):
    """Check the report against the sequences, not against itself.

    `genes_by_source` maps a source protein ID to the list of DNA strings written for
    it, in file order. `deviations` maps a DNA string to its independently computed
    deviation from the host profile.
    """
    problems = []
    wanted = report["candidates_requested"]
    threshold = report["cai_threshold"]
    policy = report["policy"]
    band = tuple(report["gc_band"])

    def fault(name, message):
        problems.append("%s: %s" % (name, message))

    for record in report["records"]:
        name, candidates = record["id"], record["candidates"]
        written = genes_by_source.get(name, [])
        if [c["dna"] for c in candidates] != written:
            fault(name, "report and FASTA disagree about its candidates")
        if len(candidates) > wanted:
            fault(name, "returned more candidates than asked for")
        if len(candidates) != record["candidate_count"]:
            fault(name, "candidate_count does not match the candidates")
        if len({c["dna"] for c in candidates}) != len(candidates):
            fault(name, "padded its set with a repeated sequence")

        # recomputed metrics, then every claim checked against them
        actual = []
        for candidate in candidates:
            dna = candidate["dna"]
            cai = engine_cai(dna, weights)
            gc = gc_fraction(dna)
            actual.append(cai)
            rank = candidate["rank"]
            if not close(candidate["cai"], cai):
                fault(name, "candidate %d claims CAI %.6f, recomputed %.6f"
                      % (rank, candidate["cai"], cai))
            if not close(candidate["gc"], gc, 1e-6):
                fault(name, "candidate %d claims GC %.6f, recomputed %.6f"
                      % (rank, candidate["gc"], gc))
            if dna in deviations and not close(candidate["deviation"], deviations[dna], 1e-4):
                fault(name, "candidate %d claims deviation %.6f, recomputed %.6f"
                      % (rank, candidate["deviation"], deviations[dna]))
            # hard bounds, against the request rather than against the claim
            if not band[0] - LIMIT_TOLERANCE <= gc <= band[1] + LIMIT_TOLERANCE:
                fault(name, "candidate %d has GC %.4f, outside the requested %.2f-%.2f"
                      % (rank, gc, band[0], band[1]))
            bound = max(candidate["dev_limit"], candidate["max_relax"])
            measured = deviations.get(dna, candidate["deviation"])
            if measured > bound + 1e-4:
                fault(name, "candidate %d deviates %.6f, above the %.6f that applied"
                      % (rank, measured, bound))

        # ranking, by the recomputed values
        for i in range(1, len(actual)):
            if actual[i] > actual[i - 1] + CLAIM_TOLERANCE:
                fault(name, "candidates are not in descending CAI: %.6f then %.6f"
                      % (actual[i - 1], actual[i]))
        best = max(actual) if actual else None
        if record["low_cai"] != (best is None or best <= threshold):
            fault(name, "low_cai does not match the recomputed best CAI")

        attempts = record["attempts"]
        if not 1 <= len(attempts) <= 2:
            fault(name, "ran %d times; the policy allows one retry" % len(attempts))
        if len(attempts) == 2:
            first = attempts[0]["best_cai"]
            if first is not None and first > threshold:
                fault(name, "retried although its first run qualified")
            for key, declared in (("max_relax", policy["retry"]["max_relax"]),
                                  ("restarts", policy["retry"]["restarts"]),
                                  ("expansion_budget", policy["retry"]["expansion_budget"])):
                if attempts[1][key] != declared:
                    fault(name, "retry %s was %r, not the declared %r"
                          % (key, attempts[1][key], declared))
        elif best is not None and best <= threshold:
            fault(name, "did not retry although nothing exceeded the threshold")
    return problems


def main(protein_path, dna_path, table_path, report_path=None, threshold=0.10):
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    from Bio import SeqIO

    from codon_harmony.util import codon_use, seq_opt

    threshold = float(threshold)
    cut, host_profile, cra = codon_use.host_codon_usage("316407", threshold, table_path)
    _, weights = host_weights(json.load(open(table_path)), threshold)

    proteins = OrderedDict(
        (r.id, str(r.seq)) for r in SeqIO.parse(protein_path, "fasta", IUPAC.protein)
    )
    genes = list(SeqIO.parse(dna_path, "fasta", IUPAC.unambiguous_dna))
    report = json.load(open(report_path)) if report_path else None

    by_source = OrderedDict((name, []) for name in proteins)
    orphans, deviations = [], {}
    for gene in genes:
        name = source_of_description(gene.description) or gene.id
        (by_source[name] if name in by_source else orphans).append(gene)

    dropped = {c for c, f in host_profile.items() if f == 0}
    ok = True
    print("%-32s %7s %10s %9s %9s %7s  %s"
          % ("candidate", "len_aa", "deviation", "CAI", "CAI(bio)", "GC", "check"))
    for name, protein in proteins.items():
        candidates = by_source[name]
        if not candidates:
            print("%-32s MISSING from the Rust output" % name)
            ok = False
            continue
        seen = set()
        for gene in candidates:
            dna = str(gene.seq)
            metrics, problems = check_candidate(dna, gene.description, weights)
            if str(Seq(dna, IUPAC.unambiguous_dna).translate()) != protein:
                problems.append("TRANSLATION CHANGED")
            if dna in seen:
                problems.append("DUPLICATE SEQUENCE")
            seen.add(dna)
            _, deviation = seq_opt.compare_profiles(
                codon_use.count_codons(Seq(dna, IUPAC.unambiguous_dna)), host_profile, 1.0
            )
            deviations[dna] = deviation
            used_dropped = sorted(set(codons_of(dna)) & dropped)
            if used_dropped:
                problems.append("uses dropped codons: " + " ".join(used_dropped))
            try:
                legacy = cra.cai_for_gene(dna)
            except Exception:
                legacy = 0.0
            print("%-32s %7d %10.4f %9.6f %9.3f %7.3f  %s"
                  % (gene.id, len(protein), deviation, metrics["cai"], legacy,
                     metrics["gc"], "; ".join(problems) if problems else "ok"))
            ok &= not problems

    for gene in orphans:
        print("%-32s names a source protein that is not in the input" % gene.id)
        ok = False

    if report is not None:
        written = OrderedDict(
            (name, [str(g.seq) for g in genes]) for name, genes in by_source.items()
        )
        for problem in audit_report(report, written, weights, deviations):
            print(problem)
            ok = False
    print("\n%s" % ("all candidates verified" if ok else "PROBLEMS FOUND"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
