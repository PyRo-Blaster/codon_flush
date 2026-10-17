#!/usr/bin/env python
"""Audit the Rust engine's output with the original Python package.

Reads the protein FASTA it was given and the DNA FASTA the Rust tool produced, then
re-derives every claim independently: that the DNA translates back to the protein, how
far the codon usage sits from the host profile (Biopython / codon_harmony's own
`compare_profiles`), and the codon adaptation index (Biopython).

    python scripts/validate.py proteins.fasta genes.fasta tests/data/e_coli_316407.json
"""
import sys

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

from codon_harmony.util import codon_use, seq_opt


def main(protein_path, dna_path, table_path, threshold=0.10):
    cut, host_profile, cra = codon_use.host_codon_usage("316407", threshold, table_path)
    proteins = {r.id: str(r.seq) for r in SeqIO.parse(protein_path, "fasta", IUPAC.protein)}
    genes = {r.id: str(r.seq) for r in SeqIO.parse(dna_path, "fasta", IUPAC.unambiguous_dna)}

    dropped = {c for c, f in host_profile.items() if f == 0}
    print("%-10s %7s %10s %7s %7s  %s" % ("id", "len_aa", "deviation", "CAI", "GC", "protein"))
    ok = True
    for name, protein in proteins.items():
        if name not in genes:
            print("%-10s MISSING from the Rust output" % name)
            ok = False
            continue
        dna = Seq(genes[name], IUPAC.unambiguous_dna)
        back = str(dna.translate())
        matches = back == protein
        ok &= matches
        _, deviation = seq_opt.compare_profiles(codon_use.count_codons(dna), host_profile, 1.0)
        codons = [str(dna[i:i + 3]) for i in range(0, len(dna), 3)]
        used_dropped = sorted({c for c in codons if c in dropped})
        try:
            cai = cra.cai_for_gene(str(dna))
        except Exception:
            cai = 0.0
        gc = sum(1 for b in str(dna) if b in "GC") / len(dna)
        print("%-10s %7d %10.4f %7.3f %7.3f  %s%s" % (
            name, len(protein), deviation, cai, gc,
            "ok" if matches else "CHANGED",
            "  uses dropped codons: " + " ".join(used_dropped) if used_dropped else "",
        ))
        if used_dropped:
            ok = False
    print("\n%s" % ("all sequences verified" if ok else "PROBLEMS FOUND"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
