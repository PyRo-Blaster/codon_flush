#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Regression tests for the sequence cleanup and harmonization steps.

These tests use a locally stored host codon usage table so that they never
reach out to the network.
"""


import os
import random
import unittest

from unittest import mock

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from codon_harmony import codon_harmony
from codon_harmony.util import codon_use, seq_opt

# locally stored host table -- "Escherichia coli" (NCBI taxonomy ID 316407)
HOST_ID = "316407"
HOST_TABLE = os.path.join(os.path.dirname(__file__), "data", "e_coli_316407.json")
HOST_THRESHOLD = 0.10

# the cleanup steps are stochastic, so each one is checked with several
# independent random seeds
SEEDS = range(20)


def _codons(dna_sequence):
    """Split a DNA sequence into a list of codons."""
    sequence = str(dna_sequence)
    return [sequence[i : i + 3] for i in range(0, len(sequence), 3)]


def _homopolymer_violations(dna_sequence, n_codons=2, homopolymer_threshold=4):
    """Find the windows in which a single nucleotide appears more often than
    `homopolymer_threshold` times.

    The windows match the ones used by
    `codon_harmony.util.seq_opt.remove_local_homopolymers`.
    """
    sequence = str(dna_sequence)
    violations = []
    for i in range(0, len(sequence), 3):
        end = i + (n_codons * 3)
        window = sequence[i:end] if end < len(sequence) else sequence[i:]
        if not window:
            continue
        counts = {letter: window.count(letter) for letter in set(window)}
        letter = max(counts, key=lambda letter: counts[letter])
        if counts[letter] > homopolymer_threshold:
            violations.append((i, window, letter, counts[letter]))
    return violations


class HostTableTestCase(unittest.TestCase):
    """Base class that loads the host codon usage table once."""

    @classmethod
    def setUpClass(cls):
        cls.codon_use_table, cls.host_profile, cls.cra = codon_use.host_codon_usage(
            HOST_ID, HOST_THRESHOLD, HOST_TABLE
        )
        cls.dropped_codons = sorted(
            codon for codon, frequency in cls.host_profile.items() if not frequency
        )


class TestRemoveRepeatingSequences(HostTableTestCase):
    """Tests for `codon_harmony.util.seq_opt.remove_repeating_sequences`"""

    # 14 distinct codons: no codon, codon pair or codon triplet is repeated
    non_repeating_dna = "ATGACCGAGTCTAAAGGCTTCCATTGGGTTCCAAACGATCTT"

    def test_sequence_without_repeats_is_unchanged(self):
        """A sequence that contains no repeats must not be mutated.

        The sliding window shrinks to a single codon at the end of the
        sequence, and a one-codon window trivially contains only identical
        codons. That must not be reported as a repeat.
        """
        dna = Seq(self.non_repeating_dna, IUPAC.unambiguous_dna)
        for seed in SEEDS:
            random.seed(seed)
            proc_dna = seq_opt.remove_repeating_sequences(dna, self.codon_use_table, 9)
            assert (
                str(proc_dna) == self.non_repeating_dna
            ), "seed {}: mutated a sequence without repeats " "({} -> {})".format(
                seed, self.non_repeating_dna, proc_dna
            )

    def test_repeats_are_still_removed(self):
        """Genuine repeats must still be detected and removed."""
        # three consecutive identical codons
        dna = Seq("TTGTTGTTGACCGAGTCTAAAGGC", IUPAC.unambiguous_dna)
        for seed in SEEDS:
            random.seed(seed)
            proc_dna = seq_opt.remove_repeating_sequences(dna, self.codon_use_table, 9)
            assert proc_dna != dna, "seed {}: repeat not removed".format(seed)
            assert proc_dna.translate() == dna.translate()

        # a repeated 9-mer
        dna = Seq("ACCGAGTCTAAAACCGAGTCTGGC", IUPAC.unambiguous_dna)
        for seed in SEEDS:
            random.seed(seed)
            proc_dna = seq_opt.remove_repeating_sequences(dna, self.codon_use_table, 9)
            assert proc_dna != dna, "seed {}: repeated 9-mer not removed".format(seed)
            assert proc_dna.translate() == dna.translate()


class TestRemoveSpliceSites(HostTableTestCase):
    """Tests for `codon_harmony.util.seq_opt.remove_splice_sites`"""

    # contains a splice donor; every codon is used by the host
    splice_site_dna = "GCCTCTTGTACGCGGGAAGTGCGGCCATCTGGAAACCTGCGT"

    def test_host_table_drops_codons(self):
        """The test only makes sense if the host threshold dropped codons."""
        assert self.dropped_codons
        assert "AGA" in self.dropped_codons

    def test_only_codons_used_by_the_host_are_substituted_in(self):
        """Codons that the host threshold set to a frequency of zero must
        never be substituted into the sequence: they drive the codon
        adaptation index of the gene to zero.
        """
        dna = Seq(self.splice_site_dna, IUPAC.unambiguous_dna)
        assert not [
            codon for codon in _codons(dna) if not self.host_profile[codon]
        ], "the input sequence must only contain codons used by the host"

        for seed in SEEDS:
            random.seed(seed)
            proc_dna = seq_opt.remove_splice_sites(dna, self.codon_use_table)
            assert proc_dna.translate() == dna.translate()

            unused = [
                codon for codon in _codons(proc_dna) if not self.host_profile[codon]
            ]
            assert (
                not unused
            ), "seed {}: substituted codons unused by the host " "({}) into {}".format(
                seed, unused, proc_dna
            )

    def test_splice_sites_are_still_removed(self):
        """Splice sites must still be removed."""
        dna = Seq(self.splice_site_dna, IUPAC.unambiguous_dna)
        for seed in SEEDS:
            random.seed(seed)
            proc_dna = seq_opt.remove_splice_sites(dna, self.codon_use_table)
            assert proc_dna != dna, "seed {}: splice site not removed".format(seed)

    def test_site_ending_in_the_last_codon_keeps_length_and_protein(self):
        """A splice site that ends in the last codon must not lengthen the
        gene or change the protein.

        The site bounds include the codon one past the site. For a site that
        ends in the last codon that codon does not exist, and it must not be
        written to the sequence.
        """
        # the donor (GGGTGAGT) ends at the last nucleotide
        dna = Seq("ACCGAACGGGTGAGT", IUPAC.unambiguous_dna)
        for seed in SEEDS:
            random.seed(seed)
            proc_dna = seq_opt.remove_splice_sites(dna, self.codon_use_table)
            assert len(proc_dna) == len(dna), "seed {}: {} -> {}".format(
                seed, dna, proc_dna
            )
            assert proc_dna.translate() == dna.translate(), "seed {}: {}".format(
                seed, proc_dna
            )


class TestRemoveLocalHomopolymers(HostTableTestCase):
    """Tests for `codon_harmony.util.seq_opt.remove_local_homopolymers`"""

    def test_all_homopolymers_are_removed(self):
        """Every window must pass, not just the last one.

        The loop has to keep going as long as any window in the pass was
        mutated: a sequence of homopolymers cannot be cleaned up in a
        single pass.
        """
        # Note: not every homopolymer can be removed. Both codons of a
        # window are always mutated together, so e.g. an odd number of Lys
        # codons (AAA/AAG) cannot reach a state in which every window holds
        # four or fewer A's.
        for name, sequence in [("poly-Phe", "TTT" * 4), ("poly-Gly", "GGG" * 4)]:
            dna = Seq(sequence, IUPAC.unambiguous_dna)
            assert _homopolymer_violations(dna), "{}: no homopolymer to remove".format(
                name
            )

            for seed in SEEDS:
                random.seed(seed)
                proc_dna = seq_opt.remove_local_homopolymers(dna, self.codon_use_table)
                assert proc_dna.translate() == dna.translate()

                violations = _homopolymer_violations(proc_dna)
                assert (
                    not violations
                ), "{}, seed {}: homopolymers remain in {} ({})".format(
                    name, seed, proc_dna, violations
                )

    def test_unfixable_window_does_not_keep_looping(self):
        """A window that cannot be fixed must not keep the loop busy.

        `mutate_codon` always returns a different codon, so TTC|TTT and
        TTT|TTC (five T's each) swap places forever. The pass cap has to be
        small and independent of the gene length, otherwise every long gene
        with an `FF` pair pays for thousands of pointless passes.
        """
        dna = Seq("TTCTTT" + "ACCGAGTCT" * 20, IUPAC.unambiguous_dna)
        assert _homopolymer_violations(dna)

        with mock.patch.object(
            seq_opt, "mutate_codon", wraps=seq_opt.mutate_codon
        ) as mutate_codon:
            proc_dna = seq_opt.remove_local_homopolymers(dna, self.codon_use_table)

        assert proc_dna.translate() == dna.translate()
        assert mutate_codon.call_count <= 100, "{} mutations".format(
            mutate_codon.call_count
        )

    def test_sequence_without_homopolymers_is_unchanged(self):
        """A sequence without homopolymers must not be mutated."""
        dna = Seq("ACCGAGTCTACCACCGAGTCTACC", IUPAC.unambiguous_dna)
        assert not _homopolymer_violations(dna)
        for seed in SEEDS:
            random.seed(seed)
            proc_dna = seq_opt.remove_local_homopolymers(dna, self.codon_use_table)
            assert proc_dna == dna, "seed {}: mutated a clean sequence".format(seed)


class TestHarmonizeSequenceAfterSpliceRemoval(HostTableTestCase):
    """Tests for `codon_harmony._harmonize_sequence`

    Splice site removal runs after the best sequence has been picked by its
    codon adaptation index (CAI) and after the deviation from the host
    profile has been checked, so both have to be revisited afterwards.
    """

    # all-His, so that the substituted sequence below is synonymous
    protein = "H" * 40
    # the rarer of the two His codons used by the host, which both lowers
    # the CAI and pushes the codon use away from the host profile
    substituted_dna = "CAC" * len(protein)

    def setUp(self):
        self.args = codon_harmony.get_parser().parse_args(
            [
                "--input",
                "unused.fasta",
                "--output",
                "unused.fasta",
                "--cycles",
                "3",
                "--inner-cycles",
                "1",
                "--no-remove-start-sites",
                "--remove-splice-sites",
                "--restriction-enzymes",
            ]
        )
        self.seq_record = SeqRecord(
            Seq(self.protein, IUPAC.protein),
            id="test_sequence",
            name="Test",
            description="all-His test sequence",
        )

    def _harmonize(self):
        """Run `_harmonize_sequence` with a stubbed splice site removal step
        and return the resulting DNA sequence and the emitted log records.
        """
        from codon_harmony.data import RestrictionEnzymes

        substituted = Seq(self.substituted_dna, IUPAC.unambiguous_dna)
        random.seed(0)
        with mock.patch.object(
            seq_opt, "remove_splice_sites", return_value=substituted
        ) as remove_splice_sites:
            with self.assertLogs(
                "codon_harmony.codon_harmony", level=codon_harmony.logging.OUTPUT
            ) as logs:
                result = codon_harmony._harmonize_sequence(
                    self.seq_record,
                    self.args,
                    self.codon_use_table,
                    self.host_profile,
                    self.cra,
                    RestrictionEnzymes([]),
                )
        assert remove_splice_sites.called
        assert "dna" in result, "no DNA sequence was produced"

        # the sequence is wrapped across several lines of the FASTA record
        dna = Seq("".join(result["dna"].strip().split("\n")[1:]), IUPAC.unambiguous_dna)
        assert str(dna) == self.substituted_dna
        return dna, logs.records

    def test_cai_is_recomputed_after_splice_site_removal(self):
        """The reported CAI must describe the sequence that is written out."""
        dna, records = self._harmonize()
        expected_cai = self.cra.cai_for_gene(str(dna))

        reported = [
            record.getMessage()
            for record in records
            if record.getMessage().startswith("The designed gene's CAI is:")
        ]
        assert len(reported) == 1, "the CAI was not reported exactly once"
        assert "{:.2f}".format(expected_cai) in reported[0], (
            "reported CAI does not describe the final sequence "
            "(expected {:.2f}): {}".format(expected_cai, reported[0])
        )

    def test_warns_when_final_deviation_exceeds_max_relax(self):
        """Splice site removal can push the sequence past `max_relax`, which
        is not re-checked. The user has to be told about it.
        """
        dna, records = self._harmonize()
        _, difference = seq_opt.compare_profiles(
            codon_use.count_codons(dna), self.host_profile, 1 + self.args.max_relax
        )
        assert difference >= self.args.max_relax, (
            "the stubbed sequence is supposed to deviate from the host "
            "profile by more than max_relax"
        )

        warnings = [
            record.getMessage() for record in records if record.levelname == "WARNING"
        ]
        assert any(
            "{:.2f}".format(difference) in warning for warning in warnings
        ), "no warning about the final deviation from the host profile: {}".format(
            warnings
        )


if __name__ == "__main__":
    unittest.main()
