"""The auditor must fail on falsified claims.

`scripts/validate.py` exists to check the Rust engine's output independently. A review
found it passing output whose every reported metric had been replaced with a flattering
number, so these tests feed it deliberately wrong claims and require each to be caught.

Only the metric and audit helpers are exercised, which import nothing from Biopython,
so this runs wherever pytest does.
"""
import copy
import importlib.util
import json
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TABLE = os.path.join(HERE, "data", "e_coli_316407.json")


def load_validate():
    path = os.path.join(ROOT, "scripts", "validate.py")
    spec = importlib.util.spec_from_file_location("validate_under_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


validate = load_validate()
WEIGHTS = validate.host_weights(json.load(open(TABLE)), 0.10)[1]

# Two encodings of MAKI. The first uses E. coli's preferred codons, the second its
# less favoured ones, so they differ in CAI and the ranking has something to order.
PREFERRED = "ATGGCGAAAATT"
LESS_USED = "ATGGCAAAAATC"


def honest_report(sequences, band=(0.30, 0.65), dev_limit=0.10, threshold=0.80):
    candidates = []
    for rank, dna in enumerate(sequences, start=1):
        candidates.append({
            "rank": rank,
            "output_id": "record_000001_candidate_%02d" % rank,
            "dna": dna,
            "cai": validate.engine_cai(dna, WEIGHTS),
            "gc": validate.gc_fraction(dna),
            "deviation": 0.0,
            "dev_limit": dev_limit,
            "max_relax": dev_limit,
            "run": 1,
            "composition": 0,
            "restart": rank - 1,
            "seed": 1,
            "moves": 0,
            "substitutions": 0,
        })
    best = max(c["cai"] for c in candidates)
    return {
        "candidates_requested": len(sequences),
        "cai_threshold": threshold,
        "gc_band": list(band),
        "policy": {
            "first": {"max_relax": dev_limit, "restarts": 4, "expansion_budget": 2000},
            "retry": {"max_relax": 0.20, "restarts": 8, "expansion_budget": 4000},
        },
        "records": [{
            "id": "p1",
            "candidates": candidates,
            "candidate_count": len(candidates),
            "low_cai": best <= threshold,
            "attempts": [{
                "run": 1,
                "max_relax": dev_limit,
                "restarts": 4,
                "expansion_budget": 2000,
                "best_cai": best,
            }],
        }],
    }


def audit(report):
    written = {"p1": [c["dna"] for c in report["records"][0]["candidates"]]}
    return validate.audit_report(report, written, WEIGHTS, {})


def test_recomputed_cai_matches_the_engines_convention():
    # Every codon of MAKI preferred except the single-codon Met, so CAI is 1.
    assert validate.engine_cai(PREFERRED, WEIGHTS) == pytest.approx(1.0)
    assert validate.engine_cai(LESS_USED, WEIGHTS) < 1.0
    assert validate.gc_fraction("GGCC") == 1.0
    assert validate.gc_fraction("AATT") == 0.0


def test_honest_output_passes():
    assert audit(honest_report([PREFERRED, LESS_USED])) == []


def test_a_falsified_cai_claim_is_caught():
    report = honest_report([PREFERRED, LESS_USED])
    report["records"][0]["candidates"][1]["cai"] = 0.999
    problems = audit(report)
    assert any("claims CAI" in p for p in problems), problems


def test_a_falsified_gc_claim_is_caught():
    report = honest_report([PREFERRED, LESS_USED])
    report["records"][0]["candidates"][0]["gc"] = 0.01
    problems = audit(report)
    assert any("claims GC" in p for p in problems), problems


def test_every_claim_falsified_at_once_still_fails():
    """The exact shape of the reported bypass: flatter every number."""
    report = honest_report([PREFERRED, LESS_USED])
    for candidate in report["records"][0]["candidates"]:
        candidate["cai"] = 0.999
        candidate["gc"] = 0.01
        candidate["deviation"] = 0.0
    report["records"][0]["attempts"][0]["best_cai"] = 0.999
    problems = audit(report)
    assert problems, "a report with no true number in it must not pass"
    assert any("claims CAI" in p for p in problems), problems


def test_a_sequence_outside_the_requested_gc_band_is_caught():
    # both sequences sit at GC 0.333, below a band asking for 0.40 at least
    report = honest_report([PREFERRED, LESS_USED], band=(0.40, 0.65))
    problems = audit(report)
    assert any("outside the requested" in p for p in problems), problems


def test_candidates_out_of_rank_order_are_caught():
    report = honest_report([LESS_USED, PREFERRED])  # ascending CAI
    problems = audit(report)
    assert any("not in descending CAI" in p for p in problems), problems


def test_a_padded_duplicate_is_caught():
    report = honest_report([PREFERRED, PREFERRED])
    problems = audit(report)
    assert any("repeated sequence" in p for p in problems), problems


def test_a_deviation_above_the_bound_that_applied_is_caught():
    report = honest_report([PREFERRED, LESS_USED])
    written = {"p1": [PREFERRED, LESS_USED]}
    # the recomputed deviation, not the claim, is what the bound is checked against
    problems = validate.audit_report(report, written, WEIGHTS, {PREFERRED: 0.5})
    assert any("deviates" in p for p in problems), problems


def test_a_report_that_disagrees_with_the_fasta_is_caught():
    report = honest_report([PREFERRED, LESS_USED])
    problems = validate.audit_report(report, {"p1": [PREFERRED]}, WEIGHTS, {})
    assert any("disagree about its candidates" in p for p in problems), problems


def test_a_missing_retry_after_a_low_best_cai_is_caught():
    report = honest_report([LESS_USED], threshold=0.999)
    report["records"][0]["low_cai"] = True
    problems = audit(report)
    assert any("did not retry" in p for p in problems), problems


def test_a_wrong_low_cai_flag_is_caught():
    report = honest_report([PREFERRED, LESS_USED])
    report["records"][0]["low_cai"] = True
    problems = audit(report)
    assert any("low_cai" in p for p in problems), problems


def test_the_fasta_description_claim_is_checked_too():
    description = "rank=1/5 cai=0.999000 gc=0.0100 source=p1"
    _, problems = validate.check_candidate(LESS_USED, description, WEIGHTS)
    assert any("CLAIMED CAI" in p for p in problems), problems
    assert any("CLAIMED GC" in p for p in problems), problems
