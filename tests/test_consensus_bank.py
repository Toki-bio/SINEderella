"""Tests for RC merge and consensus rebuild."""
import tempfile
import unittest
from pathlib import Path

from consensus_bank_lib import (
    find_rc_clusters,
    identity,
    majority_consensus,
    orient_at_rich_3prime,
    rc,
    read_fa,
    write_fa,
)


class TestRCMerge(unittest.TestCase):
    def test_rc_pair_detected(self):
        a = "GGCCGGATGGCCGAGTGGTAACGCGTTGGCGTGCCACGCAGGAGGACCCG"
        b = rc(a)
        cons = {"oma_A": a, "oma_B": b}
        clusters = find_rc_clusters(cons, 90.0)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(set(clusters[0]), {"oma_A", "oma_B"})

    def test_orient_at_rich_3prime(self):
        # poly-A tail at 3' end in forward orientation
        fwd = "A" * 20 + "G" * 20 + "T" * 30
        rev = rc(fwd)
        self.assertEqual(orient_at_rich_3prime(rev), fwd)

    def test_majority_no_n_on_clear_majority(self):
        seqs = ["ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTTCGT"]
        c = majority_consensus(seqs, tie_to_n=False)
        self.assertNotIn("N", c)
        self.assertEqual(c[4], "A")


class TestCanonicalizeScript(unittest.TestCase):
    def test_merge_writes_output(self):
        import subprocess
        import sys

        repo = Path(__file__).resolve().parents[1]
        a = "GGCCGGATGGCCGAGTGGTAACGCGTTGGCGTGCCACGCAGGAGGACCCG" + "T" * 40
        b = rc(a)
        with tempfile.TemporaryDirectory() as td:
            inp = Path(td) / "in.fa"
            out = Path(td) / "out.fa"
            write_fa(inp, {"oma_keep": a, "oma_drop": b})
            subprocess.check_call(
                [sys.executable, str(repo / "canonicalize_consensus_bank.py"),
                 str(inp), "-o", str(out), "--min-id", "90"],
            )
            merged = read_fa(out)
            self.assertEqual(len(merged), 1)
            self.assertIn("oma_drop", merged)  # lexicographically first name kept


if __name__ == "__main__":
    unittest.main()
