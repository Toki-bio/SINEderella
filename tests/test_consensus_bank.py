"""Tests for RC merge and consensus rebuild."""
import tempfile
import unittest
from pathlib import Path

from consensus_bank_lib import (
    find_rc_clusters,
    identity,
    longest_tandem_at_end,
    majority_consensus,
    orient_at_rich_3prime,
    orient_by_simple_repeat_tail,
    rc,
    read_fa,
    tail_strand_score,
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
        fwd = "A" * 20 + "G" * 20 + "T" * 30
        rev = rc(fwd)
        self.assertEqual(orient_at_rich_3prime(rev), fwd)

    def test_trinucleotide_tail_at_3prime(self):
        core = "G" * 40
        fwd = core + ("ATG" * 12)
        rev = rc(fwd)
        r = orient_by_simple_repeat_tail(rev)
        self.assertEqual(r.action, "flipped")
        self.assertEqual(r.seq, fwd)

    def test_homopolymer_tail_orient(self):
        fwd = "G" * 50 + "A" * 24
        rev = rc(fwd)
        end = longest_tandem_at_end(fwd)
        self.assertGreaterEqual(end["repeats"], 3)
        self.assertEqual(orient_by_simple_repeat_tail(rev).action, "flipped")

    def test_no_tail_keeps_input_orientation(self):
        # rsi peel r1_9seqs: no simple-repeat tail on either strand. Was flipped on
        # a 0.2 score difference; must stay as given.
        r1 = ("GGCCCGGTGGCTCAGGCGGTTGGAGCTCCATGCTCCTAACTCCGAAGGCTGCCGGTTCGATTCCCACATGGGCCAGTG"
              "GGCTCTCAACCACAAGGTTGCCAGTTCGACTCCTGATCCCGCAAGGGATGGTGGGCTGTGCCCCCTGCAACTAACAA")
        r = orient_by_simple_repeat_tail(r1)
        self.assertNotEqual(r.action, "flipped")
        self.assertEqual(r.seq, r1)

    def test_sine10_like_backwards(self):
        # T-rich/simple repeat at 5' when backwards; A-run at 3' when flipped
        backwards = "AATTTTTTTTTTTGTTA" + "G" * 80 + "CCGACTTTAT"
        r = orient_by_simple_repeat_tail(backwards)
        self.assertEqual(r.action, "flipped")
        score_fwd, _ = tail_strand_score(backwards)
        score_rev, _ = tail_strand_score(r.seq)
        self.assertGreater(score_rev, score_fwd)

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


class TestNoShorterVariantMerge(unittest.TestCase):
    """rsi r9 (105 bp) is the 5' part of r7 (154 bp) and MEG-RS (135 bp) of MEG-RL (207 bp): both separate SINEs."""

    R7 = ("GGGCGGCCGGTTAGCTCAGTTGGTTAGAGCGCGGTGCTCTTAACAACAAGGTTGCCGGTTCGATCCCCACATGGGCCACTGTGAGCTGCGCCCTCCACAACTAGATTGAAACAACTACTTGACTTGGAGCTGATGGGTCCTGGAAAAACACACT")
    R9 = ("GGGTGGCCGGTTAGCTCAGTTGGTTAGAGCGTGGTGCTAATAACACCAAGGTTGCCGGTTCGATCCCCGCATGGGCCACTGTGAGCTGCGCCCTCCTTAAAAAAA")
    R4 = "CCGGATGGCTCAGTTGGTTGGAGCGCGTGCTCTCAACCACAAGGTTGCCAGTTCGATTCCTCGACTCCCGCAAGGGATGGTGGGCTGTGCCCCCTGCAACTAGCAACGGCAACTGGACCTGGAGCTGAGCTGCGCCCTCCACAA"
    R2 = "CCGGATGGCTCAGTTGGTTGGAGCGCGGGCTCTCAACCACAAGGTTGCCAGTTCAATTCCTCGACTCCCGCAAGGGATGGTGGGCAGCGCCCCCTGCAACTAAAATTGAACACGGCACCTTGAGCTGAGCTGCCGCTGAGCTCCGG"

    def test_shorter_variant_not_merged(self):
        from consensus_bank_lib import find_rc_clusters
        self.assertEqual(find_rc_clusters({"r7": self.R7, "r9": self.R9}, 80.0), [])

    def test_same_length_direct_pair_at_83_percent_not_merged(self):
        from consensus_bank_lib import find_rc_clusters
        self.assertEqual(find_rc_clusters({"r2": self.R2, "r4": self.R4}, 80.0), [])

    def test_rc_pair_and_exact_duplicate_still_merge(self):
        from consensus_bank_lib import find_rc_clusters
        self.assertEqual(find_rc_clusters({"a": self.R7, "b": rc(self.R7)}, 80.0), [["a", "b"]])
        self.assertEqual(find_rc_clusters({"a": self.R7, "b": self.R7}, 80.0), [["a", "b"]])

    def test_length_variant_candidates_listed_not_merged(self):
        from consensus_bank_lib import find_length_variant_pairs
        pairs = find_length_variant_pairs({"r7": self.R7, "r9": self.R9, "r2": self.R2}, min_id=0.85)
        self.assertEqual([(p["short"], p["long"]) for p in pairs], [("r9", "r7")])
        self.assertEqual(pairs[0]["extra"], 49)
        self.assertEqual(pairs[0]["offset"], 0)

    def test_equal_length_pair_is_not_a_length_variant(self):
        from consensus_bank_lib import find_length_variant_pairs
        self.assertEqual(find_length_variant_pairs({"r2": self.R2, "r4": self.R4}, min_id=0.8), [])


if __name__ == "__main__":
    unittest.main()
