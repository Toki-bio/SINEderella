"""Tests for tools/consensus_audit.py: alignment counts, every verdict branch, copy splitting, and (when mafft + gawk exist)
an end-to-end toy run that must give the same rebuild for the same seed."""
import os
import random
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import consensus_audit as ca  # noqa: E402

R = random.Random(5)


def rnd(n):
    return "".join(R.choice("ACGT") for _ in range(n))


class Align(unittest.TestCase):
    def test_identical(self):
        s = rnd(120)
        self.assertEqual(ca.align(s, s), (0, 0))

    def test_substitutions_and_indel(self):
        s = rnd(100)
        t = s[:30] + ("A" if s[30] != "A" else "C") + s[31:60] + s[63:]      # 1 substitution, 3 bp deleted
        mm, gp = ca.align(s, t)
        self.assertEqual((mm, gp), (1, 3))

    def test_tail_difference_counts_as_gaps(self):
        s = rnd(100)
        self.assertEqual(ca.align(s, s + "TTTTT"), (0, 5))


class Verdict(unittest.TestCase):
    def test_branches(self):
        v = ca.verdict
        self.assertEqual(v(150, [150, 149], [0, 1], 1), "MATCH")
        self.assertEqual(v(150, [100, 101], [0, 0], 0), "SHORTER")
        self.assertEqual(v(150, [180, 178], [0, 0], 0), "LONGER")
        self.assertEqual(v(150, [150, 150], [9, 2], 3), "DIVERGED")
        self.assertEqual(v(150, [150, 150], [1, 1], 9), "UNSTABLE")
        self.assertEqual(v(150, [], [], None), "FAILED")
        self.assertEqual(v(600, [585, 590], [2, 3], 1), "MATCH")          # 5 % tolerance on a long consensus
        self.assertEqual(v(40, [34, 35], [0, 0], 0), "MATCH")             # 8 bp floor on a short one


class Split(unittest.TestCase):
    def test_split_by_label(self):
        d = tempfile.mkdtemp()
        try:
            fa = os.path.join(d, "assigned.fasta")
            with open(fa, "w") as fh:
                fh.write(">c1:1-10(+)|famA|100\nAAAA\nCCCC\n>c1:20-30(+)|famB|90\nGGGG\n>c2:1-5(-)|famA|80\nTTTT\n>c3:1-5(+)|other|70\nACAC\n")
            n = ca.split_assigned(fa, d, {"famA", "famB"})
            self.assertEqual(n, {"famA": 2, "famB": 1})
            self.assertEqual(open(os.path.join(d, "famA.fa")).read().count(">"), 2)
            self.assertIn("AAAA\nCCCC", open(os.path.join(d, "famA.fa")).read())
            self.assertFalse(os.path.exists(os.path.join(d, "other.fa")))
        finally:
            shutil.rmtree(d)


@unittest.skipUnless(shutil.which("mafft") and shutil.which("gawk") and shutil.which("bash"), "needs mafft, gawk, bash")
class EndToEnd(unittest.TestCase):
    def test_rebuild_matches_and_is_reproducible(self):
        d = tempfile.mkdtemp()
        try:
            cons = rnd(150)
            fam = "toy"
            os.makedirs(os.path.join(d, "results"))
            with open(os.path.join(d, "consensuses.clean.fa"), "w") as fh:
                fh.write(">%s\n%s\n" % (fam, cons))
            with open(os.path.join(d, "results", "assigned.fasta"), "w") as fh:
                for i in range(60):
                    c = list(cons)
                    for _ in range(R.randint(0, 6)):
                        c[R.randrange(len(c))] = R.choice("ACGT")
                    fh.write(">chr:%d-%d(+)|%s|100\n%s\n" % (i * 500, i * 500 + 150, fam, "".join(c)))
            outs = []
            for _ in range(2):
                self.assertEqual(ca.main.__module__, "consensus_audit")
                sys.argv = ["consensus_audit.py", d, "--jobs", "2"]
                ca.main()
                outs.append(open(os.path.join(d, "results", "consensus_audit", "rebuilt.fa")).read())
            self.assertEqual(outs[0], outs[1])                              # same seeds -> same rebuild
            rows = [l.split("\t") for l in open(os.path.join(d, "results", "consensus_audit", "summary.tsv")).read().splitlines()]
            r = dict(zip(rows[0], rows[1]))
            self.assertEqual(r["verdict"], "MATCH")
        finally:
            shutil.rmtree(d)


if __name__ == "__main__":
    unittest.main()
