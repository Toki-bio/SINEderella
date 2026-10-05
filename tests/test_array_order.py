"""Tests for tools/array_order.py: the regular-spacing rule over all loci, and the reordering of the top 100."""
import os
import random
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import array_order as ao  # noqa: E402

R = random.Random(3)


class Regular(unittest.TestCase):
    def test_array_is_found(self):
        loci = [("c1", 1000 + i * 2050 + R.randint(-30, 30)) for i in range(20)]
        self.assertEqual(len(ao.regular_runs(loci)), 20)

    def test_array_of_7kb_period(self):
        # rsi MEG-RS NC_142507.1:29.8-30.6 Mb: 50 copies ~7.0 kb apart, missed with the 6 kb cap
        loci = [("c1", 29819634 + i * 7020 + R.randint(-40, 40)) for i in range(50)]
        self.assertEqual(len(ao.regular_runs(loci)), 0)
        self.assertEqual(len(ao.regular_runs(loci, ao.WIDE_GAP, ao.WIDE_MIN)), 50)

    def test_array_with_unit_length_variants(self):
        gaps = [2000, 2050, 1400, 2300, 2040, 1450, 2060, 2300]
        pos, loci = 5000, []
        for g in gaps:
            loci.append(("c1", pos)); pos += g
        loci.append(("c1", pos))
        self.assertEqual(len(ao.regular_runs(loci)), len(loci))

    def test_dispersed_copies_are_not_marked(self):
        # 2 000 copies on 4 contigs, exponential gaps (mean 20 kb): chance runs must be rare
        loci = []
        for c in range(4):
            p = 0
            for _ in range(500):
                p += int(R.expovariate(1 / 20000.0)) + 1
                loci.append(("c%d" % c, p))
        frac = len(ao.regular_runs(loci)) / len(loci)
        self.assertLess(frac, 0.03, frac)

    def test_far_apart_copies(self):
        loci = [("c1", i * 40000) for i in range(10)]
        self.assertEqual(ao.regular_runs(loci), {})

    def test_two_arrays_one_contig(self):
        a = [("c1", 1000 + i * 2000) for i in range(8)]
        b = [("c1", 900000 + i * 3000) for i in range(8)]
        r = ao.regular_runs(a + b)
        self.assertEqual(len(r), 16)
        self.assertEqual(len(set(r.values())), 2)


class Order(unittest.TestCase):
    def test_independent_loci_first_even_when_arrays_rank_first(self):
        rows = []
        for i in range(120):                       # the array copies score highest
            rows.append(["fam", str(1000 - i), "c1", str(1000 + i * 2000), str(1135 + i * 2000), "+"])
        for i in range(40):                        # 40 independent copies score lower
            rows.append(["fam", str(500 - i), "c%d" % (2 + i), str(10 ** 6 + i * 777), str(10 ** 6 + i * 777 + 135), "+"])
        d = tempfile.mkdtemp()
        p = os.path.join(d, "loci.tsv")
        with open(p, "w") as fh:
            fh.write("\n".join("\t".join(r) for r in rows) + "\n")
        out = subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "..", "tools", "array_order.py"), p],
                             capture_output=True, text=True, check=True).stdout.splitlines()
        top = [l.split("\t")[2] for l in out[:100]]
        self.assertEqual(sum(1 for c in top if c != "c1"), 40)      # all independent copies are on the plate
        self.assertEqual(sum(1 for c in top if c == "c1"), 60)      # the rest filled from the array
        self.assertTrue(all(l.split("\t")[7] == "array" for l in out if l.split("\t")[2] == "c1"))


if __name__ == "__main__":
    unittest.main()
