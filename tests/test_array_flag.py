"""tools/array_flag.py end to end on a synthetic run: a dense dispersed family (one copy per 3 kb, like tbr VES) forms regular runs
by chance and must NOT be flagged ARRAY; a family whose copies sit in 2 kb arrays must. Found 2026-10-05: VES (620 000 copies) was
flagged "Tandem array" because the flag compared the raw share with 20 % and had no chance null."""
import os
import random
import subprocess
import sys
import tempfile
import unittest

TOOL = os.path.join(os.path.dirname(__file__), "..", "tools", "array_flag.py")


class ArrayFlag(unittest.TestCase):
    def test_dense_dispersed_not_flagged_array_flagged(self):
        R = random.Random(5)
        rows = []
        for c in range(6):                                   # dense dispersed: 6 contigs x 1 500 copies, mean gap 3 kb
            p = 0
            for _ in range(1500):
                p += int(R.expovariate(1 / 3000.0)) + 250
                rows.append(("DENSE", "ctg%d:%d-%d(+)" % (c, p, p + 250)))
        for a in range(20):                                  # arrays: 20 arrays of 25 units, 2 kb spacing, on their own contigs
            for i in range(25):
                s = 100000 + a * 500000 + i * 2000
                rows.append(("ARR", "arr%d:%d-%d(+)" % (a, s, s + 250)))
        for c in range(3):                                   # plus 300 dispersed ARR copies, 100 kb apart on three long contigs
            for i in range(100):
                s = 50000 + i * 100000
                rows.append(("ARR", "far%d:%d-%d(+)" % (c, s, s + 250)))
        d = tempfile.mkdtemp()
        os.makedirs(os.path.join(d, "results"))
        with open(os.path.join(d, "results", "assignment_full.tsv"), "w") as fh:
            fh.write("Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\n")
            for fam, loc in rows:
                fh.write("%s\t%s\t1000\t10\tassigned\t500\n" % (loc, fam))
        r = subprocess.run([sys.executable, TOOL, d], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        lines = [l.split("\t") for l in open(os.path.join(d, "results", "array_flag.tsv")).read().splitlines()]
        h = lines[0]
        tab = {l[0]: dict(zip(h, l)) for l in lines[1:]}
        self.assertIn("null_pct", h)
        dense, arr = tab["DENSE"], tab["ARR"]
        self.assertGreater(float(dense["pct_in_arrays"]), 10.0, dense)      # chance runs are common at this density ...
        self.assertEqual(dense["flag"], "-", dense)                         # ... and the null absorbs them
        self.assertLess(float(dense["excess_pct"]), 20.0, dense)
        self.assertEqual(arr["flag"], "ARRAY", arr)
        self.assertGreater(float(arr["excess_pct"]), 50.0, arr)


if __name__ == "__main__":
    unittest.main()
