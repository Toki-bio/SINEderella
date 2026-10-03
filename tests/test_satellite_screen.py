"""Tests for tools/satellite_screen.py: kind A (short monomers < 100 bp apart), kind B (regular spacing of full hits), controls."""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import satellite_screen as ss  # noqa: E402

R = random.Random(9)


def hits_of(spec):
    out = []
    for c, s, e, st in spec:
        out.append((c, s, e, st))
    return out


class Screen(unittest.TestCase):
    def test_kind_a_monomer_run(self):
        # 30 monomers of 120 bp, 20 bp apart: a satellite-like locus
        h = [("c1", 1000 + i * 140, 1000 + i * 140 + 120, "+") for i in range(30)]
        rows, s = ss.screen(h)
        self.assertEqual(s["loci_A"], 1)
        self.assertEqual(s["hits_in_A"], 30)
        self.assertEqual(rows[0][0], "A")
        self.assertEqual(rows[0][4], 30)

    def test_three_monomers_stay(self):
        h = [("c1", 1000 + i * 140, 1000 + i * 140 + 120, "+") for i in range(3)]       # a trimer is not a satellite (>= 4)
        rows, s = ss.screen(h)
        self.assertEqual(s["hits_in_A"], 0)

    def test_strand_flip_breaks_run(self):
        h = [("c1", 1000 + i * 140, 1000 + i * 140 + 120, "+" if i < 6 else "-") for i in range(12)]
        rows, s = ss.screen(h)
        self.assertEqual(s["loci_A"], 2)

    def test_kind_b_regular_spacing(self):
        # full-length hits 250 bp every ~2 100 bp: units much longer than the SINE
        h = [("c1", 5000 + i * 2100 + R.randint(-20, 20), 5250 + i * 2100, "+") for i in range(25)]
        rows, s = ss.screen(h)
        self.assertEqual(s["hits_in_A"], 0)
        self.assertEqual(s["hits_in_B"], 25)
        self.assertEqual(rows[0][0], "B")
        self.assertTrue(2000 <= rows[0][7] <= 2200)

    def test_dispersed_copies_are_clean(self):
        h = []
        for c in range(6):
            p = 0
            for _ in range(400):
                p += int(R.expovariate(1 / 30000.0)) + 300
                h.append(("c%d" % c, p, p + 250, R.choice("+-")))
        rows, s = ss.screen(h)
        self.assertLess((s["hits_in_A"] + s["hits_in_B"]) / len(h), 0.01)

    def test_dense_dispersed_family_chance_is_subtracted(self):
        # one hit per ~12 kb (like gja Squam3A): regular spacing occurs by chance, the null must absorb it
        h = []
        for c in range(4):
            p = 0
            for _ in range(2000):
                p += int(R.expovariate(1 / 12000.0)) + 300
                h.append(("c%d" % c, p, p + 250, "+"))
        rows, s = ss.screen(h)
        self.assertLess(s["excess_B"], 3.0, s)

    def test_real_array_in_dense_background(self):
        h = []
        for c in range(4):
            p = 0
            for _ in range(2000):
                p += int(R.expovariate(1 / 12000.0)) + 300
                h.append(("c%d" % c, p, p + 250, "+"))
        h += [("arr", 5000 + i * 2100 + R.randint(-20, 20), 5250 + i * 2100, "+") for i in range(600)]
        rows, s = ss.screen(h)
        self.assertGreater(s["excess_B"], 5.0, s)

    def test_dimer_and_tandem_pair_are_not_flagged(self):
        h = [("c1", 1000, 1250, "+"), ("c1", 1290, 1540, "+"), ("c1", 50000, 50250, "+"), ("c1", 50050 + 250, 50550, "+")]
        rows, s = ss.screen(h)
        self.assertEqual(s["hits_in_A"] + s["hits_in_B"], 0)

    def test_a_hits_not_double_counted_in_b(self):
        h = [("c1", 1000 + i * 140, 1000 + i * 140 + 120, "+") for i in range(30)]
        rows, s = ss.screen(h)
        self.assertEqual(s["hits_in_B"], 0)


if __name__ == "__main__":
    unittest.main()
