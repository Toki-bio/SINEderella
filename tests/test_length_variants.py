"""Toy tests for tools/length_variants.py: every verdict branch on synthetic copies, plus the BTOP parser."""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import length_variants as lv  # noqa: E402

R = random.Random(11)
BASES = "ACGT"


def rnd(n):
    return "".join(R.choice(BASES) for _ in range(n))


CONS = rnd(200)


def make_copy(end, short_allele, tsd, qs=1):
    mm = {}
    if short_allele:                       # diagnostic columns 25, 26: the short version carries other bases
        for col in (25, 26):
            mm[col] = "G" if CONS[col - 1] != "G" else "T"
    t = rnd(10)
    if tsd:
        up = rnd(20) + t
        dn = "A" * R.randint(5, 30) + t + rnd(60)
    else:
        up = rnd(30)
        dn = "A" * R.randint(5, 30) + rnd(70)
    return lv.Copy(qs, end, mm, up, dn)


def population(n_short, n_long, link=True, tsd_short=0.7, tsd_long=0.7, spread=None):
    cs = []
    for _ in range(n_short):
        short = True if link else R.random() < 0.5
        cs.append(make_copy(R.randint(133, 137), short and R.random() < 0.9, R.random() < tsd_short))
    for _ in range(n_long):
        short = False if link else R.random() < 0.5
        cs.append(make_copy(R.randint(199, 200), short and R.random() < 0.9, R.random() < tsd_long))
    for e in (spread or []):
        cs.append(make_copy(e, R.random() < 0.5, R.random() < 0.3))
    return cs


class TestBtop(unittest.TestCase):
    def test_columns(self):
        mm, last = lv.parse_btop("34CT1GA8", 1)
        self.assertEqual(mm, {35: "T", 37: "A"})
        self.assertEqual(last, 45)

    def test_gaps(self):
        # insertion in the copy (query '-') uses no query column; deletion in the copy (subject '-') does
        mm, last = lv.parse_btop("5-A3A-2", 10)
        self.assertEqual(mm, {18: "-"})
        self.assertEqual(last, 20)


class TestVerdicts(unittest.TestCase):
    def test_two_versions(self):
        rep = lv.analyse(population(600, 800), CONS)
        self.assertEqual(rep["verdict"], "TWO_VERSIONS", rep["why"])
        self.assertEqual(len(rep["modes"]), 2)
        self.assertIn(25, rep["tests"]["linkage"]["diag_columns"])
        self.assertGreater(rep["tests"]["linkage"]["index"], 0.7)

    def test_mixed_short_mode_still_two_versions(self):
        # the short mode also holds decayed long copies (35 % carry the long-type bases): the index drops but stays above the limit
        cs = population(600, 800)
        for c in cs:
            if c.qe < 150 and R.random() < 0.35:
                c.mm = {}
        rep = lv.analyse(cs, CONS)
        self.assertEqual(rep["verdict"], "TWO_VERSIONS", rep["why"])
        self.assertLess(rep["tests"]["linkage"]["index"], 0.9)

    def test_broad_short_mode_is_one_mode(self):
        # the short version's end scatters over 93-115 (two clusters joined by a bridge): one broad mode, still two versions
        cs = population(0, 800)
        for _ in range(700):
            e = R.choice([93, 94, 113, 114]) if R.random() < 0.7 else R.randint(95, 112)
            cs.append(make_copy(e, R.random() < 0.9, R.random() < 0.7))
        rep = lv.analyse(cs, CONS)
        self.assertEqual(rep["verdict"], "TWO_VERSIONS", rep["why"])
        self.assertEqual(len(rep["modes"]), 2)

    def test_tied_alleles_skip_column(self):
        # a column where the short mode gains a base that is also the long mode's majority cannot be diagnostic
        cs = population(600, 800)
        rep = lv.analyse(cs, CONS)
        for col, (bs, bl) in rep["tests"]["linkage"]["diag_columns"].items():
            self.assertNotEqual(bs, bl)

    def test_decay_is_single_mode(self):
        # one length (200) and ends scattered evenly over 110-195: no second mode
        spread = [R.randint(110, 195) for _ in range(600)]
        rep = lv.analyse(population(0, 800, spread=spread), CONS)
        self.assertEqual(rep["verdict"], "SINGLE_MODE", rep["why"])

    def test_unlinked_ends(self):
        # two end modes, but the internal bases do not follow the length
        rep = lv.analyse(population(600, 800, link=False), CONS)
        self.assertEqual(rep["verdict"], "UNLINKED_ENDS", rep["why"])

    def test_no_tsd_in_short_mode_is_unresolved(self):
        rep = lv.analyse(population(600, 800, tsd_short=0.0), CONS)
        self.assertEqual(rep["verdict"], "UNRESOLVED")
        self.assertTrue(any("TSD" in w for w in rep["why"]))

    def test_filled_valley_is_unresolved(self):
        # two modes joined by a dense bridge of ends: not separate versions
        spread = [R.randint(140, 196) for _ in range(2500)]
        rep = lv.analyse(population(600, 800, spread=spread), CONS)
        self.assertNotEqual(rep["verdict"], "TWO_VERSIONS")

    def test_too_few_copies(self):
        rep = lv.analyse(population(20, 30), CONS)
        self.assertEqual(rep["verdict"], "UNRESOLVED")


if __name__ == "__main__":
    unittest.main()
