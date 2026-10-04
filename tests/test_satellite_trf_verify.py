"""Tests for the pure parts of tools/satellite_trf_verify.py: gate, window merging, window extraction, TRF parsing, per-locus reduction.
The end-to-end run (needs trf and ssearch36) is done on a server with a planted toy genome (tests/make_toy_satellite.py)."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import satellite_trf_verify as sv  # noqa: E402


class Gate(unittest.TestCase):
    def test_cluster_of_two_or_more_hits(self):
        hits = [("c1", 1000, 1070), ("c1", 1134, 1204), ("c1", 5000, 5250), ("c2", 100, 170), ("c2", 200, 270), ("c2", 330, 400)]
        w = sv.gate(hits, min_hits=2)
        self.assertEqual([(x[0], x[3]) for x in w], [("c1", 2), ("c2", 3)])        # the lone full hit on c1 is not a window

    def test_gap_limit(self):
        hits = [("c1", 0, 70), ("c1", 600, 670)]
        self.assertEqual(sv.gate(hits, min_hits=2), [])
        self.assertEqual(len(sv.gate(hits, gap=600, min_hits=2)), 1)

    def test_every_hit_opens_a_window_by_default(self):
        self.assertEqual(len(sv.gate([("c1", 1000, 1250), ("c1", 9000, 9250)])), 2)

    def test_merge_overlapping_windows(self):
        w = sv.merge_windows([("c1", 0, 500, 2), ("c1", 450, 900, 3), ("c1", 5000, 5400, 2), ("c2", 0, 300, 2)])
        self.assertEqual([(x[0], x[1], x[2]) for x in w], [("c1", 0, 900), ("c1", 5000, 5400), ("c2", 0, 300)])


class Trf(unittest.TestCase):
    def test_parse_and_map_back(self):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "x.trf")
        open(p, "w").write("@c1:1000-3000\n"
                           "10 600 67 8.9 67 90 4 700 30 20 25 25 1.9 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG SEQ\n"
                           "@c2:0-500\n"
                           "5 300 134 2.2 134 80 6 300 25 25 25 25 2.0 ACGT SEQ\n")
        r = sv.parse_trf(p)
        self.assertEqual(r[0][:7], ("c1", 1010, 1600, 67, 8.9, 90, 700))
        self.assertEqual(r[1][:4], ("c2", 5, 300, 134))

    def test_best_per_locus_drops_the_dimer_record(self):
        recs = [("c1", 100, 700, 67, 9.0, 90, 800, "A" * 67), ("c1", 110, 690, 134, 4.3, 80, 500, "A" * 134), ("c1", 5000, 5500, 66, 7.5, 85, 600, "C" * 66)]
        out = sv.best_per_locus(recs)
        self.assertEqual([(r[0], r[3]) for r in out], [("c1", 67), ("c1", 66)])

    def test_extract_windows(self):
        d = tempfile.mkdtemp()
        g = os.path.join(d, "g.fa")
        open(g, "w").write(">c1 desc\n" + "ACGT" * 100 + "\n>c2\n" + "TTGA" * 100 + "\n")
        out = os.path.join(d, "w.fa")
        n = sv.extract_windows(g, [("c1", 0, 200, 2), ("c2", 100, 400, 2)], out, samtools="/no/such/samtools")   # streaming path
        self.assertEqual(n, 2)
        txt = open(out).read().split(">")[1:]
        self.assertTrue(txt[0].startswith("c1:0-200"))
        self.assertEqual(len(txt[1].split("\n", 1)[1].strip()), 300)

    def test_faidx_regions_equals_streaming(self):
        """the samtools path (audit D10) must give the streaming path's sequences, ends past the contig clamped; skipped without samtools"""
        import shutil
        sam = shutil.which("samtools")
        d = tempfile.mkdtemp()
        g = os.path.join(d, "g.fa")
        s1, s2 = "".join("ACGT"[(i * 7) % 4] for i in range(1000)), "".join("TTGAC"[(i * 3) % 5] for i in range(700))
        open(g, "w").write(">c1 desc\n" + "\n".join(s1[i:i + 60] for i in range(0, 1000, 60)) + "\n>c2\n" + s2 + "\n")
        regions = [("c1", 0, 200), ("c1", 950, 1200), ("c2", 100, 400), ("c1", 10, 20)]
        self.assertIsNone(sv.faidx_regions(g, regions, samtools="/no/such/samtools"))
        if not sam:
            self.skipTest("needs samtools")
        got = sv.faidx_regions(g, regions, samtools=sam)
        want = [s1[0:200], s1[950:1000], s2[100:400], s1[10:20]]
        self.assertEqual(got, want)
        out1, out2 = os.path.join(d, "w1.fa"), os.path.join(d, "w2.fa")
        wins = [(c, s, e, 1) for c, s, e in regions]
        sv.extract_windows(g, wins, out1, samtools="/no/such/samtools")
        sv.extract_windows(g, wins, out2, samtools=sam)
        self.assertEqual(open(out1).read(), open(out2).read())


if __name__ == "__main__":
    unittest.main()
