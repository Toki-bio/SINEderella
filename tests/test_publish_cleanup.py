#!/usr/bin/env python3
"""Verify publish cleanup: no Tal defaults, alignment link modes."""
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from urllib.parse import quote

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("step6", REPO / "step6_report.py")
step6 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(step6)


class TestBuildAlignmentSection(unittest.TestCase):
    def setUp(self):
        self.aln = Path(tempfile.mkdtemp()) / "alignments"
        self.aln.mkdir()
        for suf in ("top100", "rand100", "subfam"):
            (self.aln / f"mysp_SINE10_{suf}.aln.fa").write_text(">x\nA\n")

    def test_local_relative_links(self):
        html = step6.build_alignment_section("mysp", ["SINE10"], aln_dir=self.aln)
        self.assertIn("alignments/mysp_SINE10_top100.aln.fa", html)
        self.assertIn("Relative links", html)
        self.assertNotIn("Toki-bio/Tal", html)

    def test_remote_msa_viewer(self):
        raw = "https://raw.githubusercontent.com/org/repo/main/mysp/alignments/"
        html = step6.build_alignment_section(
            "mysp", ["SINE10"], raw_base=raw, aln_dir=self.aln
        )
        fn = "mysp_SINE10_top100.aln.fa"
        self.assertIn("toki-bio.github.io/MSA-viewer", html)
        self.assertIn(quote(raw + fn, safe=""), html)

    def test_remote_repo_root(self):
        raw = "https://raw.githubusercontent.com/org/repo/main/"
        html = step6.build_alignment_section(
            "mysp", ["SINE10"], raw_base=raw, aln_dir=self.aln
        )
        expected = raw.rstrip("/") + "/mysp/alignments/mysp_SINE10_top100.aln.fa"
        self.assertIn(quote(expected, safe=""), html)

    def test_empty_when_no_files(self):
        html = step6.build_alignment_section("mysp", ["SINE99"], aln_dir=self.aln)
        self.assertEqual(html, "")


class TestNoTalDefaults(unittest.TestCase):
    def test_source_grep(self):
        src = (REPO / "step6_report.py").read_text(encoding="utf-8")
        self.assertNotIn("Toki-bio/Tal", src)
        self.assertNotIn("_others", src)
        self.assertNotIn("raw.githubusercontent.com/Toki-bio/Tal", src)
        self.assertIn("pages_index: Optional[str]", src)


if __name__ == "__main__":
    unittest.main()
