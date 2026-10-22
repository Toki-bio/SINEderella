"""tools/continuation_check.py on a planted case (2026-10-22, the rsi r10 + r8-module case).

Bank: unit A (150 bp) and composite P = A + M + B (420 bp). Genome: 60 copies of A+B (P without its middle piece M; not in
the bank), 60 copies of A, 60 copies of P, 20 unassigned A+B copies. The A+B copies are assigned to A with coordinates that
stop at A's end, as the vote does. Expected: A's assigned copies continue on the 3' side into P around positions 291-420 in
about half of them (flagged), P's copies do not continue, and the candidate built from the continuing copies is A+B,
reported as not in the bank. Needs ssearch36 and bedtools; the build needs mafft and gawk.
"""
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools" / "continuation_check.py"
need = [t for t in ("ssearch36", "bedtools") if shutil.which(t) is None]
pytestmark = pytest.mark.skipif(bool(need), reason="needs %s" % need)


def make_run(d, ab_in_bank=False, ab_locus_full=False):
    rnd = random.Random(3)
    B4 = "ACGT"
    rs = lambda n: "".join(rnd.choice(B4) for _ in range(n))
    mut = lambda s, r=0.05: "".join(rnd.choice(B4.replace(c, "")) if rnd.random() < r else c for c in s)
    A, M, B = rs(150), rs(140), rs(130)
    P = A + M + B
    g = list(rs(400000))
    rows, pos = [], 2000

    def plant(seq):
        nonlocal pos
        s = pos
        g[s:s + len(seq)] = list(seq)
        pos += len(seq) + 1500
        return s

    for _ in range(60):     # A+B, assigned to A; coordinates stop at A's end, or (ab_locus_full) cover A+B as in rsi r10
        s = plant(mut(A + B))
        rows.append(("ctg1:%d-%d(+)" % (s, s + (280 if ab_locus_full else 150)), "A", "assigned"))
    for _ in range(60):
        s = plant(mut(A))
        rows.append(("ctg1:%d-%d(+)" % (s, s + 150), "A", "assigned"))
    for _ in range(60):
        s = plant(mut(P))
        rows.append(("ctg1:%d-%d(+)" % (s, s + len(P)), "P", "assigned"))
    for _ in range(20):     # A+B left unassigned, best match P
        s = plant(mut(A + B))
        rows.append(("ctg1:%d-%d(+)" % (s, s + 280), "P", "rejected_low_bitscore"))
    seq = "".join(g)
    (d / "results").mkdir(parents=True)
    (d / "genome.clean.fa").write_text(">ctg1\n" + seq + "\n")
    (d / "genome.clean.fa.fai").write_text("ctg1\t%d\t6\t%d\t%d\n" % (len(seq), len(seq), len(seq) + 1))
    (d / "consensuses.clean.fa").write_text(">A\n%s\n>P\n%s\n" % (A, P) + (">AB\n%s\n" % (A + B) if ab_in_bank else ""))
    with open(d / "results" / "assignment_full.tsv", "w") as fh:
        fh.write("Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\tRunnerInfo\n")
        for loc, sf, st in rows:
            fh.write("%s\t%s\t500\t10\t%s\t100\t-\n" % (loc, sf, st))
    return A, B


def read_tsv(p):
    lines = [l.rstrip("\n").split("\t") for l in open(p)]
    return [dict(zip(lines[0], l)) for l in lines[1:]]


def run(d, *extra):
    r = subprocess.run([sys.executable, str(TOOL), str(d), "--threads", "2"] + list(extra), capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    return r.stdout


def test_flags_the_missing_composite_report_only(tmp_path):
    make_run(tmp_path)
    out = run(tmp_path, "--no-build")
    assert "CONTINUATION" in out
    rows = read_tsv(tmp_path / "results" / "continuation" / "continuation.tsv")
    flag = [r for r in rows if r["flag"] == "yes"]
    assert len(flag) == 1, rows
    r = flag[0]
    assert (r["consensus"], r["group"], r["side"], r["target"]) == ("A", "assigned", "3'", "P")
    assert 0.4 <= float(r["share"]) <= 0.6
    assert abs(int(r["target_start"]) - 291) <= 5 and abs(int(r["target_end"]) - 420) <= 5
    assert int(r["junction_bp"]) <= 2
    assert r["inside_locus"].startswith("0 of")            # the locus stopped at A's end: the continuation lies outside
    assert not [x for x in rows if x["consensus"] == "P" and x["group"] == "assigned" and x["flag"] == "yes"]


def test_flags_a_long_locus_labelled_with_the_short_consensus(tmp_path):
    # the rsi r10 case as it really was: the A+B loci cover A+B (280 bp) but carry the label A. Reading past the locus end
    # (first version) found only genomic flank; the continuation is past the end of A's match, inside the locus.
    make_run(tmp_path, ab_locus_full=True)
    run(tmp_path, "--no-build")
    rows = read_tsv(tmp_path / "results" / "continuation" / "continuation.tsv")
    flag = [r for r in rows if r["flag"] == "yes"]
    assert len(flag) == 1, rows
    r = flag[0]
    assert (r["consensus"], r["side"], r["target"]) == ("A", "3'", "P")
    assert abs(int(r["target_start"]) - 291) <= 5 and abs(int(r["target_end"]) - 420) <= 5
    ins, span = r["inside_locus"].split(" of ")
    assert int(ins) >= 0.9 * int(span.split()[0])         # all of it inside the locus


@pytest.mark.skipif(shutil.which("mafft") is None or shutil.which("gawk") is None, reason="needs mafft and gawk")
def test_builds_the_candidate_and_finds_it_missing(tmp_path):
    A, B = make_run(tmp_path)
    run(tmp_path)
    rows = read_tsv(tmp_path / "results" / "continuation" / "continuation.tsv")
    r = [x for x in rows if x["flag"] == "yes"][0]
    assert r["candidate"] and 260 <= int(r["candidate_bp"]) <= 300, r
    assert r["bank_match"].startswith("not in the bank"), r["bank_match"]
    cand = open(tmp_path / "results" / "continuation" / "candidates.fa").read().split("\n")[1]
    sys.path.insert(0, str(ROOT / "tools"))
    from consensus_audit import align
    assert align(A + B, cand)[0] <= 4      # the candidate is A+B (5 % copy mutations average out)


@pytest.mark.skipif(shutil.which("mafft") is None or shutil.which("gawk") is None, reason="needs mafft and gawk")
def test_candidate_already_in_bank_points_at_the_assignment(tmp_path):
    # the bank has A+B, yet its copies are assigned to A (here by construction): the report must say the bank has it
    make_run(tmp_path, ab_in_bank=True)
    run(tmp_path)
    rows = read_tsv(tmp_path / "results" / "continuation" / "continuation.tsv")
    r = [x for x in rows if x["flag"] == "yes" and x["consensus"] == "A"][0]
    assert r["bank_match"].startswith("already in the bank as AB"), r["bank_match"]


def test_report_table_shows_the_flag(tmp_path):
    # the report block reads results/continuation/continuation.tsv; a flagged row must appear with its target and bank answer
    sys.path.insert(0, str(ROOT))
    import report_blocks
    d = tmp_path / "results" / "continuation"
    d.mkdir(parents=True)
    cols = ["consensus", "group", "side", "sampled", "target", "copies", "share", "junction_bp", "target_start", "target_end",
            "identity", "inside_locus", "flag", "candidate", "candidate_bp", "seeds_differ_mm", "bank_match"]
    row = ["r10_19seqs", "assigned", "3'", "400", "r10_r8_P18", "234", "0.585", "1", "270", "415", "84.9", "36 of 118 bp", "yes",
           "c1", "226", "2", "not in the bank (closest r10_r8_P18)"]
    (d / "continuation.tsv").write_text("\t".join(cols) + "\n" + "\t".join(row) + "\n")
    html_out = report_blocks._continuation(tmp_path)
    assert "r10_r8_P18" in html_out and "234 of 400" in html_out and "not in the bank" in html_out
    (d / "continuation.tsv").write_text("\t".join(cols) + "\n")
    assert "No consensus is flagged" in report_blocks._continuation(tmp_path)
