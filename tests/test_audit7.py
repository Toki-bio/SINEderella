"""Regression tests for the audit fixes of 2026-11-08 (docs/AUDIT_FIXES_7.md, findings in docs/AUDIT_2026-11-08.md).

Each test names its finding id. The scripts under test are taken from AUDIT7_ROOT when that is set (a copy of the tree with
the files of main put back, to show that the test fails on the old code), else from this repository.
External programs (minimap2, seqkit, ssearch36, ...) are replaced by small stubs, so the tests run without them.
"""
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(os.environ.get("AUDIT7_ROOT", str(ROOT)))
BASH = shutil.which("bash")
need_bash = pytest.mark.skipif(BASH is None, reason="needs bash")


def _stub(bindir, name, body):
    """A small bash program named `name` in bindir."""
    bindir.mkdir(parents=True, exist_ok=True)
    p = bindir / name
    p.write_text("#!/usr/bin/env bash\n" + body + "\n", encoding="utf-8", newline="\n")
    p.chmod(p.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return p


def _bash(script, cwd, bindir=None, env_extra=None, args=()):
    env = dict(os.environ)
    if bindir is not None:
        env["PATH"] = str(bindir) + os.pathsep + env["PATH"]
    env.update(env_extra or {})
    return subprocess.run([BASH, Path(script).as_posix(), *args], cwd=str(cwd), env=env,
                          capture_output=True, text=True, encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------------------------------------------------
# 41: sear, minimap2 path (queries >= 1000 bp): a crashed minimap2 read as "no raw alignments", empty BED, exit 0
# ---------------------------------------------------------------------------------------------------------------------
def _sear_minimap2(tmp_path, minimap2_body):
    b = tmp_path / "bin"
    for t in ("bedtools", "samtools"):
        _stub(b, t, "exit 0")
    _stub(b, "seqkit", "cat >/dev/null 2>&1; exit 0")
    _stub(b, "minimap2", minimap2_body)
    w = tmp_path / "w"
    w.mkdir()
    (w / "q.fa").write_text(">q\n" + "ACGT" * 300 + "\n", encoding="utf-8", newline="\n")   # 1200 bp: minimap2 mode
    (w / "gen.fa").write_text(">c1\n" + "A" * 5000 + "\n", encoding="utf-8", newline="\n")
    (w / "gen.fa.fai").write_text("c1\t5000\t4\t5000\t5001\n", encoding="utf-8", newline="\n")
    return _bash(SRC / "sear", w, b, args=("q.fa", "gen.fa", "0.8", "65", "50")), w


@need_bash
def test_41_sear_minimap2_crash_stops(tmp_path):
    r, w = _sear_minimap2(tmp_path, "echo 'minimap2: index failed' >&2; exit 1")
    assert r.returncode != 0, "a crashed minimap2 must not read as zero copies:\n" + r.stdout[-500:]


@need_bash
def test_41_sear_minimap2_no_hits_still_ok(tmp_path):
    r, w = _sear_minimap2(tmp_path, "cat >/dev/null; exit 0")
    assert r.returncode == 0, r.stdout[-500:] + r.stderr[-500:]
    assert (w / "gen-q.bed").exists() and (w / "gen-q.bed").stat().st_size == 0


# ---------------------------------------------------------------------------------------------------------------------
# 72: step2 --add, tier-1 id -> BED: the contig was cut at the FIRST ':' ("ctg1:apple:100-200(+)" -> ctg1 0 200)
# ---------------------------------------------------------------------------------------------------------------------
def _awk_block(script, end_marker):
    """The awk program of step2's tier-1 id -> BED conversion, cut from the script as it is."""
    lines = Path(script).read_text(encoding="utf-8").splitlines()
    end = next(i for i, l in enumerate(lines) if end_marker in l)
    start = max(i for i in range(end) if lines[i].strip() == "awk '{")
    prog = "\n".join(lines[start:end + 1])
    prog = prog[prog.index("'") + 1:]
    return prog[:prog.index("}' ") + 1]


@pytest.mark.skipif(shutil.which("gawk") is None and shutil.which("awk") is None, reason="needs awk")
def test_72_step2_tier1_id_with_colon_in_contig(tmp_path):
    prog = _awk_block(SRC / "step2_asSINEment.sh", '"$TMPDIR/tier1_present.ids" | sort -k1,1 -k2,2n')
    ids = tmp_path / "ids"
    ids.write_text("ctg1:apple:100-200(+)\nchr2:5-50(-)\n", encoding="utf-8", newline="\n")
    r = subprocess.run([shutil.which("gawk") or shutil.which("awk"), prog, str(ids)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = [l.split("\t") for l in r.stdout.splitlines()]
    assert rows[0][:3] == ["ctg1:apple", "99", "200"], rows
    assert rows[1][:3] == ["chr2", "4", "50"], rows


# ---------------------------------------------------------------------------------------------------------------------
# 20: SINEderella canonicalize_consensus_bank(): a crash of canonicalize_consensus_bank.py was a WARNING, the run went on
# ---------------------------------------------------------------------------------------------------------------------
def _bash_function(script, name):
    """The text of one bash function `name(){ ... }` of a script, from its first line to the first line that is '}'."""
    lines = Path(script).read_text(encoding="utf-8").splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith(name + "(){") or l.startswith(name + "() {"))
    end = next(i for i in range(start + 1, len(lines)) if lines[i] == "}")
    return "\n".join(lines[start:end + 1]) + "\n"


@need_bash
def test_20_canonicalizer_crash_stops_run(tmp_path):
    b = tmp_path / "bin"
    _stub(b, "python3", "echo 'Traceback: boom' >&2; exit 1")
    sd = tmp_path / "sd"
    sd.mkdir()
    (sd / "canonicalize_consensus_bank.py").write_text("", encoding="utf-8")
    (tmp_path / "bank.fa").write_text(">a\nACGT\n", encoding="utf-8", newline="\n")
    drv = tmp_path / "drv.sh"
    drv.write_text("set -euo pipefail\nlog(){ echo \"$*\" >&2; }\ndie(){ echo \"ERROR: $*\" >&2; exit 1; }\n"
                   "SCRIPT_DIR=%s\n" % sd.as_posix()
                   + _bash_function(SRC / "SINEderella", "canonicalize_consensus_bank")
                   + "canonicalize_consensus_bank %s\necho CONTINUED\n" % (tmp_path / "bank.fa").as_posix(),
                   encoding="utf-8", newline="\n")
    r = _bash(drv, tmp_path, b)
    assert "CONTINUED" not in r.stdout and r.returncode != 0, r.stderr


# ---------------------------------------------------------------------------------------------------------------------
# 53: tools/rollup_families.py --votes: a family with no subfamily threshold got threshold 0, every 10-vote copy assigned
# ---------------------------------------------------------------------------------------------------------------------
def _w(p, text):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8", newline="\n")


def test_53_rollup_family_without_threshold_not_assigned(tmp_path):
    sys.path.insert(0, str(SRC / "tools"))
    import importlib
    import rollup_families as rf
    rf = importlib.reload(rf)
    assert Path(rf.__file__).resolve().parent == (SRC / "tools").resolve()
    run = tmp_path / "run"
    _w(run / "step2/step2_output/assignment_full.tsv",
       "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\tRunnerInfo\n"
       "a1\tA1\t300\t10\tassigned\t100\tNA\n"            # family A has a threshold (100)
       "a2\tA2\t200\t6\tno_unanimous\t100\tNA\n"         # A1/A2 split, score 1900 >= 100 -> family-firm
       "c1\tC1\t40\t5\tno_unanimous\tNA\tNA\n"           # family C: no subfamily ever got a threshold
       )
    _w(run / "genome.clean_step1/extracted.fasta", ">a1\nACGT\n")
    _w(run / "genome.clean_step1/searches/all_hits.labeled.bed", "chr1\t0\t100\tA1\t1\t+\n")
    _w(run / "fam.tsv", "subfamily\tfamily\tparent_source\nA1\tA\tderived\nA2\tA\tderived\nC1\tC\tderived\nC2\tC\tderived\n")
    _w(run / "fam.fa", ">A\nACGT\n>C\nTTTT\n")
    _w(run / "votes.tsv", "a2\tA1\t4\t800\t0\tNA\na2\tA2\t6\t1100\t0\tNA\nc1\tC1\t5\t40\t0\tNA\nc1\tC2\t5\t38\t0\tNA\n")
    out = tmp_path / "fv"
    rf.rollup(str(run), str(run / "fam.tsv"), str(run / "fam.fa"), str(out), str(run / "votes.tsv"))
    rows = [l.split("\t") for l in (out / "step2/step2_output/assignment_full.tsv").read_text(encoding="utf-8").splitlines()]
    t = {r[0]: r for r in rows[1:]}
    assert t["a2"][4] == "assigned"                               # the tested family still promotes
    assert t["c1"][4] != "assigned", t["c1"]                      # untested: not a family-firm copy
    assert t["c1"][1] == "C"


# ---------------------------------------------------------------------------------------------------------------------
# 9: SINEderella --resume of a run whose step 2+3 is done never ran families_step (family roll-up missing or stale)
# ---------------------------------------------------------------------------------------------------------------------
@need_bash
def test_9_resume_runs_families_step(tmp_path):
    run = tmp_path / "run"
    for f in ("genome.clean.fa", "consensuses.clean.fa", "genome.clean_step1/extracted.fasta",
              "genome.clean_step1/.step1_complete", "step2/step2_output/summary.by_subfam.tsv"):
        _w(run / f, ">x\nACGT\n")
    sd = tmp_path / "sd"
    for f in ("s1", "s2", "s3"):
        _w(sd / f, "")
    calls = tmp_path / "calls.txt"
    stubs = "".join("%s(){ echo %s >> %s; }\n" % (n, n, calls.as_posix()) for n in (
        "families_step", "build_results_links", "array_flag_step", "flank_twins_step", "verdict_step",
        "length_variants_step", "consensus_audit_step", "continuation_step", "report_step", "run_step2_step3"))
    drv = tmp_path / "drv.sh"
    drv.write_text("set -euo pipefail\nlog(){ echo \"$*\" >&2; }\ndie(){ echo \"ERROR: $*\" >&2; exit 1; }\n" + stubs
                   + "SCRIPT_DIR=%s\nSTEP1_SRC=%s/s1\nASSIGN_SRC=%s/s2\nPOST_SRC=%s/s3\nTHREADS=1\nEXISTING_RUN=%s\n"
                   % ((sd.as_posix(),) * 4 + (run.as_posix(),))
                   + _bash_function(SRC / "SINEderella", "mode_resume") + "mode_resume\n",
                   encoding="utf-8", newline="\n")
    r = _bash(drv, tmp_path)
    assert r.returncode == 0, r.stderr[-1500:]
    done = calls.read_text().split()
    assert "report_step" in done and "run_step2_step3" not in done          # the "step 2+3 already complete" branch
    assert "families_step" in done, done
    assert done.index("families_step") < done.index("report_step")


# ---------------------------------------------------------------------------------------------------------------------
# 16: tools/satellite_trf_verify.py: a kind-A ssearch36 that failed twice (or exited 0 with '***' on stderr) gave an empty
# result, read by satellite_stage.py as "no SINE-derived tandem loci"
# ---------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("rc,err", [(1, "error"), (0, " *** cannot open library units2.fa\n")])
def test_16_trf_verify_failed_search_stops(tmp_path, monkeypatch, rc, err):
    sys.path.insert(0, str(SRC / "tools"))
    import importlib
    import satellite_trf_verify as V
    V = importlib.reload(V)
    assert Path(V.__file__).resolve().parent == (SRC / "tools").resolve()
    (tmp_path / "in.fa").write_text(">c\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(V.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a[0], rc, "", err))
    recs = [("c1", 0, 300, 100, 3.0, 90, 200, "ACGT" * 25)]
    with pytest.raises(SystemExit):
        V.verify(recs, str(tmp_path / "in.fa"), "ssearch36", str(tmp_path), 1)


# ---------------------------------------------------------------------------------------------------------------------
# 50: tools/kindB_flank_check.sh: the piece-median loop was bounded by the count of flank-eligible hits, not by the stretches
# ---------------------------------------------------------------------------------------------------------------------
def _shell_block(script, first_marker, last_marker):
    """Lines of a shell script from the line before `first_marker` (the command start) to the line holding `last_marker`."""
    lines = Path(script).read_text(encoding="utf-8").splitlines()
    i = next(k for k, l in enumerate(lines) if first_marker in l)
    j = next(k for k in range(i, len(lines)) if last_marker in lines[k])
    return "\n".join(lines[i - 1:j + 1]) + "\n"


@need_bash
@pytest.mark.skipif(shutil.which("gawk") is None, reason="needs gawk")
def test_50_kindB_piece_median_reads_every_stretch(tmp_path):
    # Fable's toy t50: run r, 4 hits (stretches 0..2), 3 flank-eligible hits, the only piece at stretch 2 scores 95
    for f, t in (("border.ids", "r\n"), ("pnames.txt", "r_2_P\n"), ("piece_scores.tsv", "r\t2\t95\n"),
                 ("verdicts.flank.tsv", "r\t3\t60.0\t0.33\tBORDER\tflank\n"), ("nhits.tsv", "r\t4\n")):
        _w(tmp_path / f, t)
    blk = _shell_block(SRC / "tools/kindB_flank_check.sh", 'FILENAME == "pnames.txt" { split($0, a, "_"); if (a[3] == "P")',
                       "> piece_median.tsv")
    (tmp_path / "blk.sh").write_text("set -euo pipefail\n" + blk, encoding="utf-8", newline="\n")
    r = _bash(tmp_path / "blk.sh", tmp_path)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "piece_median.tsv").read_text().split() == ["r", "95", "1"]


# ---------------------------------------------------------------------------------------------------------------------
# D16: flankscan/fs9_run.sh: a family that fails in a rerun inside the same run dir kept its old copy_status.tsv
# ---------------------------------------------------------------------------------------------------------------------
@need_bash
@pytest.mark.skipif(shutil.which("gawk") is None, reason="needs gawk")
def test_D16_fs9_rerun_failed_family_leaves_no_stale_status(tmp_path):
    fs = tmp_path / "flankscan"
    fs.mkdir()
    shutil.copy(SRC / "flankscan/fs9_run.sh", fs / "fs9_run.sh")
    _w(fs / "fs9_twins.sh", "echo 'fs9_twins: boom' >&2; exit 1\n")          # this run: the family fails
    run = tmp_path / "run"
    _w(run / "genome.clean.fa", ">c1\n" + "A" * 1000 + "\n")
    _w(run / "genome.clean.fa.fai", "c1\t1000\t4\t1000\t1001\n")
    _w(run / "results/assignment_full.tsv", "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\tRunnerInfo\n"
       "c1:100-200(+)\tF\t300\t10\tassigned\t100\tNA\nc1:500-600(-)\tF\t300\t10\tassigned\t100\tNA\n")
    stale = run / "results/flank_twins/F/copy_status.tsv"                      # the earlier run's result
    _w(stale, "copy\tstatus\nc1:100-200(+)\tunique\nc1:500-600(-)\tunique\n")
    r = _bash(fs / "fs9_run.sh", tmp_path, args=(run.as_posix(), "1"))
    assert r.returncode == 0, r.stderr
    assert "F failed" in r.stderr
    assert not stale.exists(), "the failed family's old copy_status.tsv is still read as this run's result"


# ---------------------------------------------------------------------------------------------------------------------
# 224: tools/length_variants.py: with a mode under 50 copies no linkage column is tested, yet index 0.0 gave UNLINKED_ENDS
# ---------------------------------------------------------------------------------------------------------------------
def _lv():
    sys.path.insert(0, str(SRC / "tools"))
    import importlib
    import length_variants as lv
    lv = importlib.reload(lv)
    assert Path(lv.__file__).resolve().parent == (SRC / "tools").resolve()
    return lv


def _lv_copies(lv, n_short, n_long):
    import random
    rng = random.Random(1)
    rseq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    cons = rseq(200)
    # ends spread over 7 bp around 100 and 190: two modes with an empty valley between them
    cs = [lv.Copy(1, 97 + i % 7, {}, rseq(30), rseq(120), "S") for i in range(n_short)]
    cs += [lv.Copy(1, 187 + i % 7, {}, rseq(30), rseq(120), "L") for i in range(n_long)]
    return cs, cons


def test_224_length_variants_no_tested_column_is_untested():
    # as Fable's toy: 40 short-mode copies, 160 long, no base differences -> every column skipped (short mode < 50)
    lv = _lv()
    cs, cons = _lv_copies(lv, 40, 160)
    rep = lv.analyse(cs, cons)
    assert rep["verdict"] != "UNLINKED_ENDS", rep["why"]
    assert rep["verdict"] == "UNTESTED", rep["verdict"]
    assert rep["tests"]["linkage"]["index"] is None


def test_224_length_variants_measured_columns_keep_unlinked_ends():
    # control: 60 short copies, columns are measured and none differs -> UNLINKED_ENDS stays
    lv = _lv()
    cs, cons = _lv_copies(lv, 60, 160)
    rep = lv.analyse(cs, cons)
    assert rep["verdict"] == "UNLINKED_ENDS" and rep["tests"]["linkage"]["index"] == 0.0


def test_224_report_words_untested():
    sys.path.insert(0, str(SRC))
    import importlib
    import report_blocks as rb
    rb = importlib.reload(rb)
    assert "UNTESTED" in rb.WORDS and "not tested" in rb.WORDS["UNTESTED"]


# ---------------------------------------------------------------------------------------------------------------------
# L94: step8a top-100 sort: no LC_ALL and no tie key, so plate membership at equal bitscores depended on the locale
# ---------------------------------------------------------------------------------------------------------------------
@need_bash
def test_L94_step8a_top100_sort_ties_fixed_order(tmp_path):
    lines = (SRC / "step8a_extract_alignments.sh").read_text(encoding="utf-8").splitlines()
    srt = [l for l in lines if '"$TMPDIR/loci_${idx}.tsv" > "$TMPDIR/sorted_${idx}.tsv"' in l or
           '"$TMPDIR/soft_${idx}.tsv" >> "$TMPDIR/sorted_${idx}.tsv"' in l]
    assert len(srt) == 2
    for l in srt:
        assert "LC_ALL=C" in l, l
    # equal bitscores: the explicit key orders them by contig, then start as a number (99 before 100); the whole-line fallback
    # put "100" before "99"
    _w(tmp_path / "loci_1.tsv", "SF\t500\tc2\t100\t300\t+\nSF\t500\tc2\t99\t300\t+\nSF\t900\tc9\t5\t200\t-\n")
    (tmp_path / "s.sh").write_text("TMPDIR=%s\nidx=1\n%s\n" % (tmp_path.as_posix(), srt[0].strip()), encoding="utf-8", newline="\n")
    r = _bash(tmp_path / "s.sh", tmp_path, env_extra={"LC_ALL": "en_US.UTF-8"})
    assert r.returncode == 0, r.stderr
    got = [l.split("\t")[2:4] for l in (tmp_path / "sorted_1.tsv").read_text().splitlines()]
    assert got == [["c9", "5"], ["c2", "99"], ["c2", "100"]], got


# ---------------------------------------------------------------------------------------------------------------------
# 38: publish_run.sh picked the step-2 output by name order (step2_output_run9 over _run10)
# ---------------------------------------------------------------------------------------------------------------------
@need_bash
def test_38_publish_picks_newest_step2_output(tmp_path):
    # Fable's toy t38: six step-2 output dirs, made in run order; the newest (by time, as every other step) is run11
    import os, time
    run = tmp_path / "run"
    for i, d in enumerate(("step2_output", "step2_output_run1", "step2_output_run2", "step2_output_run9",
                           "step2_output_run10", "step2_output_run11")):
        (run / "step2" / d).mkdir(parents=True)
        os.utime(run / "step2" / d, (time.time() - 1000 + 10 * i, time.time() - 1000 + 10 * i))
    blk = _shell_block(SRC / "publish_run.sh", "step2/step2_output* 2>/dev/null", 'die "No step2_output under')
    (tmp_path / "s.sh").write_text("set -euo pipefail\ndie(){ echo \"$*\" >&2; exit 1; }\nRUN_ROOT=%s\n%s\necho \"$S2\"\n"
                                   % (run.as_posix(), blk.split("\n", 1)[1]), encoding="utf-8", newline="\n")
    r = _bash(tmp_path / "s.sh", tmp_path)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().endswith("/step2_output_run11"), r.stdout


# ---------------------------------------------------------------------------------------------------------------------
# D17: publish/align_for_publish.sh: a stale consensus alignment of an earlier publish was copied into the HTML tree again
# ---------------------------------------------------------------------------------------------------------------------
@need_bash
@pytest.mark.skipif(shutil.which("awk") is None, reason="needs awk")
def test_D17_republish_drops_stale_sine_consensus_alignment(tmp_path):
    b = tmp_path / "bin"
    _stub(b, "mafft", 'for a in "$@"; do f="$a"; done; cat "$f"')          # "aligns" by printing the input
    run = tmp_path / "run"
    _w(run / "consensuses.clean.fa", ">A\nACGT\n>B\nACGA\n>C\nTTTT\n")
    # this publish: only one subfamily is interspersed -> no SINE consensus alignment is written
    _w(run / "step2/step2_output/flank_verdict.tsv", "subfamily\tverdict\nA\tinterspersed\nB\tnot_interspersed\nC\tuntested\n")
    s8, od = tmp_path / "s8", tmp_path / "html"
    for d in (s8, od):
        _w(d / "SP_consensuses_sine.aln.fa", ">A\nACGT\n>OLD\nACGT\n")       # left by an earlier publish
    # the consensus-alignment block of the script as it is, from CONS_BANK= to the final log line
    lines = (SRC / "publish/align_for_publish.sh").read_text(encoding="utf-8").splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith('CONS_BANK="$RUN_ROOT/consensuses.clean.fa"'))
    j = next(k for k in range(i, len(lines)) if lines[k].startswith('log "Done: $STEP8_OUT"'))
    blk = "\n".join(lines[i:j]) + "\n"
    (tmp_path / "s.sh").write_text("set -euo pipefail\nlog(){ echo \"$*\" >&2; }\nRUN_ROOT=%s\nSTEP8_OUT=%s\nOUT_DIR=%s\nSPECIES=SP\n%s"
                                   % (run.as_posix(), s8.as_posix(), od.as_posix(), blk), encoding="utf-8", newline="\n")
    r = _bash(tmp_path / "s.sh", tmp_path, b)
    assert r.returncode == 0, r.stderr
    assert (od / "SP_consensuses.aln.fa").exists()
    assert not (od / "SP_consensuses_sine.aln.fa").exists(), "the earlier publish's SINE consensus alignment is published again"
