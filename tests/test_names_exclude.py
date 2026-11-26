"""Regression tests for three faults seen on the apple run (2026-11-26):

1. report labels: every consensus_<N>seqs was shown as "consensus" (report_blocks.lab, report_hierarchy._short / _full / _plate,
   flankscan/fs7_hierarchy.sh short()); a short label is used only when it is unique in the bank (r6_210seqs -> r6 still).
2. --exclude: the satellite tables copied from the source run kept the excluded consensus (SINEderella satellite_drop_consensus).
3. --exclude / --add: the all-seed SubFam alignment (subfam_input/input.clw.al) was the source run's, with the excluded seed as a
   row; rebuild_from_searches now rebuilds it from this run's loci and bank (SINEderella rebuild_allseeds, consensus.path).

The code under test is taken from NAMES_TEST_ROOT when that is set (a copy of the old tree, to show the tests fail there), else
from this repository. Shell tests need bash and gawk (Git Bash on Windows); SubFam and mafft are replaced by stubs.
"""
import importlib
import os
import random
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(os.environ.get("NAMES_TEST_ROOT", str(ROOT)))
BASH, GAWK = shutil.which("bash"), shutil.which("gawk")
need_shell = pytest.mark.skipif(BASH is None or GAWK is None, reason="needs bash and gawk")

BANK = ["consensus_4seqs", "consensus_7seqs", "consensus_17seqs", "consensus_25seqs", "consensus_34seqs", "consensus_354seqs",
        "consensus_3seqs", "SINE-3_Mad", "SINE-1_Mad"]


def _w(p, text):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8", newline="\n")


def _seq(rng, n=150):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mod(name):
    if str(SRC) not in sys.path:
        sys.path.insert(0, str(SRC))
    return importlib.reload(importlib.import_module(name))


# ---------------------------------------------------------------------------------------------------------------------
# 1. labels
# ---------------------------------------------------------------------------------------------------------------------
def _run(tmp_path, names):
    rng = random.Random(7)
    run = tmp_path / "run"
    _w(run / "consensuses.clean.fa", "".join(">%s\n%s\n" % (n, _seq(rng)) for n in names))
    _w(run / "results" / "array_flag.tsv",
       "family\tcopies\tcopies_in_arrays\tpct_in_arrays\tnull_pct\texcess_pct\tarrays\tmedian_spacing_bp\tflag\n"
       + "".join("%s\t10\t0\t0.0\t0.0\t0.0\t0\t0\t-\n" % n for n in names))
    return run


def test_blocks_labels_unique_in_bank(tmp_path):
    rb = _mod("report_blocks")
    run = _run(tmp_path, ["consensus_354seqs", "consensus_25seqs", "r6_210seqs", "SINE-3_Mad"])
    h = rb.section(run)
    # matrix headers, pair-view names and the tandem-array table: the two consensus_* rows are told apart
    assert '<th class="bh">consensus_354seqs</th>' in h and '<th class="bh">consensus_25seqs</th>' in h, h[:3000]
    assert '<th class="bh">consensus</th>' not in h
    assert '"consensus_354seqs", "consensus_25seqs"' in h            # D.names of the pair view
    assert "<tr><td>consensus_354seqs</td>" in h and "<tr><td>consensus_25seqs</td>" in h
    assert "<tr><td>consensus</td>" not in h
    # the legitimate short label stays
    assert '<th class="bh">r6</th>' in h and "<tr><td>r6</td>" in h
    assert "<tr><td>SINE-3_Mad</td>" in h


def test_blocks_labels_rule():
    rb = _mod("report_blocks")
    lab = rb.labels(["consensus_354seqs", "consensus_25seqs", "r6_210seqs", "r7", "r7_12seqs", "SINE-3_Mad"])
    assert lab == {"consensus_354seqs": "consensus_354seqs", "consensus_25seqs": "consensus_25seqs", "r6_210seqs": "r6",
                   "r7": "r7", "r7_12seqs": "r7_12seqs", "SINE-3_Mad": "SINE-3_Mad"}


def test_hierarchy_short_and_full_unique():
    rh = _mod("report_hierarchy")
    names = {"consensus_354seqs", "consensus_25seqs", "r6_210seqs", "r7_3seqs"}
    assert rh._short("consensus_354seqs", names) == "consensus_354seqs"
    assert rh._short("consensus_25seqs", names) == "consensus_25seqs"
    assert rh._short("r6_210seqs", names) == "r6"
    assert rh._full("r6", names) == "r6_210seqs"
    assert rh._full("consensus", names) is None                       # fits two: neither, not the first one
    assert rh._full("consensus_25seqs", names) == "consensus_25seqs"
    url = rh._plate("consensus", "mdo", "https://x/aln/", ["consensus_354seqs", "consensus_25seqs"])
    assert url is None
    assert "mdo_r6_210seqs_top100" in rh._plate("r6", "mdo", "https://x/aln/", ["r6_210seqs", "r7_3seqs"])


def test_hierarchy_bank_composite_labels():
    rh = _mod("report_hierarchy")
    comp = {"X_comp": ["composite", [(1, 120, "consensus_354seqs", "95.0", "-"), (121, 240, "consensus_25seqs", "93.0", "-")]],
            "consensus_354seqs": ["single", [(1, 120, "consensus_354seqs", "100.0", "-")]],
            "consensus_25seqs": ["single", [(1, 120, "consensus_25seqs", "100.0", "-")]]}
    ucol = {"consensus_25seqs": "#111", "consensus_354seqs": "#222"}
    svg = rh._bank_composites(comp, ucol)
    assert ">consensus_354seqs 95.0</text>" in svg and ">consensus_25seqs 93.0</text>" in svg, svg
    assert ">consensus 95.0</text>" not in svg


def _fs7_toy(d):
    _w(d / "cons.tsv", "consensus_354seqs\t200\t201\t0\nconsensus_25seqs\t150\t151\t0\nr6_210seqs\t180\t181\t0\nr7_3seqs\t170\t171\t0\n")
    hdr = "\t".join("c%d" % i for i in range(1, 16)) + "\n"
    _w(d / "peaks.tsv", hdr + "consensus_354seqs\t3\tconsensus_25seqs\tr\tt\t9\t9\t1\t200\t1\t0\t.\t.\t.\tP1\n"
                              "r6_210seqs\t3\tr7_3seqs\tr\tt\t9\t9\t1\t180\t1\t0\t.\t.\t.\tP2\n")
    _w(d / "reassign.tsv", "cand\ta\tb\tc\tpct\tverdict\nA\t0\t0\t0\t90.0\taccept\nB\t0\t0\t0\t90.0\taccept\n")
    _w(d / "ends.tsv", "cand\ta\tu5\tb\tc\tu3\td\topen\nA\t.\t-\t.\t.\t-\t.\t-\nB\t.\t-\t.\t.\t-\t.\t-\n")
    ch = "\t".join("c%d" % i for i in range(1, 15)) + "\n"
    _w(d / "candidates.tsv", ch + "A\tP1\tpair\tx\tx\t9\tx\tx\tx\tx\tx\tx\t350\tkept\n"
                                  "B\tP2\tpair\tx\tx\t9\tx\tx\tx\tx\tx\tx\t350\tkept\n")
    _w(d / "candidates.fa", ">A\nACGT\n>B\nACGT\n")


@need_shell
def test_fs7_short_names_unique_in_bank(tmp_path):
    d = tmp_path / "fs"
    _fs7_toy(d)
    r = subprocess.run([BASH, (SRC / "flankscan" / "fs7_hierarchy.sh").as_posix(), d.as_posix()], capture_output=True,
                       text=True, encoding="utf-8", errors="replace")
    assert r.returncode == 0, r.stderr
    rows = {l.split("\t")[0]: l.split("\t") for l in (d / "hierarchy.tsv").read_text(encoding="utf-8").splitlines()[1:]}
    assert rows["A"][7] == "consensus_354seqs" and rows["A"][11] == "consensus_25seqs", rows["A"]
    assert rows["A"][1] == "consensus_354seqs_consensus_25seqs_P1"
    assert rows["A"][17] == "consensus_354seqs:1-200;consensus_25seqs:1-150"
    assert rows["B"][7] == "r6" and rows["B"][11] == "r7" and rows["B"][1] == "r6_r7_P2", rows["B"]


# ---------------------------------------------------------------------------------------------------------------------
# 2 and 3. --exclude / --add (shell functions of the orchestrator)
# ---------------------------------------------------------------------------------------------------------------------
def _stub(bindir, name, body):
    bindir.mkdir(parents=True, exist_ok=True)
    p = bindir / name
    p.write_text("#!/usr/bin/env bash\n" + body + "\n", encoding="utf-8", newline="\n")
    p.chmod(p.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


def _bash_function(script, name):
    lines = Path(script).read_text(encoding="utf-8").splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith(name + "(){") or l.startswith(name + "() {"))
    end = next(i for i in range(start + 1, len(lines)) if lines[i] == "}")
    return "\n".join(lines[start:end + 1]) + "\n"


def _drive(tmp_path, body, bindir=None):
    drv = tmp_path / "drv.sh"
    _w(drv, "set -euo pipefail\nlog(){ echo \"$*\" >&2; }\ndie(){ echo \"ERROR: $*\" >&2; exit 1; }\n" + body)
    env = dict(os.environ)
    if bindir is not None:
        env["PATH"] = str(bindir) + os.pathsep + env["PATH"]
    return subprocess.run([BASH, drv.as_posix()], cwd=str(tmp_path), env=env, capture_output=True, text=True,
                          encoding="utf-8", errors="replace")


@need_shell
def test_exclude_drops_consensus_from_satellite_tables(tmp_path):
    sat = tmp_path / "run" / "genome.clean_step1" / "satellites"
    ind_hdr = ("consensus\tcons_len\tfull_hits\tkindA_loci\tkindA_monomers\tkindA_largest\tkindB_runs\tkindB_pct\tkindB_excess_pct"
               "\tkindB_long_min\tkindB_long_runs\tkindB_verified_arrays\tflag_A\tflag_B\n")
    row = "%s\t300\t100\t1\t4\t4\t11\t0.0\t0.7\t5\t2\t2\tSAT_A\t-\n"
    _w(sat / "indication.tsv", ind_hdr + "".join(row % n for n in ("SINE-1_Mad", "SINE-3_Mad", "consensus_354seqs")))
    _w(sat / "loci.bed", "#kind\tconsensus\tcontig\tstart\tend\tperiod_or_unit\tmonomers_or_hits\tsine_part\talso\n"
                         "A\tSINE-1_Mad\tc1\t10\t500\t120\t4\t1-120\t-\nA\tSINE-3_Mad\tc2\t10\t500\t120\t4\t1-120\t-\n")
    _w(sat / "excluded_hits.bed", "SINE-1_Mad\tc1\t10\t130\nSINE-3_Mad\tc2\t10\t130\n")
    _w(sat / "exclude_regions.bed", "c1\t10\t500\tA\tSINE-1_Mad\nc2\t10\t500\tA\tSINE-3_Mad\n")
    _w(sat / "SINE-1_Mad.kindA.units.fa", ">c1:10-500|p120|c4.0\nAAAA\n")
    _w(sat / "SINE-1_Mad.kindA.loci.tsv", "contig\n")
    _w(sat / "SINE-1_Mad.cons.fa", ">SINE-1_Mad\nACGT\n")
    _w(sat / "SINE-3_Mad.kindA.units.fa", ">c2:10-500|p120|c4.0\nCCCC\n")
    _w(sat / "units.fa", ">c1:10-500|p120|c4.0\nAAAA\n>c2:10-500|p120|c4.0\nCCCC\n")
    r = _drive(tmp_path, _bash_function(SRC / "SINEderella", "satellite_drop_consensus")
               + "satellite_drop_consensus %s SINE-1_Mad\n" % sat.as_posix())
    assert r.returncode == 0, r.stderr
    for f in ("indication.tsv", "loci.bed", "excluded_hits.bed", "exclude_regions.bed"):
        t = (sat / f).read_text(encoding="utf-8")
        assert "SINE-1_Mad" not in t and "SINE-3_Mad" in t, (f, t)
    assert (sat / "indication.tsv").read_text(encoding="utf-8").startswith("consensus\tcons_len")
    assert (sat / "loci.bed").read_text(encoding="utf-8").startswith("#kind")
    assert not list(sat.glob("SINE-1_Mad.*"))
    assert (sat / "units.fa").read_text(encoding="utf-8") == ">c2:10-500|p120|c4.0\nCCCC\n"
    # and the report's satellite table no longer lists it (the run's results/satellites is a link to this directory)
    rb = _mod("report_blocks")
    shutil.copytree(sat, tmp_path / "run" / "results" / "satellites")
    t = rb._satellites(tmp_path / "run")
    assert "SINE-1_Mad" not in t and "SINE-3_Mad" in t


@need_shell
def test_rebuild_allseeds_uses_this_runs_bank(tmp_path):
    rng = random.Random(3)
    run = tmp_path / "run"
    s1 = run / "genome.clean_step1"
    kept = [n for n in BANK if n != "SINE-1_Mad"]
    _w(run / "consensuses.clean.fa", "".join(">%s\n%s\n" % (n, _seq(rng)) for n in kept))
    _w(s1 / "consensus.path", (run / "consensuses.clean.fa").as_posix() + "\n")
    _w(s1 / "extracted.fasta", "".join(">c1:%d-%d(+)\n%s\n" % (i, i + 150, _seq(rng)) for i in range(0, 3000, 300)))
    # the source run's alignment, copied with genome.clean_step1: the excluded seed is a row of it
    _w(s1 / "subfam_input" / "input.clw.al", ">input_001\nACGT\n" + "".join(">%s\nACGT\n" % n for n in BANK))
    b = tmp_path / "bin"
    _stub(b, "SubFam", 'printf ">input_001\\nACGTACGT\\n>input_002\\nACGTTCGT\\n" > input.clw')
    _stub(b, "mafft", 'cat "${@: -1}"')
    r = _drive(tmp_path, "CHUNK_BP=30000\nFLANK=50\nTHREADS=1\n" + _bash_function(SRC / "SINEderella", "rebuild_allseeds")
               + "cd %s\nrebuild_allseeds \"$(cat consensus.path)\"\n" % s1.as_posix(), b)
    assert r.returncode == 0, r.stderr
    heads = [l[1:] for l in (s1 / "subfam_input" / "input.clw.al").read_text(encoding="utf-8").splitlines() if l.startswith(">")]
    assert heads == ["input_001", "input_002"] + kept, heads
    assert (s1 / "sampled_30000.fasta").read_text(encoding="utf-8") == (s1 / "extracted.fasta").read_text(encoding="utf-8")


def test_add_exclude_point_step1_at_their_own_bank():
    """the order inside the modes: consensus.path (the seeds of the rebuilt alignment) and, for --exclude, the satellite tables are
    set before rebuild_from_searches; rebuild_from_searches rebuilds the all-seed alignment from consensus.path"""
    sd = SRC / "SINEderella"
    rb = _bash_function(sd, "rebuild_from_searches")
    assert 'rebuild_allseeds "$(cat consensus.path)"' in rb
    for mode in ("mode_add", "mode_exclude"):
        body = _bash_function(sd, mode)
        i = body.index('rebuild_from_searches "$(pwd -P)/genome.clean_step1"')
        j = body.index('echo "$RUN/consensuses.clean.fa" > "$RUN/genome.clean_step1/consensus.path"')
        assert j < i, mode
    body = _bash_function(sd, "mode_exclude")
    assert body.index('satellite_drop_consensus "$RUN/genome.clean_step1/satellites" "$EXCLUDE_NAME"') < \
        body.index('rebuild_from_searches "$(pwd -P)/genome.clean_step1"')
