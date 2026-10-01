#!/usr/bin/env python3
"""Divergence difference next to the LEAK ratio - the concrete example decision D3 asks for.

A LEAK flag (step3_postprocess.sh, section 2) means the runner-up consensus scored at
least 0.90 of the best bitscore for a copy.  A close runner-up bitscore alone does not
make a plausible second family: the runner-up alignment can be far more diverged than the
best one (RepeatMasker's isTooDiverged idea; docs/BORROWED_TRIAGE.md section 5, DECIDE 3,
and decision D3 in section 8).  This tool lists the LEAK copies where the two divergences
disagree most.

Input is the step3 output table all_sines.bedlike.ALL.tsv (built in section 1 of
step3_postprocess.sh): 12 tab-separated columns, no header -

     1 chr              2 start0           3 end              4 subfam_from_extracted
     5 assigned_subfam  6 strand           7 best_bitscore    8 leak_flag
     9 conflict_flag   10 note            11 sim_bitscore    12 sim_ratio

The runner-up lives in the col-10 note as runner=<subfam>;runner_bs=<bitscore>;
runner_ratio=<runner_bs/best_bs>; col 8 is LEAK when runner_ratio >= 0.90.  Copy ids are
built the way step3 builds SeqID for unassigned.tsv: chr:start-end(strand), 1-based.

Divergence: the table carries bitscores only - no alignment identity and no alignment
length - so the tool says so and falls back to the columns that do exist (the scores):
divergence proxy = 1 - bitscore / consensus self bitscore, with each consensus's self
bitscore recovered from its copies as sim_bitscore / sim_ratio (cols 11-12; sim_ratio is
sim_bitscore / self).  A table that does carry identity (note keys best_id= / runner_id=,
percent or fraction) gets the true divergence 1 - identity instead.

Usage:
    python tools/leak_div_example.py STEP3_TABLE [--top N] [--selftest]

--selftest builds a tiny table in the real column layout in a temp dir, runs the same
logic on it and prints SELFTEST OK.
"""

import argparse
import os
import sys
import tempfile
from statistics import median

LEAK_RATIO_MIN = 0.90   # step3_postprocess.sh: runner_ratio >= 0.90 -> col 8 = LEAK
DIFF_RULE_PTS = 5.0     # "more than 5 points" difference rule, in percentage points

# 0-based columns of all_sines.bedlike.ALL.tsv
(COL_CHR, COL_START0, COL_END, COL_HDR_SUBFAM, COL_SUBFAM, COL_STRAND,
 COL_BEST_BS, COL_LEAK, COL_CONFLICT, COL_NOTE, COL_SIM_BS, COL_SIM_RATIO) = range(12)
NCOLS = 12


def parse_float(text):
    """Numeric field or None ('', '.', 'NA' and friends mean missing)."""
    if text is None:
        return None
    text = text.strip()
    if text in ("", ".", "NA", "None", "nan"):
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_note(note):
    """The col-10 note: ';'-separated key=value tags -> dict (first occurrence wins)."""
    tags = {}
    for part in (note or "").split(";"):
        if "=" in part:
            key, _, value = part.partition("=")
            tags.setdefault(key.strip(), value.strip())
    return tags


def identity_fraction(text):
    """Identity as a fraction: 97.5 (percent) or 0.975 (fraction) -> 0.975."""
    value = parse_float(text)
    if value is None:
        return None
    if value > 1.0:
        value /= 100.0
    return value


def copy_id(chr_, start0, end, strand):
    """Copy id the way step3_postprocess.sh builds SeqID: chr:start-end(strand)."""
    try:
        start = int(start0) + 1
    except (TypeError, ValueError):
        start = start0
    if strand not in ("+", "-"):
        strand = "+"
    return "%s:%s-%s(%s)" % (chr_, start, end, strand)


def read_table(path):
    """Read the step3 table -> list of row dicts; exit with a message if unusable."""
    if not os.path.isfile(path):
        sys.exit("ERROR: no such table: %s" % path)
    rows = []
    maxfields = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            maxfields = max(maxfields, len(fields))
            if len(fields) < NCOLS:
                fields += [""] * (NCOLS - len(fields))
            tags = parse_note(fields[COL_NOTE])
            rows.append({
                "line": lineno,
                "chr": fields[COL_CHR],
                "start0": fields[COL_START0],
                "end": fields[COL_END],
                "hdr_subfam": fields[COL_HDR_SUBFAM],
                "subfam": fields[COL_SUBFAM],
                "strand": fields[COL_STRAND],
                "best_bs": parse_float(fields[COL_BEST_BS]),
                "leak_flag": fields[COL_LEAK],
                "sim_bs": parse_float(fields[COL_SIM_BS]),
                "sim_ratio": parse_float(fields[COL_SIM_RATIO]),
                "tags": tags,
                "runner": tags.get("runner"),
                "runner_bs": parse_float(tags.get("runner_bs")),
                "runner_ratio": parse_float(tags.get("runner_ratio")),
            })
    if not rows:
        sys.exit("ERROR: no data rows in %s" % path)
    if maxfields < NCOLS:
        sys.exit("ERROR: %s: expected %d tab-separated columns "
                 "(all_sines.bedlike.ALL.tsv layout), found at most %d"
                 % (path, NCOLS, maxfields))
    return rows


def leak_ratio(row):
    """Runner-up / best bitscore ratio: the note value, or recomputed from the scores."""
    if row["runner_ratio"] is not None:
        return row["runner_ratio"]
    if row["runner_bs"] is not None and row["best_bs"]:
        return row["runner_bs"] / row["best_bs"]
    return None


def is_leak(row):
    """Col-8 flag, or the ratio definition (>= 0.90) when the flag is absent."""
    if row["leak_flag"] == "LEAK":
        return True
    ratio = leak_ratio(row)
    return ratio is not None and ratio >= LEAK_RATIO_MIN


def self_bits_by_subfam(rows):
    """Consensus self bitscore per subfam, recovered as sim_bitscore / sim_ratio."""
    observed = {}
    for row in rows:
        if row["sim_bs"] is None or not row["sim_ratio"]:
            continue
        observed.setdefault(row["subfam"], []).append(row["sim_bs"] / row["sim_ratio"])
    return {subfam: median(values) for subfam, values in observed.items()}


def clamp01(value):
    return max(0.0, min(1.0, value))


def divergences(row, self_bits):
    """(div_best, div_runner, source) for one row; None where not computable.

    Identity first (note keys best_id= / runner_id=); otherwise the score fallback:
    div = 1 - bitscore / consensus self bitscore.
    """
    best_id = identity_fraction(row["tags"].get("best_id"))
    runner_id = identity_fraction(row["tags"].get("runner_id"))
    if best_id is not None and runner_id is not None:
        return clamp01(1.0 - best_id), clamp01(1.0 - runner_id), "identity"

    div_best = None
    if row["sim_ratio"] is not None:
        div_best = 1.0 - row["sim_ratio"]
    elif row["best_bs"] is not None and row["subfam"] in self_bits:
        div_best = 1.0 - row["best_bs"] / self_bits[row["subfam"]]

    div_runner = None
    if row["runner_bs"] is not None and row["runner"] in self_bits:
        div_runner = 1.0 - row["runner_bs"] / self_bits[row["runner"]]

    return (clamp01(div_best) if div_best is not None else None,
            clamp01(div_runner) if div_runner is not None else None,
            "score")


def analyze(rows):
    """LEAK copies, their divergences, the ranking and the difference-rule counts."""
    self_bits = self_bits_by_subfam(rows)
    leaks = []
    for row in rows:
        if not is_leak(row):
            continue
        div_best, div_runner, source = divergences(row, self_bits)
        diff = None
        if div_best is not None and div_runner is not None:
            diff = div_runner - div_best
        leaks.append({
            "copy_id": copy_id(row["chr"], row["start0"], row["end"], row["strand"]),
            "best_fam": row["subfam"] or row["hdr_subfam"] or ".",
            "runner_fam": row["runner"] or ".",
            "best_bs": row["best_bs"],
            "runner_bs": row["runner_bs"],
            "ratio": leak_ratio(row),
            "div_best": div_best,
            "div_runner": div_runner,
            "diff": diff,
            "source": source,
        })
    ranked = [rec for rec in leaks if rec["diff"] is not None]
    ranked.sort(key=lambda rec: (-rec["diff"], rec["copy_id"]))
    rule = DIFF_RULE_PTS / 100.0
    removed = sum(1 for rec in ranked if rec["diff"] > rule)
    total = len(leaks)
    return {
        "self_bits": self_bits,
        "leaks": leaks,
        "ranked": ranked,
        "removed": removed,
        "total": total,
        "no_data": total - len(ranked),
        "share": (removed / total) if total else 0.0,
        "identity_rows": sum(1 for rec in leaks if rec["source"] == "identity"),
    }


def _fmt(value, spec):
    return spec % value if value is not None else "."


def format_table(header, rows):
    """Plain fixed-width table, columns separated by two spaces."""
    widths = [len(cell) for cell in header]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))
    lines = ["  ".join(cell.ljust(widths[i]) for i, cell in enumerate(header))]
    for row in rows:
        lines.append("  ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)))
    return "\n".join(lines)


def report(path, result, top):
    total = result["total"]
    print("LEAK copies (runner-up/best bitscore ratio >= %.2f) in %s: %d"
          % (LEAK_RATIO_MIN, path, total))
    if not total:
        print("No LEAK copies - nothing to compare.")
        return
    fallback_rows = total - result["identity_rows"]
    if fallback_rows:
        print()
        if result["identity_rows"] == 0:
            print("NOTE: %s provides no alignment identity (and no alignment" % path)
            print("length) for the best and runner-up hits - saying so and falling")
        else:
            print("NOTE: %d of the LEAK rows carry no alignment identity - for those,"
                  % fallback_rows)
        print("back to the columns that do exist (the scores): divergence proxy =")
        print("1 - bitscore / consensus self bitscore, with each consensus's self")
        print("bitscore recovered from its copies as sim_bitscore / sim_ratio (cols")
        print("11-12; %d subfamilies here). True divergence (1 - identity) would need"
              % len(result["self_bits"]))
        print("the alignments themselves.")
    shown = result["ranked"][:max(0, top)]
    print()
    if shown:
        print("Top %d LEAK copies by divergence difference (runner-up more diverged than the best):"
              % len(shown))
        header = ["copy_id", "best_fam", "runner_fam", "best_bs", "runner_bs",
                  "ratio", "div_best%", "div_runner%", "diff_pts"]
        table = []
        for rec in shown:
            table.append([
                rec["copy_id"],
                rec["best_fam"],
                rec["runner_fam"],
                _fmt(rec["best_bs"], "%.1f"),
                _fmt(rec["runner_bs"], "%.1f"),
                _fmt(rec["ratio"], "%.4f"),
                "%.2f" % (100.0 * rec["div_best"]),
                "%.2f" % (100.0 * rec["div_runner"]),
                "%.2f" % (100.0 * rec["diff"]),
            ])
        print(format_table(header, table))
    elif result["ranked"]:
        print("Nothing to show (--top %d)." % top)
    else:
        print("No LEAK copy has both divergences - cannot rank any.")
    print()
    print("Summary: a divergence-difference rule of more than %.0f points (a runner-up"
          % DIFF_RULE_PTS)
    print("alignment that much more diverged than the best one is not a real leak)")
    print("would remove the LEAK flag from %d of %d copies, i.e. %.1f%% of the LEAK copies."
          % (result["removed"], total, 100.0 * result["share"]))
    real = sum(1 for rec in result["ranked"] if rec["diff"] < 0.0)
    if real:
        print("In %d LEAK copies the runner-up is not more diverged than the best -"
              % real)
        print("those leaks look real.")
    if result["no_data"]:
        print("%d LEAK copy%s no data for both divergences (runner-up family without a"
              % (result["no_data"], " has" if result["no_data"] == 1 else "s have"))
        print("recoverable self bitscore, or missing bitscores) and keeps the flag.")


def run(path, top):
    rows = read_table(path)
    result = analyze(rows)
    report(path, result, top)
    return result


# Tiny table in the real 12-column layout of all_sines.bedlike.ALL.tsv.
# Self bitscores recoverable from cols 11-12: S1=200, S2=200, S3=400; S9 never assigned.
SELFTEST_ROWS = [
    # chr    start0  end    hdr  assigned strand best_bs leak conflict note sim_bs sim_ratio
    ("chr1", "100", "300", "S1", "S1", "+", "180.0", "LEAK", ".",
     "status=assigned;votes=10/10;thr=81.0;runner=S3;runner_bs=170.0;runner_ratio=0.9444",
     "180.0", "0.9000"),
    ("chr1", "500", "700", "S2", "S2", "-", "190.0", "LEAK", ".",
     "status=assigned;votes=10/10;thr=85.5;runner=S1;runner_bs=185.0;runner_ratio=0.9737",
     "190.0", "0.9500"),
    ("chr2", "900", "1100", "S1", "S1", "+", "170.0", "LEAK", ".",
     "status=assigned;votes=10/10;thr=76.5;runner=S2;runner_bs=155.0;runner_ratio=0.9118",
     "170.0", "0.8500"),
    ("chr2", "1300", "1500", "S2", "S2", "+", "200.0", ".", ".",
     "status=assigned;votes=10/10;thr=90.0;runner=S1;runner_bs=150.0;runner_ratio=0.7500",
     "200.0", "1.0000"),
    ("chr3", "10", "210", "S3", "S3", "+", "380.0", "LEAK", ".",
     "status=assigned;votes=10/10;thr=171.0;runner=S9;runner_bs=370.0;runner_ratio=0.9737",
     "380.0", "0.9500"),
]


def selftest():
    """Build a tiny table in the real column layout, run the same logic, check it."""
    with tempfile.TemporaryDirectory(prefix="leak_div_example_") as tmp:
        path = os.path.join(tmp, "all_sines.bedlike.ALL.tsv")
        with open(path, "w") as fh:
            for row in SELFTEST_ROWS:
                fh.write("\t".join(row) + "\n")
        print("SELFTEST: tiny step3 table at %s" % path)
        result = run(path, top=10)

    # 4 LEAK copies (rows 1, 2, 3, 5); row 4 is below the 0.90 ratio.
    assert result["total"] == 4, result["total"]
    # self bitscores recovered from sim_bitscore / sim_ratio
    for subfam, expect in (("S1", 200.0), ("S2", 200.0), ("S3", 400.0)):
        got = result["self_bits"][subfam]
        assert abs(got - expect) < 1e-6, (subfam, got)
    assert "S9" not in result["self_bits"]          # no copy assigned to S9
    # ranking by divergence difference (runner-up more diverged than best)
    ranked = result["ranked"]
    assert len(ranked) == 3, len(ranked)             # the S9 runner-up has no self bitscore
    assert ranked[0]["copy_id"] == "chr1:101-300(+)"
    assert abs(ranked[0]["diff"] - 0.475) < 1e-6, ranked[0]["diff"]
    assert ranked[1]["copy_id"] == "chr2:901-1100(+)"
    assert abs(ranked[1]["diff"] - 0.075) < 1e-6, ranked[1]["diff"]
    assert ranked[2]["copy_id"] == "chr1:501-700(-)"
    assert abs(ranked[2]["diff"] - 0.025) < 1e-6, ranked[2]["diff"]
    # a >5 points difference rule would unflag 2 of the 4 LEAK copies
    assert result["removed"] == 2, result["removed"]
    assert result["no_data"] == 1, result["no_data"]
    assert abs(result["share"] - 0.5) < 1e-9, result["share"]

    # identity path (a table that does provide identity in the note)
    row = {"subfam": "S1", "runner": "S2", "best_bs": 180.0, "runner_bs": 170.0,
           "runner_ratio": 0.9444, "sim_bs": 180.0, "sim_ratio": 0.9,
           "tags": {"best_id": "97.5", "runner_id": "0.80"}}
    div_best, div_runner, source = divergences(row, {})
    assert source == "identity", source
    assert abs(div_best - 0.025) < 1e-9, div_best     # 97.5 percent identity
    assert abs(div_runner - 0.20) < 1e-9, div_runner  # 0.80 fraction identity

    print("SELFTEST OK")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="LEAK copies where the runner-up alignment is much more diverged "
                    "than the best one (divergence difference next to the LEAK ratio).")
    parser.add_argument("table", nargs="?",
                        help="step3 table (all_sines.bedlike.ALL.tsv)")
    parser.add_argument("--top", type=int, default=10, metavar="N",
                        help="how many disagreeing copies to show (default: %(default)s)")
    parser.add_argument("--selftest", action="store_true",
                        help="build a tiny table in the real column layout, run the "
                             "same logic on it and print SELFTEST OK")
    args = parser.parse_args(argv)
    if args.selftest:
        selftest()
        return 0
    if not args.table:
        parser.error("STEP3_TABLE is required unless --selftest is given")
    run(args.table, args.top)
    return 0


if __name__ == "__main__":
    sys.exit(main())
