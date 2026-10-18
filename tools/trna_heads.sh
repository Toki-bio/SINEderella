#!/usr/bin/env bash
# trna_heads.sh RUN_ROOT [THREADS=8]
#
# Where the heads come from: tRNA, 5S rRNA or 7SL RNA, for every consensus of the searched bank and every stretch of its
# composites (step2_output/composite_units.tsv, tools/composite_units.py), so each unit's head is tested on its own.
# Method (knowlege methods/trna-similarity-scans-for-sines.md; skill trna-similarity-scan; tested on rsi 2026-10-18):
#   1. ssearch36 against a set of tRNA genes (TRNA_REF, e.g. GtRNAdb hg19, 419 genes): best tRNA, identity, length,
#      bits; SINEBase's rule (>= 60 % over >= 60 nt); the best hit of a different amino acid (isotype margin)
#   2. null: SHUF_N (20) dinucleotide-preserving shuffles of every query against the same set (esl-shuffle -d, else a
#      plain shuffle); a head counts when its bits exceed the best shuffle
#   3. Infernal cmscan against Rfam (RFAM_CM; tRNA, 5S, 7SL, ...), E <= 1 - optional
#   4. tRNAscan-SE -E -H --detail -X 10 (isotype-model scores, pseudo / IPD flags) - optional
# Rfam and tRNAscan-SE miss degraded heads (rsi: r2-r4, r6); the pairwise search with its null does not.
#
# Env: TRNA_REF (required; without it the step is skipped), RFAM_CM, RNATOOLS_BIN (a bin dir with cmscan / tRNAscan-SE,
#      e.g. a separate conda env), SHUF_N=20
# Out: <step2_output>/trna/ (queries.fa, ss.m8, shuf.m8, cm.tbl, trnascan.out) and <step2_output>/trna/trna_heads.tsv
set -uo pipefail
RUN=$(readlink -f "${1:?RUN_ROOT}"); T=${2:-8}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -n "${TRNA_REF:-}" && -s "${TRNA_REF:-}" ]] || { echo "trna_heads: TRNA_REF not set or missing - skipped"; exit 0; }
S2=$(ls -dt "$RUN"/step2/step2_output* 2>/dev/null | head -n1)
[[ -n "$S2" && -s "$RUN/consensuses.clean.fa" ]] || { echo "trna_heads: no step2 output or bank in $RUN" >&2; exit 1; }
OUT="$S2/trna"; mkdir -p "$OUT"; cd "$OUT" || exit 1
python3 - "$RUN" "$S2" <<'PY'
import csv, os, sys
run, s2 = sys.argv[1], sys.argv[2]
seq, n = {}, None
for l in open(os.path.join(run, "consensuses.clean.fa")):
    l = l.strip()
    if l.startswith(">"): n = l[1:].split()[0]; seq[n] = ""
    elif n: seq[n] += l.upper()
with open("queries.fa", "w") as fo:
    for k, s in seq.items():
        fo.write(">%s\n%s\n" % (k, s))
    cu = os.path.join(s2, "composite_units.tsv")
    if os.path.isfile(cu):
        for r in csv.DictReader(open(cu), delimiter="\t"):
            a, b = int(r["start"]), int(r["end"])
            if r["kind"] == "composite" and b - a + 1 >= 30:
                fo.write(">%s:%d-%d|%s\n%s\n" % (r["consensus"], a, b, r["unit"], seq[r["consensus"]][a - 1:b]))
PY
cp -f "$TRNA_REF" ref.fa
ssearch36 -3 -m 8 -E 10 -T "$T" queries.fa ref.fa > ss.m8 2>/dev/null || { echo "trna_heads: ssearch36 failed" >&2; exit 1; }
N=${SHUF_N:-20}
if command -v esl-shuffle > /dev/null; then
  esl-shuffle -d -N "$N" --seed 42 queries.fa > shuf.fa 2>/dev/null
else
  python3 - "$N" <<'PY'
import random, sys
random.seed(42); n = int(sys.argv[1]); recs = []
for l in open("queries.fa"):
    l = l.strip()
    if l.startswith(">"): recs.append([l[1:], ""])
    else: recs[-1][1] += l
with open("shuf.fa", "w") as fo:
    for name, s in recs:
        for k in range(n):
            x = list(s); random.shuffle(x); fo.write(">%s-shuffled-%d\n%s\n" % (name, k + 1, "".join(x)))
PY
fi
ssearch36 -3 -m 8 -E 10 -T "$T" shuf.fa ref.fa > shuf.m8 2>/dev/null || { echo "trna_heads: ssearch36 (null) failed" >&2; exit 1; }
PATHX="${RNATOOLS_BIN:+$RNATOOLS_BIN:}$PATH"
: > cm.tbl; : > trnascan.out
if [[ -n "${RFAM_CM:-}" && -s "${RFAM_CM:-}" ]] && PATH="$PATHX" command -v cmscan > /dev/null; then
  PATH="$PATHX" cmscan --cpu "$T" -E 1 --tblout cm.tbl -o /dev/null "$RFAM_CM" queries.fa || echo "trna_heads: cmscan failed (Rfam column empty)" >&2
else
  echo "trna_heads: no RFAM_CM / cmscan - Rfam column empty" >&2
fi
if PATH="$PATHX" command -v tRNAscan-SE > /dev/null; then
  PATH="$PATHX" tRNAscan-SE -E -H --detail -X 10 -q -o trnascan.out -f trnascan.ss queries.fa > /dev/null 2>&1 \
    || echo "trna_heads: tRNAscan-SE failed (column empty)" >&2
else
  echo "trna_heads: no tRNAscan-SE - column empty" >&2
fi
python3 "$HERE/trna_heads_summary.py" > trna_heads.tsv || { echo "trna_heads: summary failed" >&2; exit 1; }
echo "TRNA_HEADS queries=$(grep -c '>' queries.fa) -> $OUT/trna_heads.tsv"
