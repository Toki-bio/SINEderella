#!/usr/bin/env bash
# composite_pieces.sh UNITS.tsv BANK.fa OUTDIR [MIN_LEN=50]
#
# The stretches of the bank's composites as search queries of their own, for step 1 (2026-10-26, docs/PARTIAL_COPIES.md).
#
# sear keeps a hit only when it covers >= 80 % of its query, so a copy of one unit of a composite (rsi: the head that
# r1_r3_P1, r1_r3_P2 and C11 share) is found only if it is also 80 % of some other bank entry. Searching every composite
# stretch of >= MIN_LEN bp that a unit explains (rows of kind "composite" in UNITS.tsv, from composite_units.py --bank)
# finds them: rsi without r1_9seqs, all 21 r1 loci, the 3 the bank missed included. Stretches no unit explains ("-") are
# not searched (2026-11-27): searched alone, the 69 bp tail of Atelerix e2-3longB (281-349, the shared 3' tail of e2
# copies) found 195,644 loci of its own, all unassigned.
#
# Writes, per stretch, OUTDIR/pieces.part_<composite>__piece_<start>-<end>.fa whose header's first word is the COMPOSITE's
# name, so every place that labels hits by the query header (step 1, the orchestrator's rebuild, --exclude) gives piece hits
# the composite's label; step 2 votes with the bank only. OUTDIR/pieces.named.fa holds the same stretches under their file
# names (<composite>__piece_<start>-<end>) for the satellite screen, which names queries by file.
# Prints the number of pieces. awk rewrite (2026-11-27) of composite_pieces.py: the same files, byte for byte.
set -euo pipefail
[[ $# -ge 3 ]] || { echo "usage: composite_pieces.sh UNITS.tsv BANK.fa OUTDIR [MIN_LEN=50]" >&2; exit 2; }
UNITS=$1; BANK=$2; OUT=$3; MINLEN=${4:-50}
[[ $MINLEN =~ ^[0-9]+$ ]] || { echo "composite_pieces.sh: MIN_LEN must be a whole number" >&2; exit 2; }
[[ -r $UNITS && -r $BANK ]] || { echo "composite_pieces.sh: cannot read $UNITS or $BANK" >&2; exit 1; }
mkdir -p "$OUT"
LC_ALL=C awk -v OUT="$OUT" -v MINLEN="$MINLEN" '
function trim(s) { sub(/^[ \t\r\n\f\v]+/, "", s); sub(/[ \t\r\n\f\v]+$/, "", s); return s }
FNR == NR {                                   # bank: name = first word of the header; a repeated name starts over
    line = trim($0)
    if (substr(line, 1, 1) == ">") { name = substr(line, 2); sub(/[ \t].*$/, "", name); seq[name] = ""; have = 1; next }
    if (have) seq[name] = seq[name] line
    next
}
FNR == 1 { named = OUT "/pieces.named.fa"; printf "" > named; next }        # header row of the units table
{
    sub(/[\r\n]+$/, "")
    nf = split($0, f, "\t")
    if (nf < 5 || f[2] != "composite" || !(f[1] in seq)) next
    s = f[3] + 0; e = f[4] + 0
    if (e - s + 1 < MINLEN) next
    if (f[5] == "-" || f[5] == "") next      # a stretch no unit explains is not a module
    piece = substr(seq[f[1]], s, e - s + 1)
    t = f[1]; gsub(/[^A-Za-z0-9_.#+-]/, "_", t)
    tag = t "__piece_" s "-" e
    file = OUT "/pieces.part_" tag ".fa"
    printf ">%s piece_%d-%d unit=%s\n%s\n", f[1], s, e, f[5], piece > file; close(file)
    printf ">%s\n%s\n", tag, piece > named
    n++
}
END { if (named == "") { named = OUT "/pieces.named.fa"; printf "" > named } close(named); print n + 0 }
' "$BANK" "$UNITS"
