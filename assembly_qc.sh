#!/usr/bin/env bash
# assembly_qc.sh — Genome assembly quality diagnostics for SINE searches
#
# Detects assemblies likely to produce false SINE hits due to:
#   - N-rich scaffold edges (padding artifacts)
#   - Low contiguity (many short scaffolds)
#   - High gap fraction (excessive unknown bases)
#
# Usage: assembly_qc.sh <genome.fa[.gz]> [output_prefix]
#
# Output: TSV with assembly stats + quality verdict
#   Columns: metric, value, flag
#   Verdict: SOLID / CAUTION / WARN
#
# One pass over the sequence. Scaffold lengths are sorted with the external
# sort. An earlier version appended every line onto one string and then
# bubble-sorted every scaffold length inside awk; both are quadratic, and on
# a fragmented assembly the sort never finishes.

set -euo pipefail

usage() {
    echo "Usage: $0 <genome.fa[.gz]> [output_prefix]" >&2
    echo "" >&2
    echo "Analyzes genome assembly quality for SINE search reliability." >&2
    echo "Flags assemblies with N-rich edges, low N50, or high gap fraction." >&2
    exit 1
}

[[ $# -lt 1 ]] && usage
GENOME="$1"
PREFIX="${2:-${GENOME%.gz}}"
PREFIX="${PREFIX%.fa}"
PREFIX="${PREFIX%.fasta}"
OUT="${PREFIX}_assembly_qc.tsv"

[[ ! -f "$GENOME" ]] && { echo "ERROR: File not found: $GENOME" >&2; exit 1; }

if [[ "$GENOME" == *.gz ]]; then
    CAT="zcat"
else
    CAT="cat"
fi

for cmd in awk sort; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not found" >&2; exit 1; }
done

echo "Analyzing assembly: $GENOME" >&2

WORK="${TMPDIR:-$HOME/tmp}"
mkdir -p "$WORK"
LENS="$WORK/assembly_qc_lengths.$$"
STATS="$WORK/assembly_qc_stats.$$"
trap 'rm -f "$LENS" "$STATS"' EXIT

$CAT "$GENOME" | awk -v stats="$STATS" '
BEGIN {
    EDGE = 500; seen = 0; largest = 0
    # Every counter printed to STATS must start at 0: an unset awk variable prints as an
    # empty field and shifts the read below (no short scaffolds -> largest read as 0).
    n_scaffolds = 0; total_bases = 0; total_N = 0; edge_N = 0; edge_total = 0; short_scaffolds = 0
}

function finish(   left_len, right_len, left, right, n_in_left, n_in_right) {
    n_scaffolds++
    total_bases += len
    total_N += ncount
    if (len < 1000) short_scaffolds++
    if (len > largest) largest = len
    print len

    left_len = (len < EDGE) ? len : EDGE
    right_len = (len < EDGE) ? 0 : ((len < 2 * EDGE) ? len - EDGE : EDGE)

    if (left_len > 0) {
        left = substr(head, 1, left_len)
        edge_total += left_len
        n_in_left = gsub(/[nN]/, "", left)
        edge_N += n_in_left
    }
    if (right_len > 0) {
        right = substr(tail, length(tail) - right_len + 1, right_len)
        edge_total += right_len
        n_in_right = gsub(/[nN]/, "", right)
        edge_N += n_in_right
    }
}

function reset() {
    len = 0
    ncount = 0
    head = ""
    tail = ""
}

/^>/ {
    if (seen) finish()
    seen = 1
    reset()
    next
}

{
    if (!seen) next
    line = $0
    L = length(line)
    if (L == 0) next
    tmp = line
    ncount += gsub(/[nN]/, "", tmp)
    len += L
    if (length(head) < EDGE) {
        head = head line
        if (length(head) > EDGE) head = substr(head, 1, EDGE)
    }
    tail = tail line
    if (length(tail) > EDGE) tail = substr(tail, length(tail) - EDGE + 1)
}

END {
    if (seen) finish()
    print n_scaffolds, total_bases, total_N, edge_N, edge_total, short_scaffolds, largest > stats
}
' > "$LENS"

read -r n_scaffolds total_bases total_N edge_N edge_total short_scaffolds largest < "$STATS"

n50=0
l50=0
if [[ "$n_scaffolds" -gt 0 && "$total_bases" -gt 0 ]]; then
  read -r l50 n50 < <(sort -nr "$LENS" | awk -v half="$total_bases" '
    BEGIN { half = half / 2 }
    { cumul += $1; if (cumul >= half) { print NR, $1; exit } }
  ')
fi

total_ACGT=$((total_bases - total_N))

awk -v n_scaffolds="$n_scaffolds" \
    -v total_bases="$total_bases" \
    -v total_ACGT="$total_ACGT" \
    -v total_N="$total_N" \
    -v edge_N="$edge_N" \
    -v edge_total="$edge_total" \
    -v short_scaffolds="$short_scaffolds" \
    -v largest="$largest" \
    -v n50="$n50" \
    -v l50="$l50" '
BEGIN {
    gap_pct = (total_bases > 0) ? 100.0 * total_N / total_bases : 0
    edge_n_pct = (edge_total > 0) ? 100.0 * edge_N / edge_total : 0
    short_pct = (n_scaffolds > 0) ? 100.0 * short_scaffolds / n_scaffolds : 0

    if (n50 >= 10000000) n50_flag = "ok"
    else if (n50 >= 1000000) n50_flag = "caution"
    else n50_flag = "warn"

    if (gap_pct < 5) gap_flag = "ok"
    else if (gap_pct < 15) gap_flag = "caution"
    else gap_flag = "warn"

    if (edge_n_pct < 5) edge_flag = "ok"
    else if (edge_n_pct < 20) edge_flag = "caution"
    else edge_flag = "warn"

    if (short_pct < 5) short_flag = "ok"
    else if (short_pct < 20) short_flag = "caution"
    else short_flag = "warn"

    warns = 0; cautions = 0
    if (n50_flag == "warn") warns++
    if (gap_flag == "warn") warns++
    if (edge_flag == "warn") warns++
    if (short_flag == "warn") warns++
    if (n50_flag == "caution") cautions++
    if (gap_flag == "caution") cautions++
    if (edge_flag == "caution") cautions++
    if (short_flag == "caution") cautions++

    if (warns > 0) verdict = "WARN"
    else if (cautions >= 2) verdict = "CAUTION"
    else verdict = "SOLID"

    print "metric\tvalue\tflag"
    print "scaffolds\t" n_scaffolds "\t-"
    printf "total_bases\t%d\t-\n", total_bases
    printf "total_ACGT\t%d\t-\n", total_ACGT
    printf "total_N\t%d\t-\n", total_N
    printf "gap_pct\t%.2f\t%s\n", gap_pct, gap_flag
    printf "N50\t%d\t%s\n", n50, n50_flag
    print "L50\t" l50 "\t-"
    printf "largest_scaffold\t%d\t-\n", largest
    printf "short_scaffolds_pct\t%.2f\t%s\n", short_pct, short_flag
    printf "edge_N_pct\t%.2f\t%s\n", edge_n_pct, edge_flag
    print "verdict\t" verdict "\t" verdict
}
' > "$OUT"

echo "Results written to: $OUT" >&2
echo "" >&2
cat "$OUT" >&2
