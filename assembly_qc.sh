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

# Determine cat command for gzipped input
if [[ "$GENOME" == *.gz ]]; then
    CAT="zcat"
else
    CAT="cat"
fi

# Check dependencies
for cmd in awk; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd not found" >&2; exit 1; }
done

echo "Analyzing assembly: $GENOME" >&2

# Single-pass awk analysis
$CAT "$GENOME" | awk '
BEGIN {
    n_scaffolds = 0
    total_bases = 0
    total_N = 0
    total_ACGT = 0
    edge_N = 0          # N bases in first/last 500bp of each scaffold
    edge_total = 0       # total bases in edge regions
    short_scaffolds = 0  # scaffolds < 1000bp
    EDGE = 500           # edge window size
}

/^>/ {
    # Process previous scaffold
    if (n_scaffolds > 0) {
        process_scaffold()
    }
    n_scaffolds++
    seq = ""
    next
}

{
    seq = seq $0
}

function process_scaffold() {
    len = length(seq)
    total_bases += len

    if (len < 1000) short_scaffolds++

    # Store lengths for N50
    scaffold_lengths[n_scaffolds] = len

    # Count Ns and ACGT in full sequence
    n_count = gsub(/[nN]/, "&", seq)
    total_N += n_count
    total_ACGT += (len - n_count)

    # Edge analysis: first EDGE bp and last EDGE bp
    left_len = (len < EDGE) ? len : EDGE
    right_len = (len < EDGE) ? 0 : ((len < 2*EDGE) ? len - EDGE : EDGE)

    left = substr(seq, 1, left_len)
    edge_total += left_len
    n_in_left = gsub(/[nN]/, "&", left)
    edge_N += n_in_left

    if (right_len > 0) {
        right = substr(seq, len - right_len + 1, right_len)
        edge_total += right_len
        n_in_right = gsub(/[nN]/, "&", right)
        edge_N += n_in_right
    }
}

END {
    # Process last scaffold
    if (n_scaffolds > 0) process_scaffold()

    # Sort scaffold lengths descending for N50
    n = n_scaffolds
    for (i = 1; i <= n; i++) sorted[i] = scaffold_lengths[i]
    for (i = 1; i <= n; i++) {
        for (j = i+1; j <= n; j++) {
            if (sorted[j] > sorted[i]) {
                tmp = sorted[i]; sorted[i] = sorted[j]; sorted[j] = tmp
            }
        }
    }

    # N50 calculation
    half = total_bases / 2
    cumul = 0
    n50 = 0
    l50 = 0
    for (i = 1; i <= n; i++) {
        cumul += sorted[i]
        if (cumul >= half) {
            n50 = sorted[i]
            l50 = i
            break
        }
    }

    # Largest scaffold
    largest = (n > 0) ? sorted[1] : 0

    # Metrics
    gap_pct = (total_bases > 0) ? 100.0 * total_N / total_bases : 0
    edge_n_pct = (edge_total > 0) ? 100.0 * edge_N / edge_total : 0
    short_pct = (n_scaffolds > 0) ? 100.0 * short_scaffolds / n_scaffolds : 0

    # Quality flags
    # N50 flag
    if (n50 >= 10000000) n50_flag = "ok"
    else if (n50 >= 1000000) n50_flag = "caution"
    else n50_flag = "warn"

    # Gap fraction flag
    if (gap_pct < 5) gap_flag = "ok"
    else if (gap_pct < 15) gap_flag = "caution"
    else gap_flag = "warn"

    # Edge N flag (the Anilios bituberculatus signature)
    if (edge_n_pct < 5) edge_flag = "ok"
    else if (edge_n_pct < 20) edge_flag = "caution"
    else edge_flag = "warn"

    # Short scaffolds flag
    if (short_pct < 5) short_flag = "ok"
    else if (short_pct < 20) short_flag = "caution"
    else short_flag = "warn"

    # Overall verdict
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

    # Output
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
