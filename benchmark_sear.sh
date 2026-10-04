#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# benchmark_sear.sh — Compare single-query (sear) vs multi-query (sear_multi)
#
# Usage:
#   benchmark_sear.sh <genomes.tsv> <consensus.fa> [max_genome_mb] [n_queries]
#
# genomes.tsv: tab-separated, col1=Species_Name, col2=/path/to/genome.fa
#              First 3 species are used (or fewer if list is shorter).
#
# Defaults: 100 MB genome cap, 5 queries from the consensus.
#
# What it does:
#   1. Truncates each genome to ~max_genome_mb
#   2. Selects n_queries from the consensus
#   3. Runs OLD approach: split consensus + sear per query (with -k)
#   4. Runs NEW approach: sear_multi with full consensus
#   5. Compares: wall time, hit counts, BED identity
#
# Requirements: sear, sear_multi, ssearch36, bedtools, seqkit, samtools in PATH
###############################################################################

GENOMES_TSV="${1:?Usage: benchmark_sear.sh <genomes.tsv> <consensus.fa> [max_mb] [n_queries]}"
CONSENSUS="${2:?Usage: benchmark_sear.sh <genomes.tsv> <consensus.fa> [max_mb] [n_queries]}"
MAX_MB="${3:-100}"
N_QUERIES="${4:-5}"

MAX_BYTES=$((MAX_MB * 1000000))

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

# Find sear: prefer SCRIPT_DIR, fall back to PATH
if [[ -f "$SCRIPT_DIR/sear" ]]; then
  SEAR="$SCRIPT_DIR/sear"
elif command -v sear >/dev/null 2>&1; then
  SEAR="$(command -v sear)"
else
  echo "ERROR: sear not found in $SCRIPT_DIR or PATH"; exit 1
fi

if [[ -f "$SCRIPT_DIR/sear_multi" ]]; then
  SEAR_MULTI="$SCRIPT_DIR/sear_multi"
elif command -v sear_multi >/dev/null 2>&1; then
  SEAR_MULTI="$(command -v sear_multi)"
else
  echo "ERROR: sear_multi not found in $SCRIPT_DIR or PATH"; exit 1
fi
chmod +x "$SEAR" "$SEAR_MULTI" 2>/dev/null || true

[[ -f "$CONSENSUS" ]] || { echo "ERROR: Consensus not found: $CONSENSUS"; exit 1; }
[[ -f "$GENOMES_TSV" ]] || { echo "ERROR: Genomes TSV not found: $GENOMES_TSV"; exit 1; }

# ── Select genomes (first 3) ────────────────────────────────────────────────
declare -a SPECIES=()
declare -a GENOME_PATHS=()

while IFS=$'\t' read -r sp gp rest; do
  [[ "$sp" =~ ^[[:space:]]*# ]]  && continue
  [[ "$sp" =~ ^[[:space:]]*$ ]]  && continue
  sp="$(echo "$sp" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  gp="$(echo "$gp" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  [[ -f "$gp" ]] || { echo "WARNING: Genome not found: $gp — skipping $sp"; continue; }
  SPECIES+=("$sp")
  GENOME_PATHS+=("$gp")
  (( ${#SPECIES[@]} >= 3 )) && break
done < "$GENOMES_TSV"

echo "============================================================"
echo "  sear vs sear_multi benchmark"
echo "============================================================"
echo "  Genomes:    ${SPECIES[*]}"
echo "  Queries:    $N_QUERIES (from $CONSENSUS)"
echo "  Genome cap: ${MAX_MB} MB"
echo "============================================================"
echo ""

# ── Prepare test directory ───────────────────────────────────────────────────
BENCHDIR="$(pwd)/benchmark_$(date '+%Y%m%d_%H%M%S')"
mkdir -p "$BENCHDIR"
echo "Working directory: $BENCHDIR"
echo ""

# ── Select N queries ────────────────────────────────────────────────────────
FULL_CONS="$(readlink -f "$CONSENSUS")"
TEST_CONS="$BENCHDIR/test_queries.fa"
seqkit head -n "$N_QUERIES" "$FULL_CONS" > "$TEST_CONS"
ACTUAL_Q=$(grep -c "^>" "$TEST_CONS" || true)
[[ "$ACTUAL_Q" -gt 0 ]] || { echo "ERROR: No sequences selected from consensus"; exit 1; }
echo "Selected $ACTUAL_Q query sequences"
grep "^>" "$TEST_CONS" | sed 's/>/  /'
echo ""

# ── Truncate genomes ────────────────────────────────────────────────────────
declare -a TEST_GENOMES=()
for i in "${!SPECIES[@]}"; do
  sp="${SPECIES[$i]}"
  gp="${GENOME_PATHS[$i]}"
  tg="$BENCHDIR/${sp}_${MAX_MB}mb.fa"

  echo -n "Truncating $sp to ~${MAX_MB}MB ... "
  awk -v max="$MAX_BYTES" '
    /^>/ { if (total >= max) exit; header=$0; next }
    { total += length; if (total <= max + length) { if (header) { print header; header="" }; print } else { exit } }
  ' "$gp" > "$tg"

  # Index
  samtools faidx "$tg" 2>/dev/null || true
  actual_mb=$(awk '{s+=$2} END{printf "%.0f", s/1000000}' "$tg".fai)
  n_seqs=$(wc -l < "$tg".fai)
  echo "${actual_mb} MB, ${n_seqs} sequences"

  TEST_GENOMES+=("$tg")
done
echo ""

# ── Log directory ────────────────────────────────────────────────────────────
LOGDIR="$BENCHDIR/logs"
mkdir -p "$LOGDIR"

# ── Results table ────────────────────────────────────────────────────────────
RESULTS="$BENCHDIR/results.tsv"
printf "Species\tMethod\tWall_s\tTotal_hits\tQueries\n" > "$RESULTS"

compare_beds() {
  local dir_old="$1" dir_new="$2" species="$3"
  local match=0 diff=0 total=0

  for bed_old in "$dir_old"/*.bed; do
    [[ -f "$bed_old" ]] || continue
    bn=$(basename "$bed_old")
    bed_new="$dir_new/$bn"
    total=$((total + 1))

    if [[ ! -f "$bed_new" ]]; then
      echo "    MISSING in new: $bn"
      diff=$((diff + 1))
      continue
    fi

    # Sort both and compare (ignore order differences)
    old_n=$(wc -l < "$bed_old")
    new_n=$(wc -l < "$bed_new")

    if [[ "$old_n" == "$new_n" ]]; then
      # Quick check: sorted content identical?
      if diff <(sort "$bed_old") <(sort "$bed_new") >/dev/null 2>&1; then
        match=$((match + 1))
      else
        echo "    DIFFER: $bn (same lines=$old_n, content differs)"
        diff=$((diff + 1))
      fi
    else
      echo "    DIFFER: $bn (old=$old_n lines, new=$new_n lines)"
      diff=$((diff + 1))
    fi
  done

  echo "  BED comparison: $match/$total identical, $diff different"
}

# ── Run benchmarks ──────────────────────────────────────────────────────────
for i in "${!SPECIES[@]}"; do
  sp="${SPECIES[$i]}"
  tg="${TEST_GENOMES[$i]}"

  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "  $sp"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  # ── OLD approach: split + sear per query ─────────────────────────────
  OLD_DIR="$BENCHDIR/${sp}_old"
  mkdir -p "$OLD_DIR"

  echo ""
  echo "  [OLD] Split consensus + sear per query ..."

  (
    cd "$OLD_DIR"
    # Split consensus
    seqkit split -i "$TEST_CONS" -O .

    CONS_BN=$(basename "$TEST_CONS")
    CONS_NAME="${CONS_BN%.*}"
    CONS_EXT="${CONS_BN##*.}"

    TSTART=$(date +%s%N)

    for query in "${CONS_NAME}".part_*."${CONS_EXT}"; do
      [[ -f "$query" ]] || continue
      echo "    sear: $(basename "$query")"
      "$SEAR" -k "$query" "$tg" 0.8 65 50 >> "$LOGDIR/${sp}_old.log" 2>&1
    done

    TEND=$(date +%s%N)
    ELAPSED=$(awk "BEGIN{printf \"%.1f\", ($TEND - $TSTART)/1000000000}")
    echo "$ELAPSED" > elapsed.txt

    # Clean split files
    rm -f *.2k.bnk *.2k.bnk.fai *.2k.part_* *_s.fai 2>/dev/null || true
  ) || {
    echo "  [OLD] ERROR: sear failed for $sp (see $LOGDIR/${sp}_old.log)"
    echo "0" > "$OLD_DIR/elapsed.txt"
  }

  OLD_SEC=$(cat "$OLD_DIR/elapsed.txt" 2>/dev/null || echo 0)
  OLD_HITS=0
  for bed in "$OLD_DIR"/*.bed; do
    [[ -f "$bed" ]] || continue
    OLD_HITS=$((OLD_HITS + $(wc -l < "$bed")))
  done
  printf "  [OLD] Time: %ss, Total hits: %d\n" "$OLD_SEC" "$OLD_HITS"
  printf "%s\told\t%s\t%d\t%d\n" "$sp" "$OLD_SEC" "$OLD_HITS" "$ACTUAL_Q" >> "$RESULTS"

  # ── NEW approach: sear_multi ─────────────────────────────────────────
  NEW_DIR="$BENCHDIR/${sp}_new"
  mkdir -p "$NEW_DIR"

  echo ""
  echo "  [NEW] sear_multi (all queries in single pass) ..."

  (
    cd "$NEW_DIR"

    TSTART=$(date +%s%N)
    "$SEAR_MULTI" "$TEST_CONS" "$tg" 0.8 65 50 >> "$LOGDIR/${sp}_new.log" 2>&1
    TEND=$(date +%s%N)
    ELAPSED=$(awk "BEGIN{printf \"%.1f\", ($TEND - $TSTART)/1000000000}")
    echo "$ELAPSED" > elapsed.txt
  ) || {
    echo "  [NEW] ERROR: sear_multi failed for $sp (see $LOGDIR/${sp}_new.log)"
    echo "0" > "$NEW_DIR/elapsed.txt"
  }

  NEW_SEC=$(cat "$NEW_DIR/elapsed.txt" 2>/dev/null || echo 0)
  NEW_HITS=0
  for bed in "$NEW_DIR"/*.bed; do
    [[ -f "$bed" ]] || continue
    NEW_HITS=$((NEW_HITS + $(wc -l < "$bed")))
  done
  printf "  [NEW] Time: %ss, Total hits: %d\n" "$NEW_SEC" "$NEW_HITS"
  printf "%s\tnew\t%s\t%d\t%d\n" "$sp" "$NEW_SEC" "$NEW_HITS" "$ACTUAL_Q" >> "$RESULTS"

  # ── Compare outputs ──────────────────────────────────────────────────
  echo ""
  echo "  Comparing BED outputs ..."

  # Normalize both OLD and NEW to common naming: <seqname>.bed
  NORM_OLD="$BENCHDIR/${sp}_old_norm"
  NORM_NEW="$BENCHDIR/${sp}_new_norm"
  mkdir -p "$NORM_OLD" "$NORM_NEW"

  # OLD: Boa-test_queries.part_Heno.bed → Heno.bed
  for bed in "$OLD_DIR"/*.bed; do
    [[ -f "$bed" ]] || continue
    bn=$(basename "$bed")
    seqname="${bn%.bed}"
    seqname="${seqname#*-}"    # strip taxname prefix (e.g. "Boa-")
    if [[ "$seqname" == *".part_"* ]]; then
      seqname="${seqname##*.part_}"    # strip consensus .part_ prefix
    fi
    cp "$bed" "$NORM_OLD/${seqname}.bed"
  done

  # NEW: Boa-Heno.bed → Heno.bed
  for bed in "$NEW_DIR"/*.bed; do
    [[ -f "$bed" ]] || continue
    bn=$(basename "$bed")
    seqname="${bn%.bed}"
    seqname="${seqname#*-}"    # strip taxname prefix
    cp "$bed" "$NORM_NEW/${seqname}.bed"
  done

  compare_beds "$NORM_OLD" "$NORM_NEW" "$sp"

  speedup=$(awk -v o="$OLD_SEC" -v n="$NEW_SEC" 'BEGIN{ if (n+0 > 0) printf "%.1f", o/n; else print "N/A" }')
  [[ "$speedup" != "N/A" ]] && echo "  Speedup: ${speedup}×"

  echo ""
done

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  SUMMARY"
echo "============================================================"
echo ""
column -t -s $'\t' "$RESULTS" 2>/dev/null || cat "$RESULTS"
echo ""
echo "Results saved to: $RESULTS"
echo "Logs saved to: $LOGDIR/"
echo "Working directory: $BENCHDIR"
echo "============================================================"
