#!/usr/bin/env bash
# ==========================================================
# step2_assign.sh - Assign SINE sequences to consensus subfamilies
# ==========================================================
# Strict assignment based on 10-cycle voting + bitscore thresholds
# No alignment steps - fast assignment only
set -euo pipefail

THREADS_DEFAULT="$(nproc 2>/dev/null || echo 1)"
SPLIT_SIZE=20000
MIN_REL_BITSCORE=0.45

TMPDIR=""
OUTDIR=""
CONSENSUS_FILE=""
SINE_FILE=""
THREADS=""
INCREMENTAL=false
OLD_ASSIGNMENT=""
NEW_BEDS_FILE=""

die(){ printf '[ERROR] %s\n' "$*" >&2; exit 1; }
log(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
abs_path(){ local p="$1"; [[ "$p" = /* ]] && printf '%s\n' "$p" || printf '%s\n' "$PWD/$p"; }
cleanup(){ [[ -n "${TMPDIR:-}" && -d "$TMPDIR" ]] && rm -rf "$TMPDIR"; }
trap cleanup EXIT

check_deps(){
  local req=(seqkit ssearch36 awk sed grep sort)
  local missing=()
  for t in "${req[@]}"; do
    command -v "$t" >/dev/null 2>&1 || missing+=("$t")
  done
  [[ ${#missing[@]} -eq 0 ]] || die "Missing tools: ${missing[*]}"
}

show_help(){
  cat <<'EOF'
Assign SINE sequences to consensus subfamilies using strict voting + thresholds.

Usage:
  step2_assign.sh <consensus.fa> <sines.fa> [threads] [output_dir]

Arguments:
  <consensus.fa>   FASTA file of consensus sequences (subfamily representatives)
  <sines.fa>       FASTA file of SINE sequences to assign
  [threads]        Optional: number of CPU threads (default: auto-detect)
  [output_dir]     Optional: output directory name (default: step2_output)

Output:
  assigned.fasta       Sequences assigned to subfamilies (headers: >seqID|subfamily|bitscore)
  unassigned.fasta     Sequences that failed assignment criteria
  assignment_stats.tsv Summary statistics per subfamily
  assignment_full.tsv  Detailed assignment info for all sequences

Assignment Rules:
  1) 10/10 unanimous votes across ssearch36 cycles (same consensus wins all 10)
  2) Bitscore >= threshold (calculated as MIN_REL_BITSCORE . TopN for that subfamily)

Examples:
  step2_assign.sh consensuses.fa my_sines.fa
  step2_assign.sh consensuses.fa my_sines.fa 32 my_results
EOF
}

parse_args(){
  if [[ $# -eq 1 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
    show_help
    exit 0
  fi

  # Parse optional flags before positional args
  while [[ $# -gt 0 && "$1" == --* ]]; do
    case "$1" in
      --incremental) INCREMENTAL=true; shift ;;
      --old-assignment)
        shift; [[ $# -gt 0 ]] || die "--old-assignment requires a file argument"
        OLD_ASSIGNMENT="$(abs_path "$1")"; shift ;;
      --new-beds)
        shift; [[ $# -gt 0 ]] || die "--new-beds requires a file argument"
        NEW_BEDS_FILE="$(abs_path "$1")"; shift ;;
      --help) show_help; exit 0 ;;
      *) die "Unknown option: $1" ;;
    esac
  done

  if [[ "$INCREMENTAL" == true ]]; then
    [[ -s "$OLD_ASSIGNMENT" ]] || die "Incremental mode requires --old-assignment <file>"
    [[ -s "$NEW_BEDS_FILE" ]] || log "WARNING: --new-beds file is empty or missing; all old assigned will be carried forward"
    command -v bedtools >/dev/null 2>&1 || die "Incremental mode requires bedtools"
  fi

  if [[ $# -lt 2 ]]; then
    echo "[ERROR] Expected at least 2 arguments. Run with --help for usage." >&2
    exit 1
  fi

  [[ -f "$1" ]] || die "Consensus file not found: $1"
  [[ -f "$2" ]] || die "SINE file not found: $2"

  CONSENSUS_FILE="$(abs_path "$1")"
  SINE_FILE="$(abs_path "$2")"
  THREADS="${3:-$THREADS_DEFAULT}"
  OUTDIR="${4:-step2_output}"

  # Create unique output directory if it exists
  if [[ -e "$OUTDIR" ]]; then
    local n=1
    local base="$OUTDIR"
    while [[ -e "${base}_run${n}" ]]; do
      n=$((n+1))
    done
    OUTDIR="${base}_run${n}"
  fi

  mkdir -p "$OUTDIR"
  OUTDIR="$(cd "$OUTDIR" && pwd)"

  log "Consensus file: $CONSENSUS_FILE"
  log "SINE file: $SINE_FILE"
  log "Threads: $THREADS"
  log "Output directory: $OUTDIR"
}

count_fasta(){
  awk 'BEGIN{c=0} /^>/{c++} END{printf("%d\n", c)}' "$1"
}

run_assignment(){
  log "=========================================="
  log " Assignment: SINE loci . Subfamilies"
  log "=========================================="

  TMPDIR="$(mktemp -d -t step2_assign_XXXX)"
  export LC_ALL=C

  # Clean headers (keep only first word after >)
  log "Preparing input files..."
  awk '/^>/{split($1,a," ");print a[1];next}{print}' "$CONSENSUS_FILE" > "$TMPDIR/consensus_clean.fa"
  seqkit seq -g "$SINE_FILE" > "$TMPDIR/sines_clean.fa" 2>"$OUTDIR/seqkit.log"

  local consensus_count sine_count
  consensus_count="$(count_fasta "$TMPDIR/consensus_clean.fa")"
  sine_count="$(count_fasta "$TMPDIR/sines_clean.fa")"

  log "Input: $sine_count SINE sequences"
  log "References: $consensus_count consensus subfamilies"

  # Split SINE file if needed
  rm -rf "$TMPDIR/splits"
  mkdir -p "$TMPDIR/splits"

  if [[ $sine_count -le $SPLIT_SIZE ]]; then
    cp "$TMPDIR/sines_clean.fa" "$TMPDIR/splits/part_0001.fasta"
    log "No splitting needed (. $SPLIT_SIZE sequences)"
  else
    log "Splitting into chunks of $SPLIT_SIZE sequences..."
    seqkit split -s "$SPLIT_SIZE" -O "$TMPDIR/splits_raw" "$TMPDIR/sines_clean.fa" 2>>"$OUTDIR/seqkit.log"

    # Collect split files (handles different seqkit versions)
    mapfile -t split_files < <(
      find "$TMPDIR/splits_raw" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fas" \) | sort
    )

    [[ ${#split_files[@]} -gt 0 ]] || die "seqkit split produced no output files"

    local i=1
    for f in "${split_files[@]}"; do
      printf -v tag "%04d" "$i"
      cp "$f" "$TMPDIR/splits/part_${tag}.fasta"
      i=$((i+1))
    done
    rm -rf "$TMPDIR/splits_raw"
  fi

  # Run 10 voting cycles
  log "Running 10 ssearch36 voting cycles..."
  : > "$TMPDIR/votes.tsv"
  : > "$OUTDIR/ssearch.log"

  for cycle in $(seq 1 10); do
    printf "  Cycle %2d/10\r" "$cycle" >&2
    : > "$TMPDIR/raw_${cycle}.m8"

    for part in "$TMPDIR"/splits/*.fasta; do
      ssearch36 -g -3 -T "$THREADS" -Q -n -z 11 -E 2 -w 95 -W 70 -m 8 \
        "$TMPDIR/consensus_clean.fa" "$part" \
        >> "$TMPDIR/raw_${cycle}.m8" 2>>"$OUTDIR/ssearch.log" || true
    done

    # Extract best hit per sequence for this cycle
    awk '{
      seq=$2
      cons=$1
      bs=$12+0
      if(!(seq in best_bs) || bs > best_bs[seq]){
        best_bs[seq] = bs
        best_cons[seq] = cons
      }
    }
    END{
      for(seq in best_cons){
        printf "%s\t%s\t%d\n", seq, best_cons[seq], int(best_bs[seq]+0.5)
      }
    }' "$TMPDIR/raw_${cycle}.m8" >> "$TMPDIR/votes.tsv"
  done
  echo >&2 "  Cycle 10/10 done."

  # Aggregate votes
  log "Aggregating votes and calculating thresholds..."
  awk 'BEGIN{FS=OFS="\t"}{
    key = $1 FS $2
    vote_count[key]++
    bitscore_sum[key] += $3
  }
  END{
    for(key in vote_count){
      split(key, a, FS)
      print a[1], a[2], vote_count[key], bitscore_sum[key]
    }
  }' "$TMPDIR/votes.tsv" | sort -k1,1 -k3,3nr -k2,2 > "$TMPDIR/vote_summary.tsv"

  # Save best-vote info for ALL sequences (for unassigned tracking)
  log "Recording per-sequence best vote info (all sequences)..."
  awk 'BEGIN{FS=OFS="\t"}{
    seq = $1
    cons = $2
    votes = $3 + 0
    bs_sum = $4 + 0

    if(!(seq in best_votes) || votes > best_votes[seq] ||
       (votes == best_votes[seq] && bs_sum > best_bs[seq])){
      best_votes[seq] = votes
      best_bs[seq] = bs_sum
      best_cons[seq] = cons
    }
  }
  END{
    for(seq in best_cons){
      printf "%s\t%s\t%d\t%d\n", seq, best_cons[seq], best_votes[seq], int(best_bs[seq]+0.5)
    }
  }' "$TMPDIR/vote_summary.tsv" | sort -k1,1 > "$OUTDIR/all_votes.tsv"

  # Extract sequences with 10/10 unanimous votes
  awk 'BEGIN{FS=OFS="\t"}{
    seq = $1
    cons = $2
    votes = $3 + 0
    bs_sum = $4 + 0

    if(votes == 10){
      if(!(seq in best_bs) || bs_sum > best_bs[seq]){
        best_bs[seq] = bs_sum
        best_cons[seq] = cons
      }
    }
  }
  END{
    for(seq in best_cons){
      printf "%s\t%s\t%d\t10\n", seq, best_cons[seq], int(best_bs[seq]+0.5)
    }
  }' "$TMPDIR/vote_summary.tsv" | sort -k2,2 -k3,3nr > "$TMPDIR/unanimous.tsv"

  local unanimous_count
  unanimous_count="$(awk 'END{print NR}' "$TMPDIR/unanimous.tsv")"
  log "Sequences with 10/10 unanimous votes: $unanimous_count"

  # Calculate thresholds per consensus
  log "Calculating bitscore thresholds per subfamily..."
  awk '{print $2}' "$TMPDIR/unanimous.tsv" | sort -u > "$TMPDIR/consensus_list.txt"
  : > "$TMPDIR/thresholds.tsv"

  while IFS= read -r cons; do
    [[ -n "$cons" ]] || continue

    # Get bitscores for this consensus, sorted descending
    awk -v c="$cons" '$2==c{print $3}' "$TMPDIR/unanimous.tsv" | sort -nr > "$TMPDIR/${cons}.bitscores"

    local cnt
    cnt="$(awk 'END{print NR}' "$TMPDIR/${cons}.bitscores")"

    if [[ ${cnt:-0} -lt 1 ]]; then
      printf "%s\t%d\t%d\n" "$cons" 0 0 >> "$TMPDIR/thresholds.tsv"
      continue
    fi

    # Use top-N (N = min(10, count)) to set threshold
    local N=$(( cnt < 10 ? cnt : 10 ))
    local nth_bitscore
    nth_bitscore="$(awk -v n="$N" 'NR==n{print; exit}' "$TMPDIR/${cons}.bitscores")"

    # Threshold = MIN_REL_BITSCORE * N-th best bitscore
    local threshold
    threshold="$(awk -v x="$nth_bitscore" -v m="$MIN_REL_BITSCORE" 'BEGIN{printf "%d", int(x*m+0.5)}')"

    printf "%s\t%d\t%d\n" "$cons" "$nth_bitscore" "$threshold" >> "$TMPDIR/thresholds.tsv"
  done < "$TMPDIR/consensus_list.txt"

  # Apply thresholds to determine final assignments
  log "Applying thresholds to filter assignments..."
  awk 'BEGIN{FS=OFS="\t"}
  NR==FNR{
    threshold[$1] = $3
    topN[$1] = $2
    next
  }
  {
    seq = $1
    cons = $2
    bs = $3 + 0
    votes = $4

    if(bs >= threshold[cons]){
      status = "assigned"
    } else {
      status = "rejected_low_bitscore"
    }

    printf "%s\t%s\t%d\t%d\t%s\t%d\n", seq, cons, bs, votes, status, threshold[cons]
  }' "$TMPDIR/thresholds.tsv" "$TMPDIR/unanimous.tsv" > "$TMPDIR/assignment_details.tsv"

  # Extract assigned sequence IDs
  awk '$5=="assigned"{print $1}' "$TMPDIR/assignment_details.tsv" > "$TMPDIR/assigned.ids"
  local assigned_count
  assigned_count="$(awk 'END{print NR}' "$TMPDIR/assigned.ids")"

  log "Sequences passing all criteria: $assigned_count"

  # Create assigned.fasta with renamed headers: >seqID|subfamily|bitscore
  log "Writing assigned.fasta..."
  if [[ $assigned_count -gt 0 ]]; then
    awk 'BEGIN{FS="\t"}
    NR==FNR{
      cons[$1] = $2
      bs[$1] = $3
      next
    }
    /^>/{
      name = substr($0, 2)
      split(name, a, " ")
      seqid = a[1]
      if(seqid in cons){
        printf ">%s|%s|%d\n", seqid, cons[seqid], bs[seqid]
        in_assigned = 1
      } else {
        in_assigned = 0
      }
      next
    }
    {
      if(in_assigned) print
    }' "$TMPDIR/assignment_details.tsv" "$TMPDIR/sines_clean.fa" > "$OUTDIR/assigned.fasta"
  else
    : > "$OUTDIR/assigned.fasta"
  fi

  # Create unassigned.fasta
  log "Writing unassigned.fasta..."
  local unassigned_count
  if [[ $assigned_count -gt 0 ]]; then
    sort -u "$TMPDIR/assigned.ids" > "$TMPDIR/assigned_sorted.ids"
    seqkit grep -v -f "$TMPDIR/assigned_sorted.ids" "$TMPDIR/sines_clean.fa" \
      > "$OUTDIR/unassigned.fasta" 2>>"$OUTDIR/seqkit.log"
  else
    cp "$TMPDIR/sines_clean.fa" "$OUTDIR/unassigned.fasta"
  fi
  unassigned_count="$(count_fasta "$OUTDIR/unassigned.fasta")"

  # Split assigned sequences into per-subfamily FASTA files
  log "Splitting assigned sequences into per-subfamily FASTA files..."
  mkdir -p "$OUTDIR/subfamilies"

  if [[ $assigned_count -gt 0 ]]; then
    # Get list of unique subfamilies
    awk '$5=="assigned"{print $2}' "$TMPDIR/assignment_details.tsv" | sort -u > "$TMPDIR/subfamily_list.txt"

    while IFS= read -r subfamily; do
      [[ -z "$subfamily" ]] && continue

      # Get sequence IDs for this subfamily
      awk -v sf="$subfamily" '$5=="assigned" && $2==sf{print $1}' "$TMPDIR/assignment_details.tsv" \
        > "$TMPDIR/${subfamily}.ids"

      # Extract sequences with renamed headers
      awk 'BEGIN{FS="\t"}
      NR==FNR{
        cons[$1] = $2
        bs[$1] = $3
        next
      }
      /^>/{
        name = substr($0, 2)
        split(name, a, " ")
        seqid = a[1]
        if(seqid in cons && cons[seqid] == SUBFAMILY){
          printf ">%s|%s|%d\n", seqid, cons[seqid], bs[seqid]
          in_subfamily = 1
        } else {
          in_subfamily = 0
        }
        next
      }
      {
        if(in_subfamily) print
      }' SUBFAMILY="$subfamily" "$TMPDIR/assignment_details.tsv" "$TMPDIR/sines_clean.fa" \
        > "$OUTDIR/subfamilies/${subfamily}.fasta"

      local count
      count="$(count_fasta "$OUTDIR/subfamilies/${subfamily}.fasta")"
      log "  ${subfamily}: $count sequences"
    done < "$TMPDIR/subfamily_list.txt"
  fi

  # Generate statistics
  log "Generating statistics..."

  # Per-subfamily stats (FIXED: sort data only, then prepend header)
  {
    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR && $5=="assigned"{
      count[$2]++
      next
    }
    NR!=FNR{
      cons = $1
      topN = $2
      thr = $3
      c = (cons in count ? count[cons] : 0)
      printf "%s\t%d\t%d\t%d\n", cons, c, topN, thr
    }' "$TMPDIR/assignment_details.tsv" "$TMPDIR/thresholds.tsv" \
    | sort -k2,2nr -k1,1
  } > "$OUTDIR/assignment_stats.data.tsv"

  {
    echo -e "Subfamily\tAssigned\tTopN_Bitscore\tThreshold"
    cat "$OUTDIR/assignment_stats.data.tsv"
  } > "$OUTDIR/assignment_stats.tsv"

  rm -f "$OUTDIR/assignment_stats.data.tsv"

  # Full assignment details: includes unanimous (assigned/rejected) + non-unanimous
  {
    echo -e "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold"
    # 1) unanimous sequences (assigned or rejected_low_bitscore)
    cat "$TMPDIR/assignment_details.tsv"
    # 2) non-unanimous sequences: best vote from all_votes.tsv minus unanimous
    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR{ seen[$1]=1; next }
    !($1 in seen){
      printf "%s\t%s\t%d\t%d\t%s\t%s\n", $1, $2, $4, $3, "no_unanimous", "NA"
    }' "$TMPDIR/unanimous.tsv" "$OUTDIR/all_votes.tsv"
  } > "$OUTDIR/assignment_full.tsv"

  # Summary report
  {
    echo "=========================================="
    echo "  SINE Assignment Summary"
    echo "=========================================="
    echo "Date: $(date)"
    echo
    echo "Input Files:"
    echo "  Consensus: $CONSENSUS_FILE ($consensus_count subfamilies)"
    echo "  SINEs: $SINE_FILE ($sine_count sequences)"
    echo
    echo "Assignment Criteria:"
    echo "  1) 10/10 unanimous votes (same consensus wins all cycles)"
    echo "  2) Bitscore >= ${MIN_REL_BITSCORE} . TopN (N=min(10,count))"
    echo
    echo "Results:"
    echo "  Total sequences: $sine_count"
    echo "  Unanimous 10/10: $unanimous_count"
    echo "  Assigned (passed threshold): $assigned_count ($(awk -v a=$assigned_count -v t=$sine_count 'BEGIN{printf "%.1f%%", 100*a/t}'))"
    echo "  Unassigned: $unassigned_count ($(awk -v u=$unassigned_count -v t=$sine_count 'BEGIN{printf "%.1f%%", 100*u/t}'))"
    echo
    echo "Assignments per Subfamily:"
    awk -F'\t' 'NR>1{printf "  %-20s %6d\n", $1, $2}' "$OUTDIR/assignment_stats.tsv"
    echo
    echo "Output Files:"
    echo "  assigned.fasta         - All assigned sequences with headers: >seqID|subfamily|bitscore"
    echo "  unassigned.fasta       - Sequences that failed assignment"
    echo "  assignment_stats.tsv   - Statistics per subfamily"
    echo "  assignment_full.tsv    - Detailed assignment info for all sequences (incl. non-unanimous)"
    echo "  all_votes.tsv          - Per-sequence best vote (all sequences)"
    echo "  subfamilies/           - Directory with per-subfamily FASTA files"
    echo
  } > "$OUTDIR/summary.txt"

  cat "$OUTDIR/summary.txt"
  log "Done! Output written to: $OUTDIR"
}

###############################################################################
# Incremental assignment: carry forward old firm assignments, re-evaluate only
# sequences that overlap new consensus BED or were previously unassigned/new.
###############################################################################
run_incremental_assignment(){
  log "=========================================="
  log " Incremental Assignment (carry-forward mode)"
  log "=========================================="

  TMPDIR="$(mktemp -d -t step2_incr_XXXX)"
  export LC_ALL=C

  log "Preparing input files..."
  awk '/^>/{split($1,a," ");print a[1];next}{print}' "$CONSENSUS_FILE" > "$TMPDIR/consensus_clean.fa"
  seqkit seq -g "$SINE_FILE" > "$TMPDIR/sines_clean.fa" 2>"$OUTDIR/seqkit.log"

  local consensus_count sine_count
  consensus_count="$(count_fasta "$TMPDIR/consensus_clean.fa")"
  sine_count="$(count_fasta "$TMPDIR/sines_clean.fa")"

  log "Input: $sine_count SINE sequences, $consensus_count consensus subfamilies"

  #---------------------------------------------------------------------------
  # Partition into Tier 1 (carry forward) and Tier 2 (re-evaluate)
  #---------------------------------------------------------------------------
  log "Partitioning sequences..."

  # Tier 1 candidates: firmly assigned (10/10 votes) in old run
  awk -F'\t' 'NR>1 && $5=="assigned" && $4+0==10{print $1}' "$OLD_ASSIGNMENT" \
    | sort -u > "$TMPDIR/tier1_candidates.ids"

  # All current sequence IDs
  awk '/^>/{sub(/^>/,""); print $1}' "$TMPDIR/sines_clean.fa" \
    | sort -u > "$TMPDIR/all_current.ids"

  # Tier 1 candidates present in current extraction
  comm -12 "$TMPDIR/tier1_candidates.ids" "$TMPDIR/all_current.ids" \
    > "$TMPDIR/tier1_present.ids"

  # Convert Tier 1 IDs to BED for overlap check (ID format: chr:start-end(strand))
  awk '{
    id=$1; chr=id; sub(/:.*/,"",chr)
    rng=id; sub(/^[^:]+:/,"",rng); sub(/\([+-]\).*/,"",rng)
    split(rng,xy,"-"); s=xy[1]+0-1; e=xy[2]+0
    if(s<0) s=0
    print chr"\t"s"\t"e"\t"id
  }' "$TMPDIR/tier1_present.ids" | sort -k1,1 -k2,2n > "$TMPDIR/tier1.bed"

  # Check overlap with new consensus BED
  if [[ -s "$NEW_BEDS_FILE" && -s "$TMPDIR/tier1.bed" ]]; then
    bedtools intersect -u -a "$TMPDIR/tier1.bed" -b "$NEW_BEDS_FILE" \
      | awk '{print $4}' | sort -u > "$TMPDIR/tier1_overlapping.ids"
  else
    : > "$TMPDIR/tier1_overlapping.ids"
  fi

  # Final Tier 1: present AND NOT overlapping with new BED
  comm -23 "$TMPDIR/tier1_present.ids" "$TMPDIR/tier1_overlapping.ids" \
    > "$TMPDIR/tier1_final.ids"

  # Tier 2: everything NOT in Tier 1 final
  comm -23 "$TMPDIR/all_current.ids" "$TMPDIR/tier1_final.ids" \
    > "$TMPDIR/tier2.ids"

  local tier1_n tier2_n overlap_n
  tier1_n=$(wc -l < "$TMPDIR/tier1_final.ids")
  tier2_n=$(wc -l < "$TMPDIR/tier2.ids")
  overlap_n=$(wc -l < "$TMPDIR/tier1_overlapping.ids")

  log "Tier 1 (carry forward): $tier1_n | Tier 2 (re-evaluate): $tier2_n | Overlap: $overlap_n"

  # Extract Tier 1 old assignment data (unanimous format: seq cons bs 10)
  awk -F'\t' 'NR==FNR{want[$1]=1; next}
    ($1 in want) && $5=="assigned" && $4+0==10{
      printf "%s\t%s\t%d\t10\n", $1, $2, $3
    }' "$TMPDIR/tier1_final.ids" "$OLD_ASSIGNMENT" > "$TMPDIR/tier1_unanimous.tsv"

  #---------------------------------------------------------------------------
  # Run 10-cycle ssearch36 on Tier 2 only
  #---------------------------------------------------------------------------
  : > "$TMPDIR/tier2_unanimous.tsv"
  : > "$TMPDIR/tier2_all_votes.tsv"
  : > "$OUTDIR/ssearch.log"

  if [[ $tier2_n -gt 0 ]]; then
    seqkit grep -f "$TMPDIR/tier2.ids" "$TMPDIR/sines_clean.fa" \
      > "$TMPDIR/tier2_sines.fa" 2>>"$OUTDIR/seqkit.log"

    local tier2_actual
    tier2_actual="$(count_fasta "$TMPDIR/tier2_sines.fa")"
    log "Running 10 ssearch36 cycles on $tier2_actual Tier 2 sequences..."

    rm -rf "$TMPDIR/splits"; mkdir -p "$TMPDIR/splits"
    if [[ $tier2_actual -le $SPLIT_SIZE ]]; then
      cp "$TMPDIR/tier2_sines.fa" "$TMPDIR/splits/part_0001.fasta"
    else
      seqkit split -s "$SPLIT_SIZE" -O "$TMPDIR/splits_raw" "$TMPDIR/tier2_sines.fa" \
        2>>"$OUTDIR/seqkit.log"
      mapfile -t split_files < <(
        find "$TMPDIR/splits_raw" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fas" \) | sort
      )
      [[ ${#split_files[@]} -gt 0 ]] || die "seqkit split produced no output"
      local i=1
      for f in "${split_files[@]}"; do
        printf -v tag "%04d" "$i"
        cp "$f" "$TMPDIR/splits/part_${tag}.fasta"
        i=$((i+1))
      done
      rm -rf "$TMPDIR/splits_raw"
    fi

    : > "$TMPDIR/votes.tsv"
    for cycle in $(seq 1 10); do
      printf "  Cycle %2d/10\r" "$cycle" >&2
      : > "$TMPDIR/raw_${cycle}.m8"
      for part in "$TMPDIR"/splits/*.fasta; do
        ssearch36 -g -3 -T "$THREADS" -Q -n -z 11 -E 2 -w 95 -W 70 -m 8 \
          "$TMPDIR/consensus_clean.fa" "$part" \
          >> "$TMPDIR/raw_${cycle}.m8" 2>>"$OUTDIR/ssearch.log" || true
      done
      awk '{
        seq=$2; cons=$1; bs=$12+0
        if(!(seq in best_bs) || bs > best_bs[seq]){
          best_bs[seq] = bs; best_cons[seq] = cons
        }
      }
      END{
        for(seq in best_cons)
          printf "%s\t%s\t%d\n", seq, best_cons[seq], int(best_bs[seq]+0.5)
      }' "$TMPDIR/raw_${cycle}.m8" >> "$TMPDIR/votes.tsv"
    done
    echo >&2 "  Cycle 10/10 done."

    # Aggregate Tier 2 votes
    awk 'BEGIN{FS=OFS="\t"}{
      key = $1 FS $2; vote_count[key]++; bitscore_sum[key] += $3
    }
    END{
      for(key in vote_count){
        split(key, a, FS)
        print a[1], a[2], vote_count[key], bitscore_sum[key]
      }
    }' "$TMPDIR/votes.tsv" | sort -k1,1 -k3,3nr -k2,2 > "$TMPDIR/tier2_vote_summary.tsv"

    # Tier 2 all_votes (best vote per sequence)
    awk 'BEGIN{FS=OFS="\t"}{
      seq=$1; cons=$2; votes=$3+0; bs_sum=$4+0
      if(!(seq in best_votes) || votes > best_votes[seq] ||
         (votes == best_votes[seq] && bs_sum > best_bs[seq])){
        best_votes[seq] = votes
        best_bs[seq] = bs_sum
        best_cons[seq] = cons
      }
    }
    END{
      for(seq in best_cons)
        printf "%s\t%s\t%d\t%d\n", seq, best_cons[seq], best_votes[seq], int(best_bs[seq]+0.5)
    }' "$TMPDIR/tier2_vote_summary.tsv" | sort -k1,1 > "$TMPDIR/tier2_all_votes.tsv"

    # Tier 2 unanimous (10/10 only)
    awk 'BEGIN{FS=OFS="\t"}{
      seq=$1; cons=$2; votes=$3+0; bs_sum=$4+0
      if(votes == 10){
        if(!(seq in best_bs) || bs_sum > best_bs[seq]){
          best_bs[seq] = bs_sum; best_cons[seq] = cons
        }
      }
    }
    END{
      for(seq in best_cons)
        printf "%s\t%s\t%d\t10\n", seq, best_cons[seq], int(best_bs[seq]+0.5)
    }' "$TMPDIR/tier2_vote_summary.tsv" | sort -k2,2 -k3,3nr > "$TMPDIR/tier2_unanimous.tsv"
  fi

  #---------------------------------------------------------------------------
  # Merge Tier 1 + Tier 2
  #---------------------------------------------------------------------------
  cat "$TMPDIR/tier1_unanimous.tsv" "$TMPDIR/tier2_unanimous.tsv" \
    | sort -k2,2 -k3,3nr > "$TMPDIR/unanimous.tsv"

  local unanimous_count
  unanimous_count="$(awk 'END{print NR}' "$TMPDIR/unanimous.tsv")"
  log "Combined unanimous: $unanimous_count (Tier1: $(wc -l < "$TMPDIR/tier1_unanimous.tsv"), Tier2: $(wc -l < "$TMPDIR/tier2_unanimous.tsv"))"

  # Combined all_votes: Tier 1 (from old assignment) + Tier 2 (fresh)
  awk -F'\t' 'NR==FNR{want[$1]=1; next}
    ($1 in want){
      printf "%s\t%s\t%d\t%d\n", $1, $2, $4, $3
    }' "$TMPDIR/tier1_final.ids" "$OLD_ASSIGNMENT" > "$TMPDIR/tier1_all_votes.tsv"

  cat "$TMPDIR/tier1_all_votes.tsv" "$TMPDIR/tier2_all_votes.tsv" \
    | sort -k1,1 > "$OUTDIR/all_votes.tsv"

  #---------------------------------------------------------------------------
  # Calculate thresholds from combined unanimous
  #---------------------------------------------------------------------------
  log "Calculating bitscore thresholds per subfamily..."
  awk '{print $2}' "$TMPDIR/unanimous.tsv" | sort -u > "$TMPDIR/consensus_list.txt"
  : > "$TMPDIR/thresholds.tsv"

  while IFS= read -r cons; do
    [[ -n "$cons" ]] || continue
    awk -v c="$cons" '$2==c{print $3}' "$TMPDIR/unanimous.tsv" | sort -nr \
      > "$TMPDIR/${cons}.bitscores"
    local cnt
    cnt="$(awk 'END{print NR}' "$TMPDIR/${cons}.bitscores")"
    if [[ ${cnt:-0} -lt 1 ]]; then
      printf "%s\t%d\t%d\n" "$cons" 0 0 >> "$TMPDIR/thresholds.tsv"
      continue
    fi
    local N=$(( cnt < 10 ? cnt : 10 ))
    local nth_bitscore
    nth_bitscore="$(awk -v n="$N" 'NR==n{print; exit}' "$TMPDIR/${cons}.bitscores")"
    local threshold
    threshold="$(awk -v x="$nth_bitscore" -v m="$MIN_REL_BITSCORE" 'BEGIN{printf "%d", int(x*m+0.5)}')"
    printf "%s\t%d\t%d\n" "$cons" "$nth_bitscore" "$threshold" >> "$TMPDIR/thresholds.tsv"
  done < "$TMPDIR/consensus_list.txt"

  # Apply thresholds
  log "Applying thresholds to filter assignments..."
  awk 'BEGIN{FS=OFS="\t"}
  NR==FNR{
    threshold[$1] = $3; topN[$1] = $2; next
  }
  {
    seq = $1; cons = $2; bs = $3 + 0; votes = $4
    if(bs >= threshold[cons]) status = "assigned"
    else status = "rejected_low_bitscore"
    printf "%s\t%s\t%d\t%d\t%s\t%d\n", seq, cons, bs, votes, status, threshold[cons]
  }' "$TMPDIR/thresholds.tsv" "$TMPDIR/unanimous.tsv" > "$TMPDIR/assignment_details.tsv"

  awk '$5=="assigned"{print $1}' "$TMPDIR/assignment_details.tsv" > "$TMPDIR/assigned.ids"
  local assigned_count
  assigned_count="$(awk 'END{print NR}' "$TMPDIR/assigned.ids")"
  log "Sequences passing all criteria: $assigned_count"

  #---------------------------------------------------------------------------
  # Generate output files (same format as full mode)
  #---------------------------------------------------------------------------
  log "Writing assigned.fasta..."
  if [[ $assigned_count -gt 0 ]]; then
    awk 'BEGIN{FS="\t"}
    NR==FNR{
      cons[$1] = $2; bs[$1] = $3; next
    }
    /^>/{
      name = substr($0, 2); split(name, a, " "); seqid = a[1]
      if(seqid in cons){
        printf ">%s|%s|%d\n", seqid, cons[seqid], bs[seqid]; in_assigned = 1
      } else { in_assigned = 0 }
      next
    }
    { if(in_assigned) print
    }' "$TMPDIR/assignment_details.tsv" "$TMPDIR/sines_clean.fa" > "$OUTDIR/assigned.fasta"
  else
    : > "$OUTDIR/assigned.fasta"
  fi

  log "Writing unassigned.fasta..."
  local unassigned_count
  if [[ $assigned_count -gt 0 ]]; then
    sort -u "$TMPDIR/assigned.ids" > "$TMPDIR/assigned_sorted.ids"
    seqkit grep -v -f "$TMPDIR/assigned_sorted.ids" "$TMPDIR/sines_clean.fa" \
      > "$OUTDIR/unassigned.fasta" 2>>"$OUTDIR/seqkit.log"
  else
    cp "$TMPDIR/sines_clean.fa" "$OUTDIR/unassigned.fasta"
  fi
  unassigned_count="$(count_fasta "$OUTDIR/unassigned.fasta")"

  log "Splitting assigned sequences into per-subfamily FASTA files..."
  mkdir -p "$OUTDIR/subfamilies"
  if [[ $assigned_count -gt 0 ]]; then
    awk '$5=="assigned"{print $2}' "$TMPDIR/assignment_details.tsv" | sort -u \
      > "$TMPDIR/subfamily_list.txt"
    while IFS= read -r subfamily; do
      [[ -z "$subfamily" ]] && continue
      awk 'BEGIN{FS="\t"}
      NR==FNR{
        cons[$1] = $2; bs[$1] = $3; next
      }
      /^>/{
        name = substr($0, 2); split(name, a, " "); seqid = a[1]
        if(seqid in cons && cons[seqid] == SUBFAMILY){
          printf ">%s|%s|%d\n", seqid, cons[seqid], bs[seqid]; in_subfamily = 1
        } else { in_subfamily = 0 }
        next
      }
      { if(in_subfamily) print
      }' SUBFAMILY="$subfamily" "$TMPDIR/assignment_details.tsv" "$TMPDIR/sines_clean.fa" \
        > "$OUTDIR/subfamilies/${subfamily}.fasta"
      local count
      count="$(count_fasta "$OUTDIR/subfamilies/${subfamily}.fasta")"
      log "  ${subfamily}: $count sequences"
    done < "$TMPDIR/subfamily_list.txt"
  fi

  # Generate statistics
  log "Generating statistics..."
  {
    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR && $5=="assigned"{
      count[$2]++; next
    }
    NR!=FNR{
      cons = $1; topN = $2; thr = $3
      c = (cons in count ? count[cons] : 0)
      printf "%s\t%d\t%d\t%d\n", cons, c, topN, thr
    }' "$TMPDIR/assignment_details.tsv" "$TMPDIR/thresholds.tsv" \
    | sort -k2,2nr -k1,1
  } > "$OUTDIR/assignment_stats.data.tsv"

  {
    echo -e "Subfamily\tAssigned\tTopN_Bitscore\tThreshold"
    cat "$OUTDIR/assignment_stats.data.tsv"
  } > "$OUTDIR/assignment_stats.tsv"
  rm -f "$OUTDIR/assignment_stats.data.tsv"

  # Full assignment details
  {
    echo -e "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold"
    cat "$TMPDIR/assignment_details.tsv"
    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR{ seen[$1]=1; next }
    !($1 in seen){
      printf "%s\t%s\t%d\t%d\t%s\t%s\n", $1, $2, $4, $3, "no_unanimous", "NA"
    }' "$TMPDIR/unanimous.tsv" "$OUTDIR/all_votes.tsv"
  } > "$OUTDIR/assignment_full.tsv"

  # Summary report
  {
    echo "=========================================="
    echo "  SINE Incremental Assignment Summary"
    echo "=========================================="
    echo "Date: $(date)"
    echo
    echo "Input Files:"
    echo "  Consensus: $CONSENSUS_FILE ($consensus_count subfamilies)"
    echo "  SINEs: $SINE_FILE ($sine_count sequences)"
    echo "  Old assignment: $OLD_ASSIGNMENT"
    echo
    echo "Incremental Mode:"
    echo "  Tier 1 (carried forward): $tier1_n"
    echo "  Tier 2 (re-evaluated):    $tier2_n"
    echo "  Overlapping old assigned:  $overlap_n"
    echo
    echo "Assignment Criteria:"
    echo "  1) 10/10 unanimous votes (same consensus wins all cycles)"
    echo "  2) Bitscore >= ${MIN_REL_BITSCORE} . TopN (N=min(10,count))"
    echo
    echo "Results:"
    echo "  Total sequences: $sine_count"
    echo "  Unanimous 10/10: $unanimous_count"
    echo "  Assigned (passed threshold): $assigned_count ($(awk -v a=$assigned_count -v t=$sine_count 'BEGIN{printf "%.1f%%", 100*a/t}'))"
    echo "  Unassigned: $unassigned_count ($(awk -v u=$unassigned_count -v t=$sine_count 'BEGIN{printf "%.1f%%", 100*u/t}'))"
    echo
    echo "Assignments per Subfamily:"
    awk -F'\t' 'NR>1{printf "  %-20s %6d\n", $1, $2}' "$OUTDIR/assignment_stats.tsv"
    echo
    echo "Output Files:"
    echo "  assigned.fasta         - All assigned sequences with headers: >seqID|subfamily|bitscore"
    echo "  unassigned.fasta       - Sequences that failed assignment"
    echo "  assignment_stats.tsv   - Statistics per subfamily"
    echo "  assignment_full.tsv    - Detailed assignment info for all sequences (incl. non-unanimous)"
    echo "  all_votes.tsv          - Per-sequence best vote (all sequences)"
    echo "  subfamilies/           - Directory with per-subfamily FASTA files"
    echo
  } > "$OUTDIR/summary.txt"

  cat "$OUTDIR/summary.txt"
  log "Done! Output written to: $OUTDIR"
}

main(){
  check_deps
  parse_args "$@"
  if [[ "$INCREMENTAL" == true ]]; then
    run_incremental_assignment
  else
    run_assignment
  fi
}

main "$@"
