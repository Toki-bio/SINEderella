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

STEP2_WORK=""   # never reuse an inherited TMPDIR: cleanup would delete it
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
cleanup(){ [[ -n "$STEP2_WORK" && -d "$STEP2_WORK" ]] && rm -rf "$STEP2_WORK"; }
trap cleanup EXIT
STEP2_TOOLS="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)/tools"

# Unit vs composite (owner 2026-11-01; rsi r6 vs r5_r6_P26, Tal rsi_composites_test/r6_P26_20261101/): a composite holds
# its unit entire, so a copy of the unit alone aligns to both with the same alignment and the cycles split by -z 11 noise
# (rsi: 1,188 standalone r6 copies split their votes with P26 and stayed unassigned; ~60 more went to P26 10/10). Rule,
# in every cycle: when the best and the runner-up consensus are a unit and a composite containing it, the composite wins
# only if the copy's alignment to it covers >= UNIT_COVER_OTHER (0.5) of the composite's sequence outside that unit;
# otherwise the unit wins. Units and their intervals come from tools/composite_units.py --bank on the bank (no declared
# table); a composite with < 30 bp outside the unit is left to the bitscore. UNIT_RULE=0 turns it off.
UNIT_RULE="${UNIT_RULE:-1}"; UNIT_COVER_OTHER="${UNIT_COVER_OTHER:-0.5}"
unit_table(){   # BANK.fa OUT.tsv: composite_units.py --bank, or an empty table with UNIT_RULE=0
  case "$UNIT_RULE" in 0) : > "$2"; return ;; 1) ;; *) die "UNIT_RULE must be 0 or 1 (got: $UNIT_RULE)" ;; esac
  python3 "$STEP2_TOOLS/composite_units.py" --bank "$1" "$2" "$THREADS" >&2 \
    || die "tools/composite_units.py --bank failed (UNIT_RULE=0 votes without the unit/composite rule)"
}
cycle_best(){   # UNITS.tsv RAW.m8 >> votes: best and runner-up consensus per copy in one cycle, unit/composite rule applied
  awk -v F="$UNIT_COVER_OTHER" -v SW="$TMPDIR/unit_rule_swaps.tsv" 'BEGIN{FS="\t"}
    FILENAME == ARGV[1] {
      if (FNR == 1 || $2 != "composite") next
      if ($4 + 0 > clen[$1]) clen[$1] = $4 + 0
      if ($5 != "-") { k = $5 SUBSEP $1; ulen[k] += $4 - $3 + 1; iv[k] = iv[k] " " $3 "-" $4 }
      next }
    { seq=$2; cons=$1; bs=$12+0
      key=seq SUBSEP cons
      # pair_sp: where the best HSP of the pair lies on the consensus (m8 columns 7-8, the query is the consensus); read by
      # tools/rollup_families.py to find copies that cover only a part of their consensuses (partial copies)
      if(!(key in pair_bs) || bs > pair_bs[key]){ pair_bs[key]=bs; pair_sp[key]=$7 "\t" $8 } }
    # a copy is the composite only if its alignment to it covers >= F of the composite outside the unit AND >= F of the
    # unit stretches themselves: a real composite copy spans both parts, a standalone copy only one. The second condition
    # matters when the unit resembles the other part of the composite (rsi r5 vs r5_r5_P48 = r5 + r6, r5 97 percent like
    # r6): a standalone r5 can align onto the r6 part of P48, covering the "other" part without touching the unit (first
    # version, rsi 2026-11-02: 1,522 of 1,927 r5 copies went to P48)
    function other_cov(u, c, sq,   k, s, e, n, i, iv2, ab, cov, covu, other) {   # -1: rule does not apply; else 1 = composite wins
      k = u SUBSEP c; if (!(k in ulen)) return -1
      other = clen[c] - ulen[k]; if (other < 30) return -1
      split(pair_sp[sq SUBSEP c], ab, "\t"); s = ab[1] + 0; e = ab[2] + 0; if (s > e) { t = s; s = e; e = t }
      cov = e - s + 1; covu = 0; n = split(iv[k], iv2, " ")
      for (i = 1; i <= n; i++) { split(iv2[i], ab, "-"); lo = (ab[1] > s ? ab[1] : s); hi = (ab[2] < e ? ab[2] : e); if (hi >= lo) { cov -= hi - lo + 1; covu += hi - lo + 1 } }
      return (cov >= F * other && covu >= F * ulen[k]) ? 1 : 0 }
    END{
      for(key in pair_bs){
        split(key, a, SUBSEP); seq=a[1]; cons=a[2]; bs=pair_bs[key]
        if(!(seq in best_bs) || bs > best_bs[seq]){
          if(seq in best_bs){ runner_bs[seq]=best_bs[seq]; runner_cons[seq]=best_cons[seq]; runner_sp[seq]=best_sp[seq] }
          best_bs[seq]=bs; best_cons[seq]=cons; best_sp[seq]=pair_sp[key]
        } else if(!(seq in runner_bs) || bs > runner_bs[seq]){
          runner_bs[seq]=bs; runner_cons[seq]=cons; runner_sp[seq]=pair_sp[key]
        }
      }
      for(seq in best_cons){
        if (seq in runner_cons) {
          b = best_cons[seq]; r = runner_cons[seq]; win = ""
          x = other_cov(b, r, seq); if (x == 1) win = r      # best = unit, runner = composite it is in, copy reaches the rest
          y = other_cov(r, b, seq); if (y == 0) win = r      # best = composite, runner = its unit, copy does not reach the rest
          if (win == r) {
            print seq "\t" b "\t" r >> SW
            t = best_bs[seq]; best_bs[seq] = runner_bs[seq]; runner_bs[seq] = t
            t = best_sp[seq]; best_sp[seq] = runner_sp[seq]; runner_sp[seq] = t
            best_cons[seq] = r; runner_cons[seq] = b
          }
        }
        rb = (seq in runner_bs ? runner_bs[seq] : 0)
        rc = (seq in runner_cons ? runner_cons[seq] : "NA")
        printf "%s\t%s\t%d\t%s\t%d\t%s\t%s\n", seq, best_cons[seq], int(best_bs[seq]+0.5), rc, int(rb+0.5), best_sp[seq], (seq in runner_sp ? runner_sp[seq] : "0\t0")
      }
    }' "$1" "$2"
}

check_deps(){
  local req=(seqkit ssearch36 awk sed grep sort python3)   # python3: tools/composite_units.py (unit/composite rule)
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

  TMPDIR="$(mktemp -d "${STEP2_TMP_BASE:-$PWD}/.step2_assign_XXXX")"; STEP2_WORK="$TMPDIR"
  export LC_ALL=C

  # Clean headers (keep only first word after >)
  log "Preparing input files..."
  awk '/^>/{split($1,a," ");print a[1];next}{print}' "$CONSENSUS_FILE" > "$TMPDIR/consensus_clean.fa"
  seqkit seq -g "$SINE_FILE" > "$TMPDIR/sines_clean.fa" 2>"$OUTDIR/seqkit.log"

  local dup; dup="$(awk '/^>/{n[substr($1,2)]++} END{for(k in n) if(n[k]>1) print k}' "$TMPDIR/consensus_clean.fa")"
  [[ -z "$dup" ]] || die "duplicate consensus names: $dup"
  # votes are kept per copy name: two copies with one name would share one verdict (GLM audit 2026-10-19)
  dup="$(awk '/^>/{if(n[substr($1,2)]++){print substr($1,2); exit}}' "$TMPDIR/sines_clean.fa")"
  [[ -z "$dup" ]] || die "duplicate copy names in $SINE_FILE: $dup"
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

  unit_table "$TMPDIR/consensus_clean.fa" "$TMPDIR/units.tsv"; : > "$TMPDIR/unit_rule_swaps.tsv"
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
        >> "$TMPDIR/raw_${cycle}.m8" 2>>"$OUTDIR/ssearch.log" \
        || die "ssearch36 failed (cycle $cycle, ${part##*/}); see $OUTDIR/ssearch.log"
    done
    # ssearch36 exits 0 when it cannot open a file; "***" in its log is then the only sign (healthy runs leave the log empty)
    ! grep -q '\*\*\*' "$OUTDIR/ssearch.log" || die "ssearch36 error in cycle $cycle: $(grep '\*\*\*' "$OUTDIR/ssearch.log" | head -1)"

    # best hit per copy for this cycle and the best-scoring distinct runner-up (LEAK / ambiguity downstream), with the
    # unit/composite rule (cycle_best, top of this file)
    cycle_best "$TMPDIR/units.tsv" "$TMPDIR/raw_${cycle}.m8" >> "$TMPDIR/votes.tsv"
  done
  echo >&2 "  Cycle 10/10 done."
  # unit/composite rule: one row per copy and cycle where it changed the winner (copy, bitscore winner, rule winner)
  cp "$TMPDIR/unit_rule_swaps.tsv" "$OUTDIR/unit_rule_swaps.tsv"; cp "$TMPDIR/units.tsv" "$OUTDIR/unit_rule_units.tsv"
  log "Unit/composite rule: $(cut -f1 "$OUTDIR/unit_rule_swaps.tsv" | sort -u | wc -l) copies, $(wc -l < "$OUTDIR/unit_rule_swaps.tsv") cycle winners changed (unit_rule_swaps.tsv)"

  # Aggregate votes
  log "Aggregating votes and calculating thresholds..."
  awk 'BEGIN{FS=OFS="\t"}{
    key = $1 FS $2
    vote_count[key]++
    bitscore_sum[key] += $3
    runner_bs_sum[key] += $5
    if(!(key in best_runner_bs) || $5+0 > best_runner_bs[key]){
      best_runner_bs[key] = $5+0
      runner_cons[key] = $4
      runner_sp[key] = $8 FS $9
    }
    if(!(key in top_bs) || $3+0 > top_bs[key]){ top_bs[key] = $3+0; sp[key] = $6 FS $7 }
  }
  END{
    # columns 7-10 (2026-10-26): span of the copy on this consensus in its best cycle, and on the runner-up consensus in the cycle
    # with the best runner-up score; columns 1-6 as before
    for(key in vote_count){
      split(key, a, FS)
      rc = (key in runner_cons ? runner_cons[key] : "NA")
      print a[1], a[2], vote_count[key], bitscore_sum[key], runner_bs_sum[key], rc, sp[key], runner_sp[key]
    }
  }' "$TMPDIR/votes.tsv" | sort -k1,1 -k3,3nr -k2,2 > "$TMPDIR/vote_summary.tsv"
  # keep the per-consensus tally (used by tools/rollup_families.py); small next to the m8 files
  gzip -c "$TMPDIR/vote_summary.tsv" > "$OUTDIR/vote_summary.tsv.gz" || true

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
    rbs_sum = $5 + 0
    rcons = $6

    if(votes == 10){
      if(!(seq in best_bs) || bs_sum > best_bs[seq]){
        best_bs[seq] = bs_sum
        best_cons[seq] = cons
        best_rbs[seq] = rbs_sum
        best_rcons[seq] = rcons
      }
    }
  }
  END{
    for(seq in best_cons){
      rc = (seq in best_rcons ? best_rcons[seq] : "NA")
      rb = (seq in best_rbs ? best_rbs[seq] : 0)
      printf "%s\t%s\t%d\t10\t%d\t%s\n", seq, best_cons[seq], int(best_bs[seq]+0.5), int(rb+0.5), rc
    }
  }' "$TMPDIR/vote_summary.tsv" | sort -k2,2 -k3,3nr > "$TMPDIR/unanimous.tsv"

  local unanimous_count
  unanimous_count="$(awk 'END{print NR}' "$TMPDIR/unanimous.tsv")"
  log "Sequences with 10/10 unanimous votes: $unanimous_count"

  # Calculate thresholds per consensus
  log "Calculating bitscore thresholds per subfamily..."
  # N-th best 10-cycle sum per consensus, N = min(10, count); one pass, no per-consensus files
  # (a name such as rnd-1_family-5#SINE/tRNA is not a valid file name)
  sort -t$'\t' -k2,2 -k3,3nr "$TMPDIR/unanimous.tsv" \
    | awk -F'\t' -v m="$MIN_REL_BITSCORE" '{ c=$2; n[c]++; if(n[c]<=10) nth[c]=$3 }
        END{ for(c in n) printf "%s\t%d\t%d\n", c, nth[c], int(nth[c]*m+0.5) }' \
    | sort -k1,1 > "$TMPDIR/thresholds.tsv"

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
    rbs = $5 + 0
    rcons = $6

    if(bs >= threshold[cons]){
      status = "assigned"
    } else {
      status = "rejected_low_bitscore"
    }

    rratio = (bs > 0 ? rbs/bs : 0)
    tag = sprintf("runner=%s;runner_bs=%d;runner_ratio=%.4f", (rcons==""?"NA":rcons), rbs, rratio)

    printf "%s\t%s\t%d\t%d\t%s\t%d\t%s\n", seq, cons, bs, votes, status, threshold[cons], tag
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
      # only copies that passed: until 2026-10-19 rejected_low_bitscore rows went in too (rsi: 1,731 of 65,500)
      if($5 == "assigned"){ cons[$1] = $2; bs[$1] = $3 }
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

      # Extract sequences with renamed headers
      awk 'BEGIN{FS="\t"}
      NR==FNR{
        if($5 == "assigned"){ cons[$1] = $2; bs[$1] = $3 }   # only copies that passed (rejected_low_bitscore excluded)
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
        > "$OUTDIR/subfamilies/${subfamily//\//_}.fasta"

      local count
      count="$(count_fasta "$OUTDIR/subfamilies/${subfamily//\//_}.fasta")"
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
    echo -e "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\tRunnerInfo"
    # 1) unanimous sequences (assigned or rejected_low_bitscore)
    cat "$TMPDIR/assignment_details.tsv"
    # 2) non-unanimous sequences: best vote from all_votes.tsv minus unanimous
    # (no per-cycle runner-up data survives past all_votes.tsv, so tag is NA here)
    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR{ seen[$1]=1; next }
    !($1 in seen){
      printf "%s\t%s\t%d\t%d\t%s\t%s\t%s\n", $1, $2, $4, $3, "no_unanimous", "NA", "NA"
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

  TMPDIR="$(mktemp -d "${STEP2_TMP_BASE:-$PWD}/.step2_incr_XXXX")"; STEP2_WORK="$TMPDIR"
  export LC_ALL=C

  log "Preparing input files..."
  awk '/^>/{split($1,a," ");print a[1];next}{print}' "$CONSENSUS_FILE" > "$TMPDIR/consensus_clean.fa"
  seqkit seq -g "$SINE_FILE" > "$TMPDIR/sines_clean.fa" 2>"$OUTDIR/seqkit.log"

  local dup; dup="$(awk '/^>/{n[substr($1,2)]++} END{for(k in n) if(n[k]>1) print k}' "$TMPDIR/consensus_clean.fa")"
  [[ -z "$dup" ]] || die "duplicate consensus names: $dup"
  # votes are kept per copy name: two copies with one name would share one verdict (GLM audit 2026-10-19)
  dup="$(awk '/^>/{if(n[substr($1,2)]++){print substr($1,2); exit}}' "$TMPDIR/sines_clean.fa")"
  [[ -z "$dup" ]] || die "duplicate copy names in $SINE_FILE: $dup"
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
  # the contig is everything before the LAST ':' ahead of start-end(strand): contig names may hold ':' (header cleaning keeps
  # it), and cutting at the first ':' gave contig "ctg1", start 0 for "ctg1:apple:100-200(+)" (GLM audit 2026-11-08, finding 72)
  awk '{
    id=$1
    if(!match(id,/^.*:[0-9]+-[0-9]+\([+-]\)/)){ print "step2: cannot read locus id " id > "/dev/stderr"; exit 1 }
    pre=substr(id,1,RLENGTH); k=0
    for(i=length(pre);i>0;i--) if(substr(pre,i,1)==":"){ k=i; break }
    chr=substr(pre,1,k-1); rng=substr(pre,k+1); sub(/\([+-]\)$/,"",rng)
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
      rb=0; rc="NA"
      if(match($7,/runner_bs=[0-9]+/)) rb=substr($7,RSTART+10,RLENGTH-10)
      if(match($7,/runner=[^;]+/))     rc=substr($7,RSTART+7,RLENGTH-7)
      printf "%s\t%s\t%d\t10\t%d\t%s\n", $1, $2, $3, rb, rc
    }' "$TMPDIR/tier1_final.ids" "$OLD_ASSIGNMENT" > "$TMPDIR/tier1_unanimous.tsv"

  #---------------------------------------------------------------------------
  # Run 10-cycle ssearch36 on Tier 2 only
  #---------------------------------------------------------------------------
  : > "$TMPDIR/tier2_unanimous.tsv"
  : > "$OUTDIR/unit_rule_swaps.tsv"   # no tier 2 copies: no swaps, and no stale file from an earlier run
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

    unit_table "$TMPDIR/consensus_clean.fa" "$TMPDIR/units.tsv"; : > "$TMPDIR/unit_rule_swaps.tsv"
    : > "$TMPDIR/votes.tsv"
    for cycle in $(seq 1 10); do
      printf "  Cycle %2d/10\r" "$cycle" >&2
      : > "$TMPDIR/raw_${cycle}.m8"
      for part in "$TMPDIR"/splits/*.fasta; do
        ssearch36 -g -3 -T "$THREADS" -Q -n -z 11 -E 2 -w 95 -W 70 -m 8 \
          "$TMPDIR/consensus_clean.fa" "$part" \
          >> "$TMPDIR/raw_${cycle}.m8" 2>>"$OUTDIR/ssearch.log" \
          || die "ssearch36 failed (cycle $cycle, ${part##*/}); see $OUTDIR/ssearch.log"
      done
      # ssearch36 exits 0 when it cannot open a file; "***" in its log is then the only sign (healthy runs leave the log empty)
      ! grep -q '\*\*\*' "$OUTDIR/ssearch.log" || die "ssearch36 error in cycle $cycle: $(grep '\*\*\*' "$OUTDIR/ssearch.log" | head -1)"
      cycle_best "$TMPDIR/units.tsv" "$TMPDIR/raw_${cycle}.m8" >> "$TMPDIR/votes.tsv"
    done
    echo >&2 "  Cycle 10/10 done."
    cp "$TMPDIR/unit_rule_swaps.tsv" "$OUTDIR/unit_rule_swaps.tsv"
    log "Unit/composite rule (tier 2): $(cut -f1 "$TMPDIR/unit_rule_swaps.tsv" | sort -u | wc -l) copies changed"

    # Aggregate Tier 2 votes
    awk 'BEGIN{FS=OFS="\t"}{
      key = $1 FS $2; vote_count[key]++; bitscore_sum[key] += $3; runner_bs_sum[key] += $5
      if(!(key in best_runner_bs) || $5+0 > best_runner_bs[key]){ best_runner_bs[key] = $5+0; runner_cons[key] = $4; runner_sp[key] = $8 FS $9 }
      if(!(key in top_bs) || $3+0 > top_bs[key]){ top_bs[key] = $3+0; sp[key] = $6 FS $7 }
    }
    END{
      for(key in vote_count){
        split(key, a, FS)
        print a[1], a[2], vote_count[key], bitscore_sum[key], runner_bs_sum[key], runner_cons[key], sp[key], runner_sp[key]   # 7-10: spans, as in tier 1
      }
    }' "$TMPDIR/votes.tsv" | sort -k1,1 -k3,3nr -k2,2 > "$TMPDIR/tier2_vote_summary.tsv"
    gzip -c "$TMPDIR/tier2_vote_summary.tsv" > "$OUTDIR/vote_summary.tsv.gz" || true

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
          best_bs[seq] = bs_sum; best_cons[seq] = cons; rbs[seq] = $5; rc[seq] = $6
        }
      }
    }
    END{
      for(seq in best_cons)
        printf "%s\t%s\t%d\t10\t%d\t%s\n", seq, best_cons[seq], int(best_bs[seq]+0.5), int(rbs[seq]+0.5), rc[seq]
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
  # N-th best 10-cycle sum per consensus, N = min(10, count); one pass, no per-consensus files
  # (a name such as rnd-1_family-5#SINE/tRNA is not a valid file name)
  sort -t$'\t' -k2,2 -k3,3nr "$TMPDIR/unanimous.tsv" \
    | awk -F'\t' -v m="$MIN_REL_BITSCORE" '{ c=$2; n[c]++; if(n[c]<=10) nth[c]=$3 }
        END{ for(c in n) printf "%s\t%d\t%d\n", c, nth[c], int(nth[c]*m+0.5) }' \
    | sort -k1,1 > "$TMPDIR/thresholds.tsv"

  # Apply thresholds
  log "Applying thresholds to filter assignments..."
  awk 'BEGIN{FS=OFS="\t"}
  NR==FNR{
    threshold[$1] = $3; topN[$1] = $2; next
  }
  {
    seq = $1; cons = $2; bs = $3 + 0; votes = $4; rbs = $5 + 0; rcons = $6
    if(bs >= threshold[cons]) status = "assigned"
    else status = "rejected_low_bitscore"
    tag = sprintf("runner=%s;runner_bs=%d;runner_ratio=%.4f", (rcons==""?"NA":rcons), rbs, (bs > 0 ? rbs/bs : 0))
    printf "%s\t%s\t%d\t%d\t%s\t%d\t%s\n", seq, cons, bs, votes, status, threshold[cons], tag
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
      if($5 == "assigned"){ cons[$1] = $2; bs[$1] = $3 }; next     # only copies that passed (rejected_low_bitscore excluded)
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
        if($5 == "assigned"){ cons[$1] = $2; bs[$1] = $3 }; next   # only copies that passed (rejected_low_bitscore excluded)
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
        > "$OUTDIR/subfamilies/${subfamily//\//_}.fasta"
      local count
      count="$(count_fasta "$OUTDIR/subfamilies/${subfamily//\//_}.fasta")"
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
    echo -e "Sequence\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\tRunnerInfo"
    cat "$TMPDIR/assignment_details.tsv"
    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR{ seen[$1]=1; next }
    !($1 in seen){
      printf "%s\t%s\t%d\t%d\t%s\t%s\t%s\n", $1, $2, $4, $3, "no_unanimous", "NA", "NA"
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
