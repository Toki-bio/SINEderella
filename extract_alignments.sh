#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# extract_alignments.sh
#
# Generate representative alignments from a SINEderella multi-run project.
# See PLAN_alignment_viewer.md for design rationale.
#
# Usage:
#   extract_alignments.sh <consensus_bank.fa> <multi_run_dir> [output_dir]
#
# Output structure (per species x subfamily):
#   <output_dir>/<Species>/<subfamily>.core.fa        -- Tier A: random 50, flanked 30L+70R
#   <output_dir>/<Species>/<subfamily>.best50.fa      -- Tier B: top 50 by bitscore, flanked
#   <output_dir>/<Species>/<subfamily>.subfam.fa      -- Tier C: SubFam consensuses (>400)
#   <output_dir>/<Species>/<subfamily>.evidence.txt   -- Tier D: zero-hit evidence
#
# Flanks (30bp left + 70bp right) are included in Tier A and B alignments
# for TSD inspection. Copies come from step2 subfamilies/*.fasta; flanked
# versions are re-extracted from the genome using coordinates from headers.
#
# MAFFT parameters (PLAN section 4):
#   --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --preservecase
#
# Dependencies: mafft, bedtools, samtools, ssearch36, shuf
# Optional: SubFam (for Tier C)
#
# TODO (future):
#   - If consensus in database is annotated with functional parts (tRNA-derived
#     region, body, tail, promoters A/B boxes, TSD in flanks), reflect these
#     features in the alignment output as annotation tracks. Could be done via
#     a BED-like annotation file per subfamily consensus that maps functional
#     regions to column ranges in the MSA. The viewer would then render colored
#     bars above the alignment showing tRNA / body / tail / promoter / TSD.
#   - Build distance matrix / tree / clustering based on SINE subfamily copy
#     count profiles across taxa (presence/absence + copy number vector per
#     taxon). Enable user-selectable subset of taxa and/or SINE subfamilies
#     for tree construction. Could use Jaccard, Bray-Curtis, or Euclidean
#     distance on log-transformed copy counts, then NJ or UPGMA clustering.
###############################################################################

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ echo "ERROR: $*" >&2; exit 1; }

# -- Arguments ----------------------------------------------------------------
CONS_BANK="${1:-}"
MULTI_DIR="${2:-}"
OUT_BASE="${3:-}"

[[ -n "$CONS_BANK" ]] || die "Usage: $0 <consensus_bank.fa> <multi_run_dir> [output_dir]"
[[ -n "$MULTI_DIR" ]] || die "Usage: $0 <consensus_bank.fa> <multi_run_dir> [output_dir]"

[[ -f "$CONS_BANK" ]] || die "Consensus bank not found: $CONS_BANK"
[[ -d "$MULTI_DIR" ]] || die "Multi-run directory not found: $MULTI_DIR"

CONS_BANK="$(readlink -f "$CONS_BANK")"
MULTI_DIR="$(readlink -f "$MULTI_DIR")"
[[ -n "$OUT_BASE" ]] && OUT_BASE="$(readlink -f "$OUT_BASE")" || OUT_BASE="$MULTI_DIR/alignments"

THREADS=$(nproc 2>/dev/null || echo 1)

# -- Find SubFam --------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBFAM_BIN=""
if [[ -f "$SCRIPT_DIR/SubFam" ]]; then
  SUBFAM_BIN="$SCRIPT_DIR/SubFam"
  chmod +x "$SUBFAM_BIN" 2>/dev/null || true
elif command -v SubFam >/dev/null 2>&1; then
  SUBFAM_BIN="$(command -v SubFam)"
fi

# -- Find sear (for Tier D evidence) ------------------------------------------
SEAR_BIN=""
if [[ -f "$SCRIPT_DIR/sear" ]]; then
  SEAR_BIN="$SCRIPT_DIR/sear"
  chmod +x "$SEAR_BIN" 2>/dev/null || true
elif command -v sear >/dev/null 2>&1; then
  SEAR_BIN="$(command -v sear)"
fi

# -- Verify dependencies ------------------------------------------------------
for cmd in mafft bedtools samtools ssearch36 shuf; do
  command -v "$cmd" >/dev/null 2>&1 || die "$cmd not in PATH"
done

# -- MAFFT wrapper (PLAN section 4 params) ------------------------------------
run_mafft(){
  local input="$1" output="$2"
  mafft --thread "$THREADS" --threadtb "$THREADS" --threadit "$THREADS" \
        --localpair --maxiterate 1000 --ep 0.123 \
        --nuc --reorder --preservecase --quiet \
        "$input" > "$output" 2>/dev/null
}

# -- Flank post-processor ------------------------------------------------------
# After MAFFT alignment: the consensus (1st input seq) defines the SINE body.
# Columns before the first non-gap consensus base = left flank.
# Columns after  the last  non-gap consensus base = right flank.
# In each flank region per copy:
#   - remove internal gaps (collapse bases together)
#   - lowercase the bases
#   - pad with gaps on the outer side so column count stays the same
# The consensus row is never the first in MAFFT --reorder output, so we
# match it by name (first field of the header line in the original input).
postprocess_flanks(){
  local aln_file="$1" cons_name="$2"
  [[ -s "$aln_file" ]] || return 0

  local tmp="${aln_file}.ppf"
  awk -v cons_name="$cons_name" '
    # ── pass 1: slurp all sequences ──
    /^>/ {
      if (NR>1) seqs[n] = seq
      n++; hdr[n] = $0; seq = ""
      # check if this is the consensus (match first word after >)
      h = $0; sub(/^>/, "", h); sub(/[[:space:]].*/, "", h)
      if (h == cons_name) cons_idx = n
      next
    }
    { seq = seq $0 }
    END {
      seqs[n] = seq
      if (cons_idx == 0) {
        # consensus not found – output unchanged
        for (i=1; i<=n; i++) {
          print hdr[i]
          s = seqs[i]; L = length(s)
          for (j=1; j<=L; j+=80) print substr(s, j, 80)
        }
        exit
      }

      cs = seqs[cons_idx]
      alen = length(cs)

      # find left/right body boundaries (1-based column indices)
      lbound = 0; rbound = 0
      for (j=1; j<=alen; j++) {
        c = substr(cs, j, 1)
        if (c != "-" && c != ".") { if (!lbound) lbound = j; rbound = j }
      }
      if (lbound == 0) lbound = 1
      if (rbound == 0) rbound = alen

      for (i=1; i<=n; i++) {
        print hdr[i]
        s = seqs[i]

        if (i == cons_idx || alen == 0) {
          for (j=1; j<=alen; j+=80) print substr(s, j, 80)
          continue
        }

        # ── left flank: columns 1..lbound-1 ──
        lf_bases = ""
        lf_len = lbound - 1
        for (j=1; j<lbound; j++) {
          c = substr(s, j, 1)
          if (c != "-" && c != ".") lf_bases = lf_bases tolower(c)
        }
        # pad with leading gaps so ungapped bases are right-justified
        lf_pad = lf_len - length(lf_bases)
        lf_out = ""
        for (j=1; j<=lf_pad; j++) lf_out = lf_out "-"
        lf_out = lf_out lf_bases

        # ── body: columns lbound..rbound (unchanged) ──
        body = substr(s, lbound, rbound - lbound + 1)

        # ── right flank: columns rbound+1..alen ──
        rf_bases = ""
        rf_len = alen - rbound
        for (j=rbound+1; j<=alen; j++) {
          c = substr(s, j, 1)
          if (c != "-" && c != ".") rf_bases = rf_bases tolower(c)
        }
        # pad with trailing gaps so ungapped bases are left-justified
        rf_pad = rf_len - length(rf_bases)
        rf_out = rf_bases
        for (j=1; j<=rf_pad; j++) rf_out = rf_out "-"

        full = lf_out body rf_out
        for (j=1; j<=length(full); j+=80) print substr(full, j, 80)
      }
    }
  ' "$aln_file" > "$tmp" && mv "$tmp" "$aln_file"
}

# -- Helpers -------------------------------------------------------------------
count_fasta(){ awk 'BEGIN{c=0}/^>/{c++}END{printf "%d\n",c}' "$1"; }

# Sample N random FASTA sequences (linearise, shuf, delinearise)
sample_fasta(){
  local n="$1" infile="$2" outfile="$3"
  awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next }
       { printf("%s",$0) }
       END  { printf("\n") }' "$infile" > "$outfile.lin"
  shuf "$outfile.lin" | head -n "$n" > "$outfile.sel" || true
  awk 'BEGIN{FS="\t"}{printf("%s\n%s\n",$1,$2)}' "$outfile.sel" > "$outfile"
  rm -f "$outfile.lin" "$outfile.sel"
}

# Top N sequences by bitscore (parsed from FASTA header: >loc|subfam|BITS)
# Note: sort writes to intermediate file to avoid SIGPIPE with pipefail
top_n_by_bitscore(){
  local n="$1" infile="$2" outfile="$3"
  awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next }
       { printf("%s",$0) }
       END  { printf("\n") }' "$infile" > "$outfile.lin"
  # Sort by 3rd pipe-field (bitscore) descending — write to file to avoid SIGPIPE
  awk 'BEGIN{FS="\t"}{
    hdr=$1; split(hdr, p, "|");
    bs = p[3]+0;
    printf "%d\t%s\n", bs, $0
  }' "$outfile.lin" | sort -k1,1nr > "$outfile.sorted"
  head -n "$n" "$outfile.sorted" | cut -f2- \
  | awk 'BEGIN{FS="\t"}{printf("%s\n%s\n",$1,$2)}' > "$outfile"
  rm -f "$outfile.lin" "$outfile.sorted"
}

# Extract single consensus for a subfamily from the bank
extract_consensus(){
  local subfam="$1" bank="$2" outfile="$3"
  awk -v sf="$subfam" '
    /^>/{
      name = substr($1,2)
      active = (name == sf) ? 1 : 0
    }
    active{print}
  ' "$bank" > "$outfile"
}

# Parse FASTA headers to BED (0-based coords from bedtools getfasta headers)
# Header format: >ctg:start-end(strand)|Subfam|Bits
fasta_headers_to_bed(){
  local infile="$1" outfile="$2"
  grep "^>" "$infile" | sed 's/^>//' | awk '{
    id = $1
    n = split(id, parts, "|")
    loc = parts[1]

    strand = "+"
    if (sub(/\([+-]\)$/, "", loc)) {
      i = index(parts[1], "(")
      strand = substr(parts[1], i+1, 1)
    }

    # Find last colon (handles chr names containing colons)
    last_colon = 0; rest = loc
    while ((p = index(rest, ":")) > 0) {
      last_colon += p; rest = substr(rest, p+1)
    }
    if (last_colon == 0) next

    ctg = substr(loc, 1, last_colon-1)
    coords = substr(loc, last_colon+1)
    split(coords, c, "-")
    if (c[1]+0 >= 0 && c[2]+0 > 0) {
      printf "%s\t%s\t%s\t%s\t0\t%s\n", ctg, c[1], c[2], id, strand
    }
  }' > "$outfile"
}

# Re-extract sequences with 30L+70R flanks from genome.
# Input: FASTA file (with coordinate headers), genome path, output path, tmpdir
# Falls back to original sequences if genome extraction fails.
extract_flanked(){
  local infile="$1" genome="$2" outfile="$3" tmpd="$4"

  fasta_headers_to_bed "$infile" "$tmpd/flank.bed"
  if [[ ! -s "$tmpd/flank.bed" ]]; then
    cp "$infile" "$outfile"
    return
  fi

  awk -v OFS='\t' '{print $1,$2}' "${genome}.fai" > "$tmpd/genome.sizes"

  bedtools slop -s -l 30 -r 70 -g "$tmpd/genome.sizes" \
    -i "$tmpd/flank.bed" > "$tmpd/flanked.bed" 2>/dev/null || true

  if [[ -s "$tmpd/flanked.bed" ]]; then
    bedtools getfasta -s -nameOnly -fi "$genome" \
      -bed "$tmpd/flanked.bed" > "$outfile" 2>/dev/null || true
  fi

  # Fall back to original if extraction failed
  if [[ ! -s "$outfile" ]]; then
    cp "$infile" "$outfile"
  fi
}

# -- Main loop -----------------------------------------------------------------
log "============================================================"
log "extract_alignments.sh"
log "  Consensus bank: $CONS_BANK"
log "  Multi-run dir:  $MULTI_DIR"
log "  Output dir:     $OUT_BASE"
log "  Threads:        $THREADS"
log "  SubFam:         ${SUBFAM_BIN:-NOT FOUND}"
log "  sear:           ${SEAR_BIN:-NOT FOUND (Tier D will be skipped)}"
log "============================================================"

mkdir -p "$OUT_BASE"

n_species=0
n_total=0
n_tierA=0
n_tierB=0
n_tierC=0
n_tierD=0

for sp_dir in "$MULTI_DIR"/*/; do
  [[ -d "$sp_dir" ]] || continue
  sp="$(basename "$sp_dir")"
  [[ "$sp" == "summary" || "$sp" == "alignments" ]] && continue

  # Find latest run_* directory
  latest="$(ls -dt "$sp_dir"/run_* 2>/dev/null | head -n1 || true)"
  [[ -n "$latest" && -d "$latest" ]] || continue

  # Find step2 output directory
  step2_out="$(ls -dt "$latest"/step2/step2_output* 2>/dev/null | head -n1 || true)"

  # Find subfamilies directory (from step2)
  subfam_dir=""
  [[ -n "$step2_out" && -d "$step2_out/subfamilies" ]] && subfam_dir="$step2_out/subfamilies"

  # Find unassigned files for soft-assigned copy inclusion
  unassigned_tsv=""
  unassigned_fasta=""
  [[ -n "$step2_out" && -f "$step2_out/unassigned.tsv" ]] && unassigned_tsv="$step2_out/unassigned.tsv"
  [[ -n "$step2_out" && -f "$step2_out/unassigned.fasta" ]] && unassigned_fasta="$step2_out/unassigned.fasta"

  # Find genome for flank re-extraction
  genome=""
  for candidate in "$latest/genome.clean.fa" "$latest/genome.fa"; do
    [[ -f "$candidate" ]] && genome="$candidate" && break
  done

  # Index genome if needed
  if [[ -n "$genome" && ! -f "${genome}.fai" ]]; then
    samtools faidx "$genome" 2>/dev/null || true
  fi

  n_species=$((n_species + 1))
  sp_out="$OUT_BASE/$sp"
  mkdir -p "$sp_out"

  # Shared sear working directory per species: genome splits are created once
  # and reused for all Tier D calls (avoids re-splitting 1.5 GB genome per subfamily)
  sear_workdir=""
  if [[ -n "$genome" && -n "$SEAR_BIN" ]]; then
    sear_workdir=$(mktemp -d "${TMPDIR:-/tmp}/ea_sear_${sp}.XXXXXX")
  fi

  log "[$sp] Processing..."

  # Build list of ALL subfamilies from consensus bank
  mapfile -t all_subfams < <(grep "^>" "$CONS_BANK" | sed 's/^>//' | awk '{print $1}')

  # Build set of subfamilies that have actual copies (from step2 output)
  declare -A has_copies=()
  declare -A has_soft_copies=()
  if [[ -n "$subfam_dir" && -d "$subfam_dir" ]]; then
    for sf_fasta in "$subfam_dir"/*.fasta; do
      [[ -f "$sf_fasta" ]] || continue
      sf_name="$(basename "$sf_fasta" .fasta)"
      nc=$(count_fasta "$sf_fasta")
      [[ $nc -ge 1 ]] && has_copies["$sf_name"]="$nc"
    done
  fi

  # Count soft-assigned copies per subfamily from unassigned.tsv
  if [[ -n "$unassigned_tsv" && -f "$unassigned_tsv" ]]; then
    while IFS=$'\t' read -r _sf _ct; do
      [[ -n "$_sf" && "$_ct" -gt 0 ]] 2>/dev/null && has_soft_copies["$_sf"]="$_ct"
    done < <(awk -F'\t' 'NR>1 && $6!="." {count[$6]++} END{for(k in count) print k, count[k]}' "$unassigned_tsv")
  fi

  for sf_name in "${all_subfams[@]}"; do
    nseqs="${has_copies[$sf_name]:-0}"
    nsoft="${has_soft_copies[$sf_name]:-0}"
    ntotal=$((nseqs + nsoft))

    tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/ea_${sp}_${sf_name}.XXXXXX")

    # Extract this subfamily's consensus
    extract_consensus "$sf_name" "$CONS_BANK" "$tmpdir/consensus.fa"

    if [[ $ntotal -eq 0 ]]; then
      # -- Tier D: Zero hits -- produce evidence -----------------------------
      log "  [$sf_name] 0 copies -> Tier D (evidence)"

      if [[ -n "$genome" && -n "$SEAR_BIN" && -n "$sear_workdir" ]]; then
        # Use sear with relaxed thresholds (0.5 id, 30 len, 20 best)
        # sear does chunked ssearch36 with fragment scanning
        # Runs in shared sear_workdir so genome splits (*.2k.part_*.bnk) are
        # created once and reused for all Tier D calls within this species
        (
          cd "$sear_workdir"
          "$SEAR_BIN" "$tmpdir/consensus.fa" "$genome" 0.5 30 20 >/dev/null 2>&1
        ) || true

        # Collect BED results (sear outputs <taxname>-<queryname>.bed)
        bed_hits=""
        bed_hits=$(cat "$sear_workdir"/*.bed 2>/dev/null | head -n 10 || true)

        {
          echo "# Tier D: no SINE copies assigned – sear evidence (relaxed: id>=0.5, len>=30)"
          echo "# chrom	start	end	name	score	strand"
          if [[ -n "$bed_hits" ]]; then
            echo "$bed_hits"
          else
            echo "# No hits found"
          fi
        } > "$sp_out/${sf_name}.evidence.txt"

        n_tierD=$((n_tierD + 1))
      elif [[ -n "$genome" ]]; then
        # Fallback: no sear available
        echo "# Tier D: sear not available, no evidence generated" \
          > "$sp_out/${sf_name}.evidence.txt"
        n_tierD=$((n_tierD + 1))
      else
        echo "No genome available for evidence search" > "$sp_out/${sf_name}.evidence.txt"
      fi

      rm -rf "$tmpdir"
      continue
    fi

    # Build combined FASTA with firm + soft copies for this subfamily
    # Firm copies come from subfamilies/*.fasta (headers: >seqid|subfam|bitscore)
    # Soft copies come from unassigned.fasta, filtered by Soft_Subfamily in unassigned.tsv
    sf_fasta="$subfam_dir/${sf_name}.fasta"

    # Start with firm copies (if any)
    : > "$tmpdir/combined_pool.fa"
    if [[ $nseqs -gt 0 && -f "$sf_fasta" ]]; then
      cp "$sf_fasta" "$tmpdir/combined_pool.fa"
    fi

    # Append soft copies (if any)
    if [[ $nsoft -gt 0 && -n "$unassigned_tsv" && -n "$unassigned_fasta" ]]; then
      # Extract SeqIDs of soft-assigned copies for this subfamily
      awk -F'\t' -v sf="$sf_name" 'NR>1 && $6==sf {print $1}' "$unassigned_tsv" \
        > "$tmpdir/soft_ids.txt"
      if [[ -s "$tmpdir/soft_ids.txt" ]]; then
        # Extract matching sequences from unassigned.fasta and rename headers
        # Soft copies get bitscore=0 marker and |SOFT suffix for identification
        awk 'BEGIN{FS="\t"}
          NR==FNR{ ids[$1]=1; next }
          /^>/{
            name=substr($0,2); split(name,a," "); seqid=a[1]
            if(seqid in ids){
              printf ">%s|%s|0|SOFT\n", seqid, SF
              active=1
            } else { active=0 }
            next
          }
          { if(active) print }
        ' SF="$sf_name" "$tmpdir/soft_ids.txt" "$unassigned_fasta" \
          >> "$tmpdir/combined_pool.fa"
      fi
    fi

    ntotal=$(count_fasta "$tmpdir/combined_pool.fa")
    if [[ $ntotal -eq 0 ]]; then
      rm -rf "$tmpdir"
      continue
    fi

    log "  [$sf_name] $ntotal copies (firm:$nseqs soft:$nsoft)"

    # -----------------------------------------------------------------------
    # Tier A: Random 50 copies with 30L+70R flanks
    # Source: combined pool (firm + soft)
    # Flanks: re-extract from genome with bedtools slop -l 30 -r 70
    # -----------------------------------------------------------------------
    if [[ $ntotal -le 50 ]]; then
      cp "$tmpdir/combined_pool.fa" "$tmpdir/core_sample.fa"
    else
      sample_fasta 50 "$tmpdir/combined_pool.fa" "$tmpdir/core_sample.fa"
    fi

    # Re-extract with flanks if genome is available
    if [[ -n "$genome" && -f "${genome}.fai" ]]; then
      extract_flanked "$tmpdir/core_sample.fa" "$genome" "$tmpdir/tierA_seqs.fa" "$tmpdir"
    else
      cp "$tmpdir/core_sample.fa" "$tmpdir/tierA_seqs.fa"
    fi

    cat "$tmpdir/consensus.fa" "$tmpdir/tierA_seqs.fa" > "$tmpdir/tierA_input.fa"
    run_mafft "$tmpdir/tierA_input.fa" "$sp_out/${sf_name}.core.fa" || {
      log "    WARNING: mafft failed for Tier A $sf_name"
    }
    postprocess_flanks "$sp_out/${sf_name}.core.fa" "$sf_name"
    n_tierA=$((n_tierA + 1))

    # -----------------------------------------------------------------------
    # Tier B: Best 50 by bitscore with 30L+70R flanks
    # Source: combined pool (firm + soft); soft have bitscore=0 so rank last
    # -----------------------------------------------------------------------
    top_n_by_bitscore 50 "$tmpdir/combined_pool.fa" "$tmpdir/best50_core.fa"

    if [[ -s "$tmpdir/best50_core.fa" ]]; then
      if [[ -n "$genome" && -f "${genome}.fai" ]]; then
        extract_flanked "$tmpdir/best50_core.fa" "$genome" "$tmpdir/tierB_seqs.fa" "$tmpdir"
      else
        cp "$tmpdir/best50_core.fa" "$tmpdir/tierB_seqs.fa"
      fi

      cat "$tmpdir/consensus.fa" "$tmpdir/tierB_seqs.fa" > "$tmpdir/tierB_input.fa"
      run_mafft "$tmpdir/tierB_input.fa" "$sp_out/${sf_name}.best50.fa" || {
        log "    WARNING: mafft failed for Tier B $sf_name"
      }
      postprocess_flanks "$sp_out/${sf_name}.best50.fa" "$sf_name"
      n_tierB=$((n_tierB + 1))
    fi

    # -----------------------------------------------------------------------
    # Tier C: SubFam consensuses (only if >400 copies including soft)
    # Subsample 2500 copies, BnkSz=50 -> ~50 consensuses (PLAN Tier C)
    # -----------------------------------------------------------------------
    if [[ $ntotal -gt 400 && -n "$SUBFAM_BIN" ]]; then
      log "    $sf_name: $ntotal copies (firm:$nseqs soft:$nsoft) -> Tier C (SubFam, subsample 2500)"

      if [[ $ntotal -le 2500 ]]; then
        cp "$tmpdir/combined_pool.fa" "$tmpdir/subfam_input.fa"
      else
        sample_fasta 2500 "$tmpdir/combined_pool.fa" "$tmpdir/subfam_input.fa"
      fi

      subfam_workdir="$tmpdir/subfam_work"
      mkdir -p "$subfam_workdir"
      cp "$tmpdir/subfam_input.fa" "$subfam_workdir/input.fa"

      (
        cd "$subfam_workdir"
        "$SUBFAM_BIN" "input.fa" 50
      ) || {
        log "    WARNING: SubFam failed for $sf_name"
        rm -rf "$tmpdir"
        continue
      }

      clw_file="$subfam_workdir/input.clw"
      if [[ -s "$clw_file" ]]; then
        cat "$tmpdir/consensus.fa" "$clw_file" > "$tmpdir/tierC_input.fa"
        run_mafft "$tmpdir/tierC_input.fa" "$sp_out/${sf_name}.subfam.fa" || {
          log "    WARNING: mafft failed for Tier C $sf_name"
        }
        n_tierC=$((n_tierC + 1))
      fi
    fi

    n_total=$((n_total + 1))
    rm -rf "$tmpdir"
  done

  # Clean up shared sear working directory for this species
  [[ -n "$sear_workdir" ]] && rm -rf "$sear_workdir"

  unset has_copies
done

# -- Summary -------------------------------------------------------------------
log "============================================================"
log "  DONE"
log "  Species:       $n_species"
log "  Subfamilies:   $n_total"
log "  Tier A (core):   $n_tierA"
log "  Tier B (best50): $n_tierB"
log "  Tier C (SubFam): $n_tierC"
log "  Tier D (evid.):  $n_tierD"
log "  Output: $OUT_BASE/"
log "============================================================"
