#!/usr/bin/env bash
# ==========================================================
# step5_align_subfamilies.sh
#
# For each species in a multi-run project, align subfamily
# copies against the consensus bank.
#
# Usage:
#   step5_align_subfamilies.sh <consensi_bank.fa> <multi_run_dir>
#
# For each subfamily FASTA in each species' subfamilies/ folder:
#   - <400 copies  → sample 200, cat with consensi bank, mafft align
#   - >=400 copies → subsample 10000, run SubFam, convert .clw to
#                     aligned FASTA, then mafft --add consensi bank
#                     (skips redundant full re-alignment)
#
# Output:
#   <multi_run_dir>/alignments/<Species>/<subfamily>.aln.fa
# ==========================================================
set -euo pipefail

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ echo "ERROR: $*" >&2; exit 1; }

# ── Arguments ────────────────────────────────────────────────────────────────
CONS_BANK="${1:-}"
MULTI_DIR="${2:-}"

[[ -n "$CONS_BANK" ]] || die "usage: $0 <consensi_bank.fa> <multi_run_dir>"
[[ -n "$MULTI_DIR" ]] || die "usage: $0 <consensi_bank.fa> <multi_run_dir>"

[[ -f "$CONS_BANK" ]] || die "consensus bank not found: $CONS_BANK"
[[ -d "$MULTI_DIR" ]] || die "multi-run directory not found: $MULTI_DIR"

CONS_BANK="$(readlink -f "$CONS_BANK")"
MULTI_DIR="$(readlink -f "$MULTI_DIR")"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "$SCRIPT_DIR/SubFam" ]]; then
  SUBFAM_BIN="$SCRIPT_DIR/SubFam"
  chmod +x "$SUBFAM_BIN" 2>/dev/null || true
elif command -v SubFam >/dev/null 2>&1; then
  SUBFAM_BIN="$(command -v SubFam)"
else
  die "SubFam not found in $SCRIPT_DIR or in PATH"
fi

command -v mafft >/dev/null 2>&1 || die "mafft not in PATH"
command -v shuf  >/dev/null 2>&1 || die "shuf not in PATH"

ALIGN_DIR="$MULTI_DIR/alignments"
mkdir -p "$ALIGN_DIR"

# ── Helper: count sequences in a FASTA ──────────────────────────────────────
count_fasta(){
  awk 'BEGIN{c=0}/^>/{c++}END{printf "%d\n",c}' "$1"
}

# ── Helper: sample N random FASTA sequences from a bank ─────────────────────
# Adapted from user's custom sampler
sample_fasta(){
  local n="$1" infile="$2" outfile="$3"
  # Linearise FASTA to one-seq-per-line, shuffle, take N
  awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next }
       { printf("%s",$0) }
       END  { printf("\n") }' "$infile" > "$outfile.lin"
  shuf "$outfile.lin" | head -n "$n" > "$outfile.sel" || true
  awk 'BEGIN{FS="\t"}{printf("%s\n%s\n",$1,$2)}' "$outfile.sel" > "$outfile"
  rm -f "$outfile.lin" "$outfile.sel"
}

###############################################################################
# Discover species directories
###############################################################################
log "Consensus bank: $CONS_BANK"
log "Multi-run dir:  $MULTI_DIR"

n_species=0
n_aligned=0

for sp_dir in "$MULTI_DIR"/*/; do
  [[ -d "$sp_dir" ]] || continue
  sp="$(basename "$sp_dir")"
  [[ "$sp" == "summary" || "$sp" == "alignments" ]] && continue

  # Find latest run_* directory
  latest="$(ls -dt "$sp_dir"/run_* 2>/dev/null | head -n1 || true)"
  [[ -n "$latest" && -d "$latest" ]] || continue

  # Find subfamilies directory (step2 output or results)
  subfam_dir=""
  for candidate in \
    "$(ls -dt "$latest"/step2/step2_output* 2>/dev/null | head -n1 || true)/subfamilies" \
    "$latest/results/subfamilies"; do
    if [[ -d "$candidate" ]]; then
      subfam_dir="$candidate"
      break
    fi
  done
  [[ -n "$subfam_dir" && -d "$subfam_dir" ]] || {
    log "[$sp] SKIP — no subfamilies directory found"
    continue
  }

  n_species=$((n_species + 1))
  sp_align_dir="$ALIGN_DIR/$sp"
  mkdir -p "$sp_align_dir"

  log "[$sp] subfamilies: $subfam_dir"

  mapfile -t sf_files < <(find "$subfam_dir" -maxdepth 1 -name "*.fasta" -size +0c | sort)
  nsf=${#sf_files[@]}
  [[ $nsf -gt 0 ]] || { log "[$sp] SKIP — no non-empty subfamily FASTAs"; continue; }

  isf=0
  for sf_fasta in "${sf_files[@]}"; do
    sf_name="$(basename "$sf_fasta" .fasta)"
    isf=$((isf + 1))
    nseqs=$(count_fasta "$sf_fasta")
    [[ $nseqs -ge 1 ]] || continue

    log "  [$isf/$nsf] $sf_name ($nseqs copies)"

    # Create a temporary working directory for this subfamily
    tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/s5_${sp}_${sf_name}.XXXXXX")
    trap_cleanup(){ rm -rf "$tmpdir"; }

    if [[ $nseqs -lt 400 ]]; then
      # ── Small subfamily: sample 200 (or all if fewer), cat + align ──
      sample_n=$((nseqs < 200 ? nseqs : 200))
      log "    Sampling $sample_n copies, cat with consensi, aligning..."

      if [[ $nseqs -le 200 ]]; then
        cp "$sf_fasta" "$tmpdir/sampled.fa"
      else
        sample_fasta "$sample_n" "$sf_fasta" "$tmpdir/sampled.fa"
      fi

      # Concatenate: consensi bank + sampled copies
      cat "$CONS_BANK" "$tmpdir/sampled.fa" > "$tmpdir/to_align.fa"

      # Align
      mafft --thread "$(nproc 2>/dev/null || echo 1)" \
            --threadtb "$(nproc 2>/dev/null || echo 1)" \
            --threadit "$(nproc 2>/dev/null || echo 1)" \
            --nuc --auto --reorder --quiet \
            "$tmpdir/to_align.fa" \
            > "$sp_align_dir/${sf_name}.aln.fa" 2>/dev/null || {
        log "    WARNING: mafft failed for $sf_name"
        rm -rf "$tmpdir"
        continue
      }

    else
      # ── Large subfamily: subsample 10000, run SubFam, use --add ─────
      sample_n=$((nseqs < 10000 ? nseqs : 10000))
      log "    Subsampling $sample_n copies for SubFam..."

      if [[ $nseqs -le 10000 ]]; then
        cp "$sf_fasta" "$tmpdir/${sf_name}.fa"
      else
        sample_fasta "$sample_n" "$sf_fasta" "$tmpdir/${sf_name}.fa"
      fi

      # Run SubFam in a dedicated subfolder
      subfam_workdir="$tmpdir/subfam_work"
      mkdir -p "$subfam_workdir"
      cp "$tmpdir/${sf_name}.fa" "$subfam_workdir/${sf_name}.fa"
      (
        cd "$subfam_workdir"
        "$SUBFAM_BIN" "${sf_name}.fa"
      ) || {
        log "    WARNING: SubFam failed for $sf_name"
        rm -rf "$tmpdir"
        continue
      }

      clw_file="$subfam_workdir/${sf_name}.clw"
      [[ -s "$clw_file" ]] || {
        log "    WARNING: SubFam produced no .clw for $sf_name"
        rm -rf "$tmpdir"
        continue
      }

      # Convert .clw to aligned FASTA (preserve gaps) if Clustal format;
      # if already FASTA, use as-is.
      if head -1 "$clw_file" | grep -q "^CLUSTAL"; then
        log "    Converting .clw from Clustal to aligned FASTA..."
        awk '
          /^CLUSTAL/      { next }
          /^$/            { next }
          /^[[:space:]]/  { next }
          {
            if (NF >= 2) {
              name = $1; sequence = $2
              if (!(name in seq)) order[++count] = name
              seq[name] = seq[name] sequence
            }
          }
          END {
            for (i = 1; i <= count; i++)
              printf ">%s\n%s\n", order[i], seq[order[i]]
          }' "$clw_file" > "$tmpdir/subfam_aln.fa"
      else
        cp "$clw_file" "$tmpdir/subfam_aln.fa"
      fi

      log "    Adding consensi bank to SubFam alignment (mafft --add)..."
      mafft --add "$CONS_BANK" \
            --thread "$(nproc 2>/dev/null || echo 1)" \
            --keeplength \
            --reorder --quiet \
            "$tmpdir/subfam_aln.fa" \
            > "$sp_align_dir/${sf_name}.aln.fa" 2>/dev/null || {
        log "    WARNING: mafft --add failed for $sf_name"
        rm -rf "$tmpdir"
        continue
      }
    fi

    # Clean up temp dir for this subfamily
    rm -rf "$tmpdir"
    n_aligned=$((n_aligned + 1))
    log "    -> $sp_align_dir/${sf_name}.aln.fa"
  done
done

log "============================================================"
log "Done. $n_aligned alignments across $n_species species."
log "Output: $ALIGN_DIR/"
log "============================================================"
