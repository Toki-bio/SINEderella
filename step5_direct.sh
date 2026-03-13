#!/usr/bin/env bash
# ==========================================================
# step5_direct.sh — Simple direct alignment for single runs
#
# Usage:
#   step5_direct.sh <run_dir>
#
# Outputs alignments directly to: <run_dir>/alignments/
# No symlinks, no wrapper folders, just results.
# ==========================================================
set -euo pipefail

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ echo "ERROR: $*" >&2; exit 1; }

# ── Arguments ────────────────────────────────────────────────────────────────
RUN_DIR="${1:-}"
[[ -n "$RUN_DIR" ]] || die "Usage: $0 <run_dir>"
[[ -d "$RUN_DIR" ]] || die "Run directory not found: $RUN_DIR"
RUN_DIR="$(readlink -f "$RUN_DIR")"

CONS_BANK="$RUN_DIR/consensuses.clean.fa"
[[ -f "$CONS_BANK" ]] || die "Consensus bank not found: $CONS_BANK"

# Find subfamilies directory
SUBFAM_DIR=""
for candidate in \
  "$RUN_DIR/step2/step2_output/subfamilies" \
  "$RUN_DIR/results/subfamilies"; do
  [[ -d "$candidate" ]] && SUBFAM_DIR="$candidate" && break
done
[[ -n "$SUBFAM_DIR" ]] || die "Subfamilies directory not found in run"

# Output directory
OUT_DIR="$RUN_DIR/alignments"
mkdir -p "$OUT_DIR"

log "Run: $RUN_DIR"
log "Subfamilies: $SUBFAM_DIR"
log "Output: $OUT_DIR"

# ── Helper functions ─────────────────────────────────────────────────────────
count_fasta(){ awk 'BEGIN{c=0}/^>/{c++}END{printf "%d\n",c}' "$1"; }

sample_fasta(){
  local n="$1" infile="$2" outfile="$3"
  awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next }
       { printf("%s",$0) }
       END  { printf("\n") }' "$infile" > "$outfile.lin"
  shuf "$outfile.lin" | head -n "$n" > "$outfile.sel" || true
  awk 'BEGIN{FS="\t"}{printf("%s\n%s\n",$1,$2)}' "$outfile.sel" > "$outfile"
  rm -f "$outfile.lin" "$outfile.sel"
}

# ── Process each subfamily ───────────────────────────────────────────────────
mapfile -t sf_files < <(find "$SUBFAM_DIR" -maxdepth 1 -name "*.fasta" -size +0c | sort)
nsf=${#sf_files[@]}
[[ $nsf -gt 0 ]] || die "No subfamily FASTAs found"

log "Found $nsf subfamilies"
n_aligned=0

for sf_fasta in "${sf_files[@]}"; do
  sf_name="$(basename "$sf_fasta" .fasta)"
  nseqs=$(count_fasta "$sf_fasta")
  [[ $nseqs -ge 1 ]] || continue

  log "  $sf_name ($nseqs copies)"

  tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/s5_${sf_name}.XXXXXX")
  trap "rm -rf '$tmpdir'" EXIT

  if [[ $nseqs -lt 400 ]]; then
    # Small: sample 200, cat with consensus, align
    sample_n=$((nseqs < 200 ? nseqs : 200))
    [[ $nseqs -le 200 ]] && cp "$sf_fasta" "$tmpdir/sampled.fa" || \
      sample_fasta "$sample_n" "$sf_fasta" "$tmpdir/sampled.fa"
    cat "$CONS_BANK" "$tmpdir/sampled.fa" > "$tmpdir/to_align.fa"
    mafft --thread "$(nproc 2>/dev/null || echo 1)" --nuc --auto --reorder --quiet \
          "$tmpdir/to_align.fa" > "$OUT_DIR/${sf_name}.aln.fa" 2>/dev/null || \
      log "    WARNING: mafft failed"
  else
    # Large: subsample 10000, use SubFam + mafft --add
    sample_n=$((nseqs < 10000 ? nseqs : 10000))
    [[ $nseqs -le 10000 ]] && cp "$sf_fasta" "$tmpdir/${sf_name}.fa" || \
      sample_fasta "$sample_n" "$sf_fasta" "$tmpdir/${sf_name}.fa"
    
    # Run SubFam if available
    if command -v SubFam >/dev/null 2>&1; then
      subfam_workdir="$tmpdir/subfam_work"
      mkdir -p "$subfam_workdir"
      cp "$tmpdir/${sf_name}.fa" "$subfam_workdir/${sf_name}.fa"
      (cd "$subfam_workdir" && SubFam "${sf_name}.fa") || {
        log "    WARNING: SubFam failed, falling back to direct mafft"
        cat "$CONS_BANK" "$tmpdir/${sf_name}.fa" > "$tmpdir/to_align.fa"
        mafft --thread "$(nproc 2>/dev/null || echo 1)" --nuc --auto --reorder --quiet \
              "$tmpdir/to_align.fa" > "$OUT_DIR/${sf_name}.aln.fa" 2>/dev/null || true
        rm -rf "$tmpdir"
        continue
      }
      
      clw_file="$subfam_workdir/${sf_name}.clw"
      [[ -s "$clw_file" ]] || {
        log "    WARNING: SubFam produced no .clw"
        rm -rf "$tmpdir"
        continue
      }
      
      # Convert .clw to FASTA if needed
      if head -1 "$clw_file" | grep -q "^CLUSTAL"; then
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
      
      # Add consensus with mafft
      mafft --add "$CONS_BANK" --thread "$(nproc 2>/dev/null || echo 1)" \
            --keeplength --reorder --quiet "$tmpdir/subfam_aln.fa" \
            > "$OUT_DIR/${sf_name}.aln.fa" 2>/dev/null || \
        log "    WARNING: mafft --add failed"
    else
      # Fallback if SubFam not available
      cat "$CONS_BANK" "$tmpdir/${sf_name}.fa" > "$tmpdir/to_align.fa"
      mafft --thread "$(nproc 2>/dev/null || echo 1)" --nuc --auto --reorder --quiet \
            "$tmpdir/to_align.fa" > "$OUT_DIR/${sf_name}.aln.fa" 2>/dev/null || \
        log "    WARNING: mafft failed"
    fi
  fi

  rm -rf "$tmpdir"
  [[ -f "$OUT_DIR/${sf_name}.aln.fa" ]] && n_aligned=$((n_aligned + 1))
done

log "============================================================"
log "Done! $n_aligned alignments in: $OUT_DIR/"
log "============================================================"
