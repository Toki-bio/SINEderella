#!/usr/bin/env bash
# extract_subfam_only.sh
# Tier C (SubFam) generation only -- skips Tier A/B mafft runs.
# Run this when extract_alignments.sh already completed Tier A+B but
# SubFam was not available at that time.
#
# Usage:
#   extract_subfam_only.sh <consensus_bank.fa> <multi_run_dir> [output_dir]
#
# Dependencies: mafft, SubFam (auto-detected in SCRIPT_DIR or PATH)
set -euo pipefail

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ echo "ERROR: $*" >&2; exit 1; }

CONS_BANK="${1:-}"
MULTI_DIR="${2:-}"
OUT_BASE="${3:-}"

[[ -n "$CONS_BANK" ]] || die "Usage: $0 <consensus_bank.fa> <multi_run_dir> [output_dir]"
[[ -n "$MULTI_DIR"  ]] || die "Usage: $0 <consensus_bank.fa> <multi_run_dir> [output_dir]"
[[ -f "$CONS_BANK"  ]] || die "Consensus bank not found: $CONS_BANK"
[[ -d "$MULTI_DIR"  ]] || die "Multi-run directory not found: $MULTI_DIR"

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
[[ -n "$SUBFAM_BIN" ]] || die "SubFam not found in $SCRIPT_DIR or PATH"
log "SubFam: $SUBFAM_BIN"

# -- count_fasta / sample_fasta / extract_consensus helpers -------------------
count_fasta(){ grep -c '^>' "$1" 2>/dev/null || echo 0; }

sample_fasta(){
  local n="$1" input="$2" output="$3"
  # Pair header+seq on single line, sample n, then split back
  # Use shuf -n to avoid SIGPIPE killing script under set -eo pipefail
  awk '/^>/{if(seq)print hdr"\t"seq; hdr=$0; seq=""; next}{seq=seq $0} END{if(seq)print hdr"\t"seq}' "$input" \
    | shuf -n "$n" | awk -F'\t' '{print $1"\n"$2}' > "$output"
}

extract_consensus(){
  local name="$1" cons_bank="$2" out="$3"
  awk -v n="$name" '/^>/{found=0; if($1==">"n || $0==">"n) found=1} found{print}' "$cons_bank" > "$out"
  # Fallback: match first word after >
  if [[ ! -s "$out" ]]; then
    awk -v n="$name" 'BEGIN{f=0} /^>/{f=0; split($0,a," "); if(a[1]==">"n){f=1}} f{print}' "$cons_bank" > "$out"
  fi
}

# -- mafft wrapper ------------------------------------------------------------
run_mafft(){
  local input="$1" output="$2"
  mafft --thread "$THREADS" --threadtb "$THREADS" --threadit "$THREADS" \
        --localpair --maxiterate 1000 --ep 0.123 \
        --nuc --reorder --preservecase --quiet \
        "$input" > "$output" 2>/dev/null
}

# -- Main loop ----------------------------------------------------------------
n_species=0
n_tierC=0
n_skip_exists=0

for sp_dir in "$MULTI_DIR"/*/; do
  [[ -d "$sp_dir" ]] || continue
  sp="$(basename "$sp_dir")"
  [[ "$sp" == "summary" || "$sp" == "alignments" ]] && continue

  latest="$(ls -dt "$sp_dir"/run_* 2>/dev/null | head -n1 || true)"
  [[ -n "$latest" && -d "$latest" ]] || continue

  step2_out="$(ls -dt "$latest"/step2/step2_output* 2>/dev/null | head -n1 || true)"
  subfam_dir=""
  [[ -n "$step2_out" && -d "$step2_out/subfamilies" ]] && subfam_dir="$step2_out/subfamilies"
  [[ -n "$subfam_dir" ]] || continue

  n_species=$((n_species + 1))
  sp_out="$OUT_BASE/$sp"
  mkdir -p "$sp_out"

  log "[$sp] Scanning for Tier C candidates..."

  for sf_fasta in "$subfam_dir"/*.fasta; do
    [[ -f "$sf_fasta" ]] || continue
    sf_name="$(basename "$sf_fasta" .fasta)"
    nseqs=$(count_fasta "$sf_fasta")
    [[ $nseqs -gt 400 ]] || continue

    out_file="$sp_out/${sf_name}.subfam.fa"
    if [[ -f "$out_file" ]]; then
      log "  [$sf_name] SKIP (already exists: $out_file)"
      n_skip_exists=$((n_skip_exists + 1))
      continue
    fi

    log "  [$sf_name] $nseqs copies -> Tier C (SubFam)"

    tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/esf_${sp}_${sf_name}.XXXXXX")

    extract_consensus "$sf_name" "$CONS_BANK" "$tmpdir/consensus.fa"
    if [[ ! -s "$tmpdir/consensus.fa" ]]; then
      log "    WARNING: consensus not found for $sf_name, skipping"
      rm -rf "$tmpdir"
      continue
    fi

    if [[ $nseqs -le 2500 ]]; then
      cp "$sf_fasta" "$tmpdir/subfam_input.fa"
    else
      sample_fasta 2500 "$sf_fasta" "$tmpdir/subfam_input.fa"
    fi

    subfam_workdir="$tmpdir/subfam_work"
    mkdir -p "$subfam_workdir"
    cp "$tmpdir/subfam_input.fa" "$subfam_workdir/input.fa"

    (
      cd "$subfam_workdir"
      bash "$SUBFAM_BIN" "input.fa" 50
    ) || {
      log "    WARNING: SubFam failed for $sf_name, skipping"
      rm -rf "$tmpdir"
      continue
    }

    clw_file="$subfam_workdir/input.clw"
    if [[ -s "$clw_file" ]]; then
      cat "$tmpdir/consensus.fa" "$clw_file" > "$tmpdir/tierC_input.fa"
      run_mafft "$tmpdir/tierC_input.fa" "$out_file" || {
        log "    WARNING: mafft failed for Tier C $sf_name"
      }
      [[ -s "$out_file" ]] && n_tierC=$((n_tierC + 1))
      log "    -> $out_file"
    else
      log "    WARNING: SubFam produced no .clw for $sf_name"
    fi

    rm -rf "$tmpdir"
  done
done

log "============================================================"
log "  DONE"
log "  Species processed: $n_species"
log "  Tier C generated:  $n_tierC"
log "  Skipped (exists):  $n_skip_exists"
log "  Output: $OUT_BASE/"
log "============================================================"
