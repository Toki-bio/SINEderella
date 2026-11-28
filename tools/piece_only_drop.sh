#!/usr/bin/env bash
# piece_only_drop.sh RUN_ROOT
#
# Loci found only by a composite piece and not resolved by step 2 are not loci (owner 2026-11-28, docs/PARTIAL_COPIES.md).
#
# Step 1 searches every composite stretch with a unit as a query of its own (tools/composite_pieces.sh). A piece's hits that
# no whole bank consensus hit are kept only when step 2 assigned the copy, or when the family step made it a partial copy
# (family_view/partial_copies.tsv, Status "partial"). Every other such copy is taken out of step2_output/loci_verdict.tsv
# and listed in step2_output/piece_only_dropped.tsv. Atelerix bank10: the 69 bp tail of e2-3longB (281-349) found 108,429
# loci of its own, 129 assigned; rsi_final: 1,511 piece-only loci, 281 assigned, 771 partial (20 of the 21 r1 loci).
#
# Run by publish/verdict_chain.sh after tools/flank_verdict.py. PIECE_ONLY_DROP=0 skips it.
# Piece-only = the copy (coordinates from its name) overlaps no hit in <step1>/searches/gen-consensuses.clean.part_*.bed
# (whole bank consensuses, after the satellite screen). Needs bedtools.
set -euo pipefail
RUN=${1:?usage: piece_only_drop.sh RUN_ROOT}
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
[[ "${PIECE_ONLY_DROP:-1}" == "0" ]] && { log "piece-only drop: skipped (PIECE_ONLY_DROP=0)"; exit 0; }
S2=$(ls -dt "$RUN"/step2/step2_output* 2>/dev/null | head -n1 || true)
LV=$S2/loci_verdict.tsv
[[ -n "$S2" && -s "$LV" ]] || { echo "ERROR: piece_only_drop: no loci_verdict.tsv in $RUN/step2" >&2; exit 1; }
OUT=$S2/piece_only_dropped.tsv
# already applied to this loci_verdict.tsv (flank_verdict.py rewrites it, which makes it newer again): keep the list
if [[ -f "$OUT" && "$OUT" -nt "$LV" ]]; then log "PIECE_ONLY already applied ($OUT); loci_verdict.tsv unchanged"; exit 0; fi
shopt -s nullglob
SEARCHES=("$RUN"/*_step1/searches)
(( ${#SEARCHES[@]} == 1 )) || { echo "ERROR: piece_only_drop: expected one <step1>/searches in $RUN, found ${#SEARCHES[@]}" >&2; exit 1; }
PIECES=("${SEARCHES[0]}"/gen-pieces.part_*.bed)
WHOLE=("${SEARCHES[0]}"/gen-consensuses.clean.part_*.bed)
printf 'copy\tsubfamily\tassignment_status\treason\n' > "$OUT.tmp"
if (( ${#PIECES[@]} == 0 )); then
  mv "$OUT.tmp" "$OUT"; touch "$OUT"; log "PIECE_ONLY no pieces searched; loci_verdict.tsv unchanged"; exit 0
fi
(( ${#WHOLE[@]} > 0 )) || { echo "ERROR: piece_only_drop: no whole-consensus hit files in ${SEARCHES[0]}" >&2; exit 1; }
command -v bedtools > /dev/null || { echo "ERROR: piece_only_drop: bedtools not on PATH" >&2; exit 1; }
T=$(mktemp -d "${TMPDIR:-/tmp}/piece_only.XXXXXX"); trap 'rm -rf "$T"' EXIT

# copy name "chr:start-end(+)" -> bed (chr may hold ':'; start-end as in the name)
LC_ALL=C awk -F'\t' 'NR > 1 { s = $1; sub(/\([+-]\)$/, "", s); n = split(s, a, ":"); split(a[n], c, "-")
    print substr(s, 1, length(s) - length(a[n]) - 1) "\t" c[1] "\t" c[2] "\t" $1 }' "$LV" | LC_ALL=C sort -k1,1 -k2,2n > "$T/loci.bed"
cat "${WHOLE[@]}" | cut -f1-3 | LC_ALL=C sort -k1,1 -k2,2n > "$T/whole.bed"
bedtools intersect -v -sorted -a "$T/loci.bed" -b "$T/whole.bed" | cut -f4 > "$T/pieceonly.ids"

PART=$RUN/family_view/partial_copies.tsv
[[ -s "$PART" ]] || PART=/dev/null
# partial_copies.tsv: Sequence ... Status (column named Status); keep Status == partial
LC_ALL=C awk -F'\t' -v OUT="$OUT.tmp" -v LVNEW="$T/lv.new" '
  FILENAME == ARGV[1] { po[$1] = 1; next }
  FILENAME == ARGV[2] { if (FNR == 1) { for (i = 1; i <= NF; i++) if ($i == "Status") sc = i; next }
                        if (sc && $sc == "partial") part[$1] = 1; next }
  FNR == 1 { print > LVNEW; next }
  { if (($1 in po) && $3 != "assigned" && !($1 in part)) {
        print $1 "\t" $2 "\t" $3 "\tfound only by a composite piece; not assigned, not a partial copy" >> OUT; d++; next }
    if ($1 in po) { if ($3 == "assigned") ka++; else kp++ }
    print > LVNEW }
  END { printf "PIECE_ONLY piece_only=%d dropped=%d kept_assigned=%d kept_partial=%d\n", d + ka + kp, d, ka, kp > "/dev/stderr" }
' "$T/pieceonly.ids" "$PART" "$LV"
mv "$T/lv.new" "$LV.tmp" && mv "$LV.tmp" "$LV"
mv "$OUT.tmp" "$OUT"; touch "$OUT"     # newer than loci_verdict.tsv: marks the rule as applied
