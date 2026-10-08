#!/usr/bin/env bash
# SINEderella-route vs COSEG on a hand-labelled copy set. Run on KIT.
#   run_compare.sh SPECIES [WORKDIR]       SPECIES = saq | ccr | teu | dmo  (POS__SPECIES__* in aln_c)
# Steps: labelled copies -> reference = consensus of ALL copies (no subfamily is favoured) ->
# COSEG input (drop = Price's rule, keep = all copies) -> COSEG over -m sweep ->
# SubFam chunks on the same copies -> score everything against the owner's groups.
set -euo pipefail
SP=${1:?species}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
W=${2:-/data/V/toki/coseg_cmp/$SP}
ALN=${ALN:-/data/W/toki/SINE_disc/aln_c}
COSEG=${COSEG:-/data/V/toki/coseg/src}
SUBFAM=${SUBFAM:-$HERE/../../SubFam}
PY=${PYTHON:-python3.12}
THREADS=${THREADS:-16}
mkdir -p "$W"
cd "$W" || exit 1

$PY "$HERE/prep_labeled.py" "$ALN" "$SP" copies
N=$(grep -c '^>' copies.fa)

# reference: plurality consensus of an alignment of all copies, degapped
mafft --thread "$THREADS" --auto --quiet copies.fa > copies.aln
cons -sequence copies.aln -plurality $(( N / 4 )) -outseq ref.raw.fa -name ref -auto 2>/dev/null || \
    cons -filter -plurality $(( N / 4 )) -name ref < copies.aln > ref.raw.fa
awk '/^>/{print; next} {gsub(/[-nN]/,""); printf "%s", $0} END{print ""}' ref.raw.fa | awk 'NR==1{print;next}{printf "%s",$0} END{print ""}' > ref.fa
echo "reference: $(tail -n +2 ref.fa | tr -d '\n' | wc -c) bp"

for mode in drop keep; do
    $PY "$HERE/to_coseg.py" ref.fa copies.fa "cs_$mode" "$mode"
    tr 'acgt' 'ACGT' < ref.fa | awk 'NR==1{print;next}{printf "%s",$0} END{print ""}' > "cs_$mode.cons"
done

SPECS=()
for mode in drop keep; do
    for m in 20 10 5; do
        d="cs_${mode}_m$m"; rm -rf "$d"; mkdir "$d"
        ( cd "$d" || exit 1
          cp "../cs_$mode.seqs" "../cs_$mode.ins" "../cs_$mode.cons" .
          perl "$COSEG/runcoseg.pl" -k -d -m "$m" -c "cs_$mode.cons" -s "cs_$mode.seqs" -i "cs_$mode.ins" > run.log 2>&1 ) || echo "COSEG failed: $d (see $d/run.log)"
        if [ -f "$d/cs_$mode.seqs.assign" ]; then SPECS+=("$d/cs_$mode.seqs.assign@cs_$mode.names:COSEG_${mode}_m$m"); fi
    done
done

for n in 50 20 10; do
    d="sf_n$n"; rm -rf "$d"; mkdir "$d"
    ( cd "$d" || exit 1; bash "$SUBFAM/SubFam" ../copies.fa "$n" > log.txt 2>&1 ) || echo "SubFam failed: $d"
    [ -f "$d/copies.chunks.tsv" ] && SPECS+=("$d/copies.chunks.tsv:SubFam_n$n")
done

$PY "$HERE/score.py" copies.labels "${SPECS[@]}" | tee scores.txt
