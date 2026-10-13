#!/usr/bin/env bash
# fast_split.sh RUN FAMILIES.clean.tsv SUBFAMILY_BANK.fa [THREADS] - the second stage of the fast families mode.
#
# The run was made with the FAMILY consensuses as the bank (one query per family: fast, but biased towards the subfamilies the family
# consensus resembles). Here every family with several subfamilies has its firmly assigned copies assigned again, by the ordinary
# step 2, against the consensuses of its own subfamilies only. Subfamilies of different families never compete.
#
# Out: RUN/subfamily_split/<family>/ (a step 2 output dir each) and RUN/subfamily_split/subfamily_assignment.tsv
#      Sequence  Family  Subfamily  Bitscore  Votes  Status  Threshold
# Single-member families are copied through (subfamily = family). Copies of a family that were not firm at family level keep their
# family status and get Subfamily "NA".
set -euo pipefail
RUN=$(readlink -f "$1"); TAB=$(readlink -f "$2"); BANK=$(readlink -f "$3"); THREADS=${4:-${THREADS:-8}}
S2=$(ls -dt "$RUN"/step2/step2_output* 2>/dev/null | head -n1)
[[ -n "$S2" && -s "$S2/assignment_full.tsv" ]] || { echo "fast_split: no step 2 output in $RUN" >&2; exit 1; }
OUT="$RUN/subfamily_split"; mkdir -p "$OUT"; cd "$OUT" || exit 1
: > subfamily_assignment.tsv
printf 'Sequence\tFamily\tSubfamily\tBitscore\tVotes\tStatus\tThreshold\n' > subfamily_assignment.tsv
# family -> members
awk -F'\t' 'NR > 1 { print $2 "\t" $1 }' "$TAB" | sort -k1,1 > members.tsv
for fam in $(cut -f1 members.tsv | sort -u); do
    n=$(awk -F'\t' -v f="$fam" '$1 == f' members.tsv | wc -l)
    awk -F'\t' -v f="$fam" '$2 == f && $5 == "assigned" { print $1 }' "$S2/assignment_full.tsv" > "ids.$fam.txt"
    if (( n < 2 )); then
        sub=$(awk -F'\t' -v f="$fam" '$1 == f { print $2 }' members.tsv)
        awk -F'\t' -v f="$fam" -v s="$sub" 'BEGIN { OFS = "\t" } $2 == f { print $1, f, ($5 == "assigned" ? s : "NA"), $3, $4, $5, $6 }' "$S2/assignment_full.tsv" >> subfamily_assignment.tsv
        rm -f "ids.$fam.txt"
        continue
    fi
    mkdir -p "$fam"; cd "$fam" || exit 1
    awk -F'\t' -v f="$fam" '$1 == f { print $2 }' ../members.tsv > members.txt
    seqkit grep -f members.txt "$BANK" > members.fa
    seqkit grep -f "../ids.$fam.txt" "$RUN/genome.clean_step1/extracted.fasta" > copies.fa
    if [[ -s copies.fa ]]; then
        bash "$RUN/step2_asSINEment.sh" members.fa copies.fa "$THREADS" step2_output > step2.log 2>&1 \
            || { echo "fast_split: step 2 failed for family $fam (see $OUT/$fam/step2.log)" >&2; exit 1; }
        awk -F'\t' -v f="$fam" 'BEGIN { OFS = "\t" } NR > 1 { print $1, f, $2, $3, $4, $5, $6 }' step2_output/assignment_full.tsv >> ../subfamily_assignment.tsv
    fi
    cd ..
    # copies of this family that were not firm at family level
    awk -F'\t' -v f="$fam" 'BEGIN { OFS = "\t" } $2 == f && $5 != "assigned" { print $1, f, "NA", $3, $4, $5, $6 }' "$S2/assignment_full.tsv" >> subfamily_assignment.tsv
done
echo "fast_split: $(( $(wc -l < subfamily_assignment.tsv) - 1 )) copies, $(cut -f1 members.tsv | sort -u | wc -l) families -> $OUT/subfamily_assignment.tsv" >&2
