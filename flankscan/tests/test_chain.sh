#!/usr/bin/env bash
# test_chain.sh - stage 6c on a toy with one planted three-unit chain (make_toy_chain.sh): the whole flankscan
# (fs_all.sh, stages 1-7) must find the pairs TD-TE and TE-TF, see that each candidate's end is open, build the
# chain TD + TE + TF once, accept it, and leave its ends closed. PASS/FAIL per check; exit 1 on a FAIL.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; FS="$(dirname "$HERE")"
T=${TOYC:-$HOME/tmp/fs_chain}; O=$T/fs
rm -rf "$T"; bash "$HERE/make_toy_chain.sh" "$T" >&2
bash "$FS/fs_all.sh" "$T" "$O" 8 50 2>&1 | grep -E "^fs|FAILED|==" >&2
fails=0
chk() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fails=$((fails+1)); fi; }
K() { gawk -F'\t' "NR>1 && ($1)" "$O/candidates.tsv" | wc -l; }
chk "chain: stage 4 finds the pairs TD-TE and TE-TF (>= 50 copies each)" \
    "[ \$(gawk -F'\t' 'NR>1 && \$6>=50 && ((\$1==\"TD\"&&\$3==\"TE\")||(\$1==\"TE\"&&\$3==\"TD\")||(\$1==\"TE\"&&\$3==\"TF\")||(\$1==\"TF\"&&\$3==\"TE\"))' $O/peaks.tsv | wc -l) -ge 2 ]"
chk "chain: each pair candidate has an OPEN end into the third unit (stage 6b)" \
    "[ \$(gawk -F'\t' 'NR>1 && \$1 ~ /^TD__TE/ && \$8 ~ /3/ && \$6==\"TF\"' $O/ends.tsv | wc -l) -ge 1 ] && [ \$(gawk -F'\t' 'NR>1 && \$1 ~ /^TE__TF/ && \$8 ~ /5/ && \$3==\"TD\"' $O/ends.tsv | wc -l) -ge 1 ]"
chk "chain: the pair candidates are extended (status extended:...)" \
    "[ \$(K '\$14 ~ /^extended:/') -ge 2 ]"
chk "chain: exactly one chain candidate is kept (both routes fold into one)" "[ \$(K '\$3==\"chain\" && \$14==\"kept\"') -eq 1 ]"
chk "chain: its consensus = the planted TD + linker + TE + linker + TF (best hit >= 95 % id, >= 95 % length)"     "PL=\$(seqkit fx2tab -n -l $T/planted.fa | cut -f2); ssearch36 -m 8 -E 1e-5 -z 11 -Z 1000 $T/planted.fa $O/candidates.fa 2> /dev/null | sort -t\$'	' -k12,12gr | head -1 | gawk -F'	' -v PL=\$PL '{ exit !(\$3 >= 95 && (\$8 - \$7 + 1) >= 0.95 * PL) }'"
chk "chain: hierarchy lists ONE element of three units TD, TE, TF with the two linkers" \
    "[ \$(gawk -F'\t' 'NR>1' $O/hierarchy.tsv | wc -l) -eq 1 ] && gawk -F'\t' 'NR==2{ n=split(\$18,p,\";\"); exit !(n==5 && p[1] ~ /^TD:1-1(4|5)/ && p[2]==\"gap:30\" && p[3] ~ /^TE:1-16/ && p[4]==\"gap:40\" && p[5] ~ /^TF:1-18[0-3]/)}' $O/hierarchy.tsv"
chk "chain: verdict accept (>= 70 % of its copies read as one full chain)" \
    "gawk -F'\t' 'NR==2{exit !(\$6==\"accept\" && \$5>=70)}' $O/hierarchy.tsv"
chk "chain: its own ends are closed" "gawk -F'\t' 'NR==2{exit !(\$15==\"-\")}' $O/hierarchy.tsv"
chk "chain: accepted.fa holds the chain, named from its units" \
    "[ \$(grep -c '^>' $O/accepted.fa) -eq 1 ] && grep -q '^>TD_TE_TF_C' $O/accepted.fa"
echo "fails=$fails"; exit $((fails > 0))
