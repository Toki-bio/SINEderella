#!/usr/bin/env bash
# test_fs9.sh [WORKDIR] - toy genome with every flank case planted (make_toy_flank.py), fs9 run, scored against the truth
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; W=${1:-${TMPDIR:-/tmp}/fs9_toy}
rm -rf "$W"; mkdir -p "$W"
python3 "$HERE/make_toy_flank.py" "$W/toy" 1
bash "$HERE/../fs9_twins.sh" "$W/toy/genome.fa" "$W/toy/copies.bed" "$W/out" 4 > "$W/run.log" 2>&1 || { tail -20 "$W/run.log"; exit 1; }
gawk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] { cls[$1] = $2; grp[$1] = $3; id[$1] = $4; next }
    FILENAME == ARGV[2] { st[$1] = $2; next }
    { }
    END {
        fail = 0
        for (c in cls) { k = cls[c]; s = st[c]
            if (k == "U" || k == "HOT" || k == "LC" || k == "R") { n[k]++; if (s ~ /twin/) { bad[k]++; fail++ } }
            else if (k == "END") { n[k]++; if (s != "untestable") { bad[k]++; fail++ } }
            else if (k == "ARR") { n[k]++; if (s != "array") { bad[k]++; fail++ } }
            else if (k == "S") { g = grp[c]; gn[g]++; gid[g] = id[c]; if (s == "twin1" || s == "twin2") gok[g]++; gs[g] = gs[g] " " s } }
        for (k in n) printf "%-4s copies %3d  wrong %d\n", k, n[k], bad[k] + 0
        for (g in gn) { want = (gid[g] >= 85) ? "twin1" : "twin2"; printf "%-3s identity %s%%  found %d/%d  status:%s\n", g, gid[g], gok[g] + 0, gn[g], gs[g] }
        print (fail == 0 ? "NOTE: no false classification" : "FALSE CLASSIFICATIONS: " fail) }' "$W/toy/truth.tsv" "$W/out/copy_status.tsv" "$W/out/copy_status.tsv" | sort
