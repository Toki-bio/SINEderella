#!/usr/bin/env bash
# fs7_hierarchy.sh OUT_DIR
#
# Stage 7 of flankscan: the hierarchy of the elements for the report - which kept candidate is built
# from which parts, where each part starts and ends in its own consensus, the gap between them, how
# many copies carry the layout and the stage-6 verdict. Also writes the accepted candidates under short
# bank names, ready for `SINEderella --add` (adding them stays the user's decision).
#
# In : OUT/candidates.tsv, OUT/peaks.tsv, OUT/reassign.tsv (stage 6), OUT/cons.tsv (stage 3),
#      OUT/ends.tsv (stage 6b), OUT/chains.tsv (stage 6c), OUT/candidates.fa
# Out: OUT/hierarchy.tsv   element bank_name type copies pct_full verdict cons_len
#                          part1 part1_start part1_end gap part2 part2_start part2_end
#                          open_ends open5_unit open3_unit parts
#                          parts = the whole element, any number of units, in its orientation:
#                          unit:start-end entries (positions in the unit's own consensus) joined by ";",
#                          with gap:N entries for linkers / unexplained sequence between units, e.g.
#                          r1:1-154;gap:39;r3:1-201. part1 / part2 = its first two units (older readers).
#                          bank_name = the units joined by "_" + the peak id, family suffixes like _58seqs
#                          dropped; open_ends = 5 / 3 / 53 / - (stage 6b / 6c: the end still continues into
#                          a bank unit, open5_unit / open3_unit); a candidate accepted by stage 6 with an
#                          open end gets verdict "open" and is left out of accepted.fa. Candidates that
#                          were extended into a chain (status extended:...) are not listed - the chain is.
#      OUT/accepted.fa     the accepted candidates, renamed to bank_name
set -euo pipefail
OUT=${1:?OUT_DIR}
cd "$OUT"
EF=ends.tsv; [[ -s "$EF" ]] || EF=/dev/null      # stage 6b output; without it no end is checked
CF=chains.tsv; [[ -s "$CF" ]] || CF=/dev/null    # stage 6c output
gawk -F'\t' -v OFS='\t' '
    function short(f) { sub(/_[0-9]+seqs$/, "", f); return f }
    FILENAME == ARGV[1] { len[$1] = $2; next }                                   # cons.tsv: name len
    FILENAME == ARGV[2] { if (FNR > 1) pk[$15] = $0; next }                      # peaks.tsv by peak id
    FILENAME == ARGV[3] { if (FNR > 1) { pf[$1] = $5; vd[$1] = $6 }; next }      # reassign.tsv
    FILENAME == ARGV[4] { if (FNR > 1) { eo[$1] = $8; e5[$1] = $3; e3[$1] = $6 }; next }   # ends.tsv
    FILENAME == ARGV[5] { if (FNR > 1) { cparent[$1] = $2; cu5[$1] = $4; cq5[$1] = $5 "-" $6; cg5[$1] = $7
                                         cu3[$1] = $8; cq3[$1] = $9 "-" $10; cg3[$1] = $11 }; next }   # chains.tsv
    FNR == 1 { print "element", "bank_name", "type", "copies", "pct_full", "verdict", "cons_len",
                     "part1", "part1_start", "part1_end", "gap", "part2", "part2_start", "part2_end",
                     "open_ends", "open5_unit", "open3_unit", "parts"; next }
    {
        # the parts of every candidate (kept or not: a chain builds on its parent)
        if ($3 == "chain") {
            P = parts[cparent[$1]]
            if (cu5[$1] != "-") P = short(cu5[$1]) ":" cq5[$1] ((cg5[$1] + 0 > 0) ? ";gap:" cg5[$1] : "") ";" P
            if (cu3[$1] != "-") P = P ((cg3[$1] + 0 > 0) ? ";gap:" cg3[$1] : "") ";" short(cu3[$1]) ":" cq3[$1]
        } else {
            split(pk[$2], p, "\t")                    # family side partner rel type n n_group exp main_j p_j gap
            if (p[2] == 3) { u = p[1]; ue = p[9];  d = p[3]; ds = p[10] }     # the copy is the upstream part
            else           { u = p[3]; ue = p[10]; d = p[1]; ds = p[9] }
            gap = (p[11] > 0) ? p[11] : 0
            P = short(u) ":1-" ue (gap > 0 ? ";gap:" gap : "") ";" short(d) ":" ds "-" len[d]
        }
        parts[$1] = P
        if ($14 != "kept") next
        np = split(P, pp, ";"); nu = 0; bank = ""
        for (i = 1; i <= np; i++) if (pp[i] !~ /^gap:/) { nu++; split(pp[i], q, ":"); unit[nu] = q[1]; span[nu] = q[2]; bank = bank (bank == "" ? "" : "_") q[1]
                                                             if (nu == 1) g1 = (pp[i + 1] ~ /^gap:/) ? substr(pp[i + 1], 5) : 0 }
        bank = bank "_" $2
        split(span[1], s1, "-"); split(span[2], s2, "-")
        vdd = ($1 in vd ? vd[$1] : ".")
        op = ($1 in eo ? eo[$1] : "-")
        if (vdd == "accept" && op != "-") vdd = "open"        # an end still continues into a bank unit: not closed
        print $1, bank, $3, $6, ($1 in pf ? pf[$1] : "."), vdd, $13,
              unit[1], s1[1], s1[2], g1, unit[2], s2[1], s2[2], op,
              ((op ~ /5/) ? short(e5[$1]) : "-"), ((op ~ /3/) ? short(e3[$1]) : "-"), P
    }' cons.tsv peaks.tsv reassign.tsv "$EF" "$CF" candidates.tsv > hierarchy.tsv
# accepted candidates under their bank names
gawk -F'\t' 'FNR == NR { if (FNR > 1 && $6 == "accept") nm[$1] = $2; next }
             /^>/ { n = substr($1, 2); keep = (n in nm); if (keep) print ">" nm[n]; next }
             keep' hierarchy.tsv candidates.fa > accepted.fa
echo "fs7: $(($(wc -l < hierarchy.tsv) - 1)) elements in $OUT/hierarchy.tsv, $(grep -c '^>' accepted.fa) accepted -> $OUT/accepted.fa" >&2
