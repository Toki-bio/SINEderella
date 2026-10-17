#!/usr/bin/env bash
# fs9b_partners.sh RUN_DIR [THREADS=8]
#
# Flankscan stage 9b: UNKNOWN PARTNERS from one-sided shared flanks. Runs after stage 9 (fs9_run.sh) on its output.
#
# Stage 3 finds a composite's partner only if the partner is in the bank (or in the fs0 partner library), and stage 4 needs
# >= 50 copies at a junction. A partner that is in no consensus therefore never shows. Stage 9 does see such copies: their
# flanks are similar from the junction outward, on ONE side, and their bodies are independent. Its body check (2026-10-14)
# then calls the pairs "flank_only" (the hot-spot signature) and drops them. Found on rsi (2026-10-17): 21 r6 / P26 copies of a
# review alignment shared ~300 bp upstream that matched nothing in the bank; 29 copies genome-wide carry it at the junction;
# fs9's flank-only pass had grouped them, the body check set them all to "unique".
#
# A group is reported when, within one family:
#   - its copies are linked by flank_only pairs (fs9: shared flank from the junction outward, bodies at family background) on
#     one side S;
#   - members that are duplicated loci are REMOVED first: a member that also shares the OTHER side with another member, or
#     that stage 9 calls a twin / array on any side (an insertion hot spot inside an older repeat and a segmental duplication
#     share both sides). A group-level share limit was wrong: it let rsi X + r6 keep 6 of 30 copies that were one
#     duplicated block on NC_142512 (both flanks, the body and an 87 bp insertion shared), owner's check 2026-10-17;
#   - what is left still has >= MING (10) copies on >= MINCTG (3) contigs;
# For each group: the partner consensus is built from the members' flanks on side S (MAFFT, plurality), and its extent is how far
# the members agree from the junction outward (>= AGREE 60 % of the bases, >= COVER 50 % of the members, in windows of 10);
# the partner is searched against the bank and the partner library (PARTNERS=fa, optional): a known name or "unknown"; the
# copies of ALL families that carry it at the junction (its junction-proximal PROBE (60) bp ending within JX (5) bp of the
# element, same orientation) are counted per family.
# Nothing is added to the bank: every group is a candidate for the user's judgement (verdict "check").
#
# In : RUN_DIR/results/flank_twins/<family>/{flank_groups.tsv,copies.bed} (stage 9), RUN_DIR/genome.clean.fa(.fai),
#      RUN_DIR/consensuses.clean.fa, RUN_DIR/results/assignment_full.tsv
# Out: RUN_DIR/results/unknown_partners.tsv  id family side copies contigs removed_duplicates partner_len partner_match carriers
#                                            carriers_by_family verdict
#      RUN_DIR/results/unknown_partners/<id>.members.txt, <id>.partner.fa, <id>.composite.fa (partner + family consensus, in
#      the element's orientation), <id>.review.aln.fa (members laid on the candidate composite, fs9b_review.py: anchored
#      columns, outer flanks lower case and packed against each member's own ends)
# Env: MING=10 MINCTG=3 LF=600 AGREE=0.6 COVER=0.5 PROBE=60 JX=5 MAXREV=60 PARTNERS=<fasta>
# Needs: gawk, bedtools, samtools, mafft, ssearch36.
set -euo pipefail
RUN=$(readlink -f "${1:?RUN_DIR}"); T=${2:-8}
MING=${MING:-10}; MINCTG=${MINCTG:-3}; LF=${LF:-600}; AGREE=${AGREE:-0.6}; COVER=${COVER:-0.5}
PROBE=${PROBE:-60}; JX=${JX:-5}; MAXREV=${MAXREV:-60}; PARTNERS=${PARTNERS:-}
G="$RUN/genome.clean.fa"; FT="$RUN/results/flank_twins"; AF="$RUN/results/assignment_full.tsv"; CONS="$RUN/consensuses.clean.fa"
OUT="$RUN/results/unknown_partners"
for f in "$G" "$AF" "$CONS"; do [[ -s "$f" ]] || { echo "fs9b: missing $f" >&2; exit 1; }; done
[[ -d "$FT" ]] || { echo "fs9b: no $FT (run stage 9, fs9_run.sh, first)" >&2; exit 1; }
[[ -s "$G.fai" ]] || samtools faidx "$G"
mkdir -p "$OUT"; export LC_ALL=C
W=$(mktemp -d "${TMPDIR:-/tmp}/fs9b.XXXXXX"); trap 'rm -rf "$W"' EXIT

# ---- 1. one-sided flank_only groups per family ------------------------------------------------------------------------------
: > "$W/groups.tsv"                                  # id family side member
for D in "$FT"/*/; do
    FAM=$(basename "$D")
    [[ -s "$D/flank_groups.tsv" && -s "$D/copies.bed" ]] || continue
    gawk -F'\t' -v OFS='\t' -v F="$FAM" -v MING="$MING" -v MINCTG="$MINCTG" '
        function root(x) { while (up[x] != x) x = up[x]; return x }
        NR == 1 { next }
        { a = $1; b = $2; s = $3
          pair[s, a, b] = 1; pair[s, b, a] = 1
          if ($7 == "twin" || $7 == "array") { dup[a] = 1; dup[b] = 1 }   # duplicated with its neighbourhood (stage 9)
          if ($7 == "flank_only") {                     # shared flank, independent bodies
              ka = s SUBSEP a; kb = s SUBSEP b
              if (!(ka in up)) up[ka] = ka; if (!(kb in up)) up[kb] = kb
              ra = root(ka); rb = root(kb); if (ra != rb) up[ra] = rb } }
        END {
          for (k in up) { r = root(k); n[r]++; mem[r] = mem[r] SUBSEP k }
          gid = 0
          for (r in n) {
              if (n[r] < MING) continue
              split(substr(mem[r], 2), m, SUBSEP)        # m[i] = side, copy alternating
              delete cs; side = m[1]; other = (side == 5) ? 3 : 5
              for (i = 1; i <= length(m); i += 2) cs[m[i + 1]] = 1
              # a member that also shares the OTHER side with another member, or that stage 9 calls a twin / array, is a
              # duplicated locus, not an independent carrier of the partner: it is removed, the group is judged on the rest
              # (rsi X + r6: 6 of 30 copies were one duplicated block on NC_142512, both flanks and the body shared)
              delete rm; nrm = 0
              for (c in cs) { if (c in dup) { rm[c] = 1; continue }
                  for (d in cs) if (d != c && ((other, c, d) in pair)) { rm[c] = 1; break } }
              for (c in rm) nrm++
              delete ctg; nc = 0; nk = 0
              for (c in cs) if (!(c in rm)) { nk++; split(c, q, ":"); if (!(q[1] in ctg)) { ctg[q[1]] = 1; nc++ } }
              if (nk < MING || nc < MINCTG) continue
              gid++
              for (c in cs) if (!(c in rm)) print F "_" side "p" gid, F, side, c, nc, nrm }
        }' "$D/flank_groups.tsv" >> "$W/groups.tsv"
done
NG=$(cut -f1 "$W/groups.tsv" | sort -u | wc -l)
echo "fs9b: $NG one-sided flank group(s) of >= $MING copies on >= $MINCTG contigs"

# ---- 2. flanks of every assigned copy (both sides, PROBE + 40 bp, element orientation) for the carrier count ------------------
gawk -F'\t' -v OFS='\t' 'NR > 1 && $5 == "assigned" && match($1, /^(.+):([0-9]+)-([0-9]+)\(([+-])\)$/, m) {
        print m[1], m[2], m[3], $1 "|" $2, 0, m[4] }' "$AF" | sort -k1,1 -k2,2n > "$W/all.bed"
FW=$((PROBE + 40))
bedtools flank -i "$W/all.bed" -g "$G.fai" -l $FW -r 0 -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null |
    sed 's/([+-])$//' > "$W/all5.fa"
bedtools flank -i "$W/all.bed" -g "$G.fai" -l 0 -r $FW -s | bedtools getfasta -fi "$G" -bed - -s -nameOnly 2> /dev/null |
    sed 's/([+-])$//' > "$W/all3.fa"
cat "$CONS" ${PARTNERS:+"$PARTNERS"} > "$W/known.fa"

printf "id\tfamily\tside\tcopies\tcontigs\tremoved_duplicates\tpartner_len\tpartner_match\tcarriers\tcarriers_by_family\tverdict\n" > "$OUT.tsv"
for ID in $(cut -f1 "$W/groups.tsv" | sort -u); do
    FAM=$(gawk -F'\t' -v I="$ID" '$1 == I { print $2; exit }' "$W/groups.tsv")
    S=$(gawk -F'\t' -v I="$ID" '$1 == I { print $3; exit }' "$W/groups.tsv")
    NC=$(gawk -F'\t' -v I="$ID" '$1 == I { print $5; exit }' "$W/groups.tsv")
    OS=$(gawk -F'\t' -v I="$ID" '$1 == I { print $6; exit }' "$W/groups.tsv")
    gawk -F'\t' -v I="$ID" '$1 == I { print $4 }' "$W/groups.tsv" | sort > "$OUT/$ID.members.txt"
    N=$(wc -l < "$OUT/$ID.members.txt")
    # ---- 3. partner: members' LF-bp flanks on side S, read in the element's orientation, aligned; extent from the junction -----
    gawk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] { m[$1]; next } ($4 in m)' "$OUT/$ID.members.txt" "$FT/$FAM/copies.bed" > "$W/m.bed"
    if [[ $S == 5 ]]; then bedtools flank -i "$W/m.bed" -g "$G.fai" -l "$LF" -r 0 -s > "$W/f.bed"
    else bedtools flank -i "$W/m.bed" -g "$G.fai" -l 0 -r "$LF" -s > "$W/f.bed"; fi
    bedtools getfasta -fi "$G" -bed "$W/f.bed" -s -nameOnly 2> /dev/null | sed 's/([+-])$//' > "$W/f.fa"
    mafft --quiet --thread "$T" --threadit 0 --auto "$W/f.fa" > "$W/f.aln" 2> /dev/null
    # plurality per column; walk from the junction (right end for a 5' flank, left end for a 3' flank) in windows of 10 columns
    gawk -v S="$S" -v AG="$AGREE" -v CV="$COVER" -v N="$N" '
        /^>/ { r++; next } { seq[r] = seq[r] toupper($0) }
        END { L = length(seq[1]); for (j = 1; j <= L; j++) { delete c; nb = 0
                  for (i = 1; i <= r; i++) { x = substr(seq[i], j, 1); if (x != "-") { c[x]++; nb++ } }
                  top = ""; tn = 0; for (x in c) if (c[x] > tn) { tn = c[x]; top = x }
                  base[j] = top; good[j] = (nb >= CV * N && tn >= AG * nb) ? 1 : 0; occ[j] = (nb >= CV * N) }
              # walk outward from the junction: stop at the first window of 10 occupied columns with < 7 good ones
              lo = 1; hi = L
              if (S == 5) { j = L; while (j >= 10) { g = 0; o = 0; for (k = j - 9; k <= j; k++) { g += good[k]; o += occ[k] }
                                if (o >= 5 && g < 0.7 * o) break; j -= 10 }; lo = j + 1 }
              else       { j = 1; while (j + 9 <= L) { g = 0; o = 0; for (k = j; k <= j + 9; k++) { g += good[k]; o += occ[k] }
                                if (o >= 5 && g < 0.7 * o) break; j += 10 }; hi = j - 1 }
              p = ""; for (j = lo; j <= hi; j++) if (occ[j] && base[j] != "") p = p base[j]
              print ">partner\n" p }' "$W/f.aln" > "$OUT/$ID.partner.fa"
    PL=$(gawk '!/^>/ { n += length($0) } END { print n + 0 }' "$OUT/$ID.partner.fa")
    MATCH="unknown"; CAR=0; BYF="-"
    if (( PL >= 30 )); then
        # ---- 4. known? partner against the bank (+ partner library) ----
        cp "$OUT/$ID.partner.fa" "$W/partner.fa"          # bare names: ssearch36 cannot open paths over ~120 characters
        ( cd "$W" && ssearch36 -m 8 -z 11 -E 1e-5 -T "$T" partner.fa known.fa 2> /dev/null ) > "$W/k.m8" || true
        MATCH=$(gawk -F'\t' -v L="$PL" 'NR == 1 && $4 >= 0.5 * L { printf "%s (%.0f %% over %d of %d bp)", $2, $3, $4, L; exit }' "$W/k.m8")
        [[ -n "$MATCH" ]] || MATCH="unknown"
        # ---- 5. carriers: the junction-proximal PROBE bp, ending within JX of the element, in any family ----
        if [[ $S == 5 ]]; then gawk -v P="$PROBE" '!/^>/ { s = s $0 } END { print ">probe\n" substr(s, length(s) - P + 1) }' "$OUT/$ID.partner.fa" > "$W/probe.fa"; LIB=all5.fa
        else gawk -v P="$PROBE" '!/^>/ { s = s $0 } END { print ">probe\n" substr(s, 1, P) }' "$OUT/$ID.partner.fa" > "$W/probe.fa"; LIB=all3.fa; fi
        ( cd "$W" && ssearch36 -3 -m 8 -z 11 -E 1e-3 -T "$T" probe.fa "$LIB" 2> /dev/null ) > "$W/c.m8" || true
        gawk -F'\t' -v S="$S" -v FW="$FW" -v JX="$JX" '
            !seen[$2]++ { e = (S == 5) ? $10 : $9; if ((S == 5 && e >= FW - JX) || (S == 3 && e <= 1 + JX)) { split($2, q, "|"); n[q[2]]++; t++ } }
            END { printf "%d\t", t + 0; o = ""; for (f in n) o = o sprintf("%s:%d,", f, n[f]); sub(/,$/, "", o); print (o == "" ? "-" : o) }' "$W/c.m8" > "$W/car.txt"
        CAR=$(cut -f1 "$W/car.txt"); BYF=$(cut -f2 "$W/car.txt")
    fi
    # ---- 6. candidate composite (element orientation) and a review alignment ----
    UC=$(gawk -v F="$FAM" '/^>/ { p = (substr($1, 2) == F); next } p' "$CONS" | tr -d '\n')
    PS=$(grep -v '^>' "$OUT/$ID.partner.fa" | tr -d '\n')
    if [[ $S == 5 ]]; then printf ">%s_composite partner(%d)+%s\n%s%s\n" "$ID" "$PL" "$FAM" "$PS" "$UC"; else printf ">%s_composite %s+partner(%d)\n%s%s\n" "$ID" "$FAM" "$PL" "$UC" "$PS"; fi > "$OUT/$ID.composite.fa"
    head -n "$MAXREV" "$OUT/$ID.members.txt" > "$W/rv.txt"
    # members cut with LF + 100 bp on the partner side and 100 bp on the other side, in the element's orientation; laid on the
    # candidate composite by fs9b_review.py (anchored columns, outer flanks lower case and packed against each member's ends)
    gawk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] { m[$1]; next } ($4 in m)' "$W/rv.txt" "$FT/$FAM/copies.bed" |
        gawk -F'\t' -v OFS='\t' -v S="$S" -v LF="$LF" '{ up = (S == 5) ? LF + 100 : 100; dn = (S == 5) ? 100 : LF + 100
            if ($6 == "+") { $2 = ($2 > up ? $2 - up : 0); $3 = $3 + dn } else { $2 = ($2 > dn ? $2 - dn : 0); $3 = $3 + up }; print }' > "$W/rv.bed"
    bedtools getfasta -fi "$G" -bed "$W/rv.bed" -s -nameOnly 2> /dev/null | sed 's/([+-])$//' > "$W/rv.fa"
    python3 -I "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/fs9b_review.py" "$OUT/$ID.composite.fa" "$W/rv.fa" "$OUT/$ID.review.aln.fa" "$T" > /dev/null 2>> "$W/review.err" ||
        echo "fs9b: review alignment of $ID failed (see the run's stderr)" >&2
    printf "%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\tcheck\n" "$ID" "$FAM" "$S" "$N" "$NC" "$OS" "$PL" "$MATCH" "$CAR" "$BYF" >> "$OUT.tsv"
done
# ---- 7. summary: drop groups without a common partner, merge the same partner found in several families, family specificity ----
# family specificity of a group: own_share = carriers in its own family / all carriers; enrichment = own_share / the family's
# share of all copies. "family-specific" = own_share >= 0.5 and enrichment >= 2 (a big family cannot reach more than 1 / its
# share: rsi P1 holds 22 % of the copies, so at most 4.5); a sequence next to copies of every family in proportion to their
# size (rsi: a ~120 bp sequence next to ~1 % of the copies of every family) stays near 1: "background".
gawk -F'\t' '{ split($4, q, "|"); n[q[2]]++; t++ } END { for (f in n) print f "\t" n[f]; print "__all__\t" t }' "$W/all.bed" > "$W/famsize.tsv"
: > "$W/partners.fa"
gawk -F'\t' -v M=30 'NR > 1 && $7 >= M { print $1 }' "$OUT.tsv" | while read -r I; do
    printf ">%s\n%s\n" "$I" "$(grep -v '^>' "$OUT/$I.partner.fa" | tr -d '\n')" >> "$W/partners.fa"; done
( cd "$W" && ssearch36 -m 8 -z -1 -E 1 -T "$T" partners.fa partners.fa 2> /dev/null ) > "$W/pp.m8" || true
gawk -F'\t' -v OFS='\t' -v MIN=30 '
    FILENAME == ARGV[1] { if (/^>/) { id = substr($1, 2) } else { len[id] += length($0) } ; next }
    FILENAME == ARGV[2] { size[$1] = $2; next }
    FILENAME == ARGV[3] { if ($1 != $2 && $3 >= 80 && $4 >= 0.6 * (len[$1] > len[$2] ? len[$1] : len[$2])) link[$1, $2] = 1; next }
    function root(x) { while (up[x] != x) x = up[x]; return x }
    FNR == 1 { next }
    $7 >= MIN { row[$1] = $0; up[$1] = $1; ids[++n] = $1 }
    END {
      for (k in link) { split(k, p, SUBSEP); if ((p[1] in up) && (p[2] in up)) { a = root(p[1]); b = root(p[2]); if (a != b) up[a] = b } }
      print "cluster", "groups", "families", "copies_in_groups", "partner_len", "partner_match", "carriers", "own_share", "enrichment", "specific_family", "kind"
      for (i = 1; i <= n; i++) { r = root(ids[i]); members[r] = members[r] "," ids[i] }
      for (r in members) {
          split(substr(members[r], 2), g, ","); fams = ""; cop = 0; best = 0; bf = "-"; L = 0; car = 0; mt = "unknown"; bsh = 0
          for (j in g) { split(row[g[j]], f, "\t"); fams = fams (index("," fams ",", "," f[2] ",") ? "" : f[2] ","); cop += f[4]
              if (f[7] > L) L = f[7]; if (f[9] > car) car = f[9]; if (f[8] != "unknown") mt = f[8]
              # specificity of this group: carriers in its own family vs expected from family size
              own = 0; nb = split(f[10], bb, ","); for (k = 1; k <= nb; k++) { split(bb[k], kv, ":"); if (kv[1] == f[2]) own = kv[2] }
              if (f[9] > 0 && size[f[2]] > 0) { sh = own / f[9]; sp = sh / (size[f[2]] / size["__all__"]); if (sh >= 0.5 && sp > best) { best = sp; bf = f[2]; bsh = sh } } }
          sub(/,$/, "", fams)
          kind = (mt != "unknown") ? "known partner" : (best >= 2 ? "family-specific: composite candidate" : "all families: background")
          printf "%s\t%s\t%s\t%d\t%d\t%s\t%d\t%.2f\t%.1f\t%s\t%s\n", r, substr(members[r], 2), fams, cop, L, mt, car, bsh, best, bf, kind } }' \
    "$W/partners.fa" "$W/famsize.tsv" "$W/pp.m8" "$OUT.tsv" | { read -r h; echo "$h"; sort -t$'\t' -k11,11 -k9,9gr; } > "$OUT.summary.tsv"
echo "fs9b: $(( $(wc -l < "$OUT.tsv") - 1 )) candidate partner group(s) -> $OUT.tsv"
echo "fs9b: $(( $(wc -l < "$OUT.summary.tsv") - 1 )) distinct partner(s) with >= 30 bp -> $OUT.summary.tsv"
column -t -s$'\t' "$OUT.summary.tsv" 2> /dev/null || cat "$OUT.summary.tsv"
exit 0
