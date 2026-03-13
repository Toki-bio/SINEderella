#!/usr/bin/env bash
set -euo pipefail

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }

RUN_ROOT="${1:-}"
[[ -n "$RUN_ROOT" ]] || die "usage: $0 RUN_ROOT [THREADS]"
RUN_ROOT="$(readlink -f "$RUN_ROOT")"

# Threads: optional second arg or auto-detect
THREADS="${2:-$(nproc 2>/dev/null || echo 1)}"

# latest step2 output dir
OUT="$(ls -dt "$RUN_ROOT"/step2/step2_output* 2>/dev/null | head -n1 || true)"
[[ -n "$OUT" && -d "$OUT" ]] || die "cannot find $RUN_ROOT/step2/step2_output*"

ASSIGN="$OUT/assignment_full.tsv"
EXTRACTED="$RUN_ROOT/genome.clean_step1/extracted.fasta"
CONS="$RUN_ROOT/consensuses.clean.fa"
HITS_LABELED="$RUN_ROOT/genome.clean_step1/searches/all_hits.labeled.bed"
CONFLICT_BED="$RUN_ROOT/genome.clean_step1/searches/regions.by_subfam.bed"

log "RUN_ROOT=$RUN_ROOT"
log "OUT=$OUT"
log "ASSIGN=$ASSIGN"
log "EXTRACTED=$EXTRACTED"
log "CONS=$CONS"
log "HITS_LABELED=$HITS_LABELED"
log "CONFLICT_BED=$CONFLICT_BED"

# hard requirements
[[ -s "$ASSIGN" ]]    || die "missing/empty: $ASSIGN"
[[ -s "$EXTRACTED" ]] || die "missing/empty: $EXTRACTED"
[[ -s "$CONS" ]]      || die "missing/empty: $CONS"
[[ -s "$CONFLICT_BED" ]] || die "missing/empty: $CONFLICT_BED"
command -v bedtools  >/dev/null 2>&1 || die "bedtools not in PATH"
command -v ssearch36 >/dev/null 2>&1 || die "ssearch36 not in PATH"
command -v seqkit    >/dev/null 2>&1 || die "seqkit not in PATH"

SELF_BITS="$OUT/self_bits.tsv"
REAL_SELF_BITS="$OUT/self_bits_real.tsv"
SIM_SCORES="$OUT/sim_scores.tsv"
ALL="$OUT/all_sines.bedlike.ALL.tsv"
SUMMARY="$OUT/summary.by_subfam.tsv"

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/postprocess_tmp.XXXXXX")
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"' EXIT

###############################################################################
# 0) self_bits.tsv (minimal): derive per-subfam max(thr) from assignment_full
#    (does NOT run ssearch36; uses what's already in ASSIGN)
###############################################################################
if [[ ! -s "$SELF_BITS" ]]; then
  log "Building $SELF_BITS (from ASSIGN: max thr per subfam)"
  awk 'BEGIN{FS="[ \t]+"; OFS="\t"}
       NF>=3{
         sf=$2; thr=$3+0
         if(!(sf in mx) || thr>mx[sf]) mx[sf]=thr
       }
       END{for(sf in mx) print sf,mx[sf]}' "$ASSIGN" \
  | sort -k1,1 > "$SELF_BITS"
else
  log "Keeping existing $SELF_BITS"
fi

###############################################################################
# 0b) Single ssearch36: real self_bits + per-copy similarity score
#     - self_bits_real.tsv:  subfam -> self-alignment bitscore
#     - sim_scores.tsv:      seqID  -> sim_bs  self_bs  sim_ratio
###############################################################################
if [[ ! -s "$SIM_SCORES" ]]; then
  log "Computing single ssearch36 similarity scores (threads=$THREADS)..."

  # --- A) Real self_bits: ssearch36 consensus vs itself ---
  log "  Running consensus self-alignment for real self_bits..."
  ssearch36 -Q -n -z 11 -E 100 -T "$THREADS" -m 8 \
    "$CONS" "$CONS" > "$tmpdir/self.m8" 2>/dev/null || true

  awk -F'\t' '{
    q=$1; s=$2; bs=$12+0
    if(q==s && (!(q in mx) || bs>mx[q])) mx[q]=bs
  } END {
    for(c in mx) printf "%s\t%.1f\n", c, mx[c]
  }' "$tmpdir/self.m8" | sort -k1,1 > "$REAL_SELF_BITS"

  log "  Self-bits computed for $(wc -l < "$REAL_SELF_BITS") subfamilies"

  # --- B) Per-copy single alignment grouped by subfamily ---
  awk -F'\t' 'NR>1 && $5=="assigned"{print $1 "\t" $2}' "$ASSIGN" \
    > "$tmpdir/assigned_map.tsv"

  awk -F'\t' '{print $2}' "$tmpdir/assigned_map.tsv" \
    | sort -u > "$tmpdir/sim_subfam_list.txt"

  nsf=$(wc -l < "$tmpdir/sim_subfam_list.txt")
  log "  Running per-copy ssearch36 for $nsf subfamilies..."

  : > "$tmpdir/all_sim.m8"
  isf=0

  while IFS= read -r sf; do
    [[ -n "$sf" ]] || continue
    isf=$((isf+1))

    # Extract single consensus for this subfamily
    awk -v sf="$sf" '
      /^>/{name=$1; sub(/^>/,"",name); found=(name==sf)}
      found{print}
    ' "$CONS" > "$tmpdir/cons_single.fa"

    # Collect copy IDs assigned to this subfamily
    awk -F'\t' -v sf="$sf" '$2==sf{print $1}' "$tmpdir/assigned_map.tsv" \
      > "$tmpdir/copy_ids.txt"

    ncopies=$(wc -l < "$tmpdir/copy_ids.txt")
    [[ $ncopies -gt 0 ]] || continue

    # Extract copies from EXTRACTED fasta
    seqkit grep -f "$tmpdir/copy_ids.txt" "$EXTRACTED" \
      > "$tmpdir/copies.fa" 2>/dev/null || true

    # Single ssearch36: consensus query vs copies library
    if [[ -s "$tmpdir/copies.fa" && -s "$tmpdir/cons_single.fa" ]]; then
      ssearch36 -Q -n -z 11 -E 100 -T "$THREADS" -m 8 \
        "$tmpdir/cons_single.fa" "$tmpdir/copies.fa" \
        >> "$tmpdir/all_sim.m8" 2>/dev/null || true
    fi

    printf "  %d/%d subfamilies (%s: %d copies)\r" "$isf" "$nsf" "$sf" "$ncopies" >&2
  done < "$tmpdir/sim_subfam_list.txt"
  echo >&2

  # Parse: best bitscore per copy from all subfamilies
  awk -F'\t' '{
    seq=$2; bs=$12+0
    if(!(seq in mx) || bs>mx[seq]) mx[seq]=bs
  } END {
    for(s in mx) printf "%s\t%.1f\n", s, mx[s]
  }' "$tmpdir/all_sim.m8" | sort -k1,1 > "$tmpdir/per_copy_bs.tsv"

  # Build sim_scores.tsv: seqID, single_bs, self_bits, sim_ratio
  awk -F'\t' 'BEGIN{OFS="\t"}
    FILENAME==ARGV[1]{sb[$1]=$2; next}
    FILENAME==ARGV[2]{sf[$1]=$2; next}
    FILENAME==ARGV[3]{
      seq=$1; bs=$2+0
      subfam=(seq in sf ? sf[seq] : "")
      self=(subfam!="" && subfam in sb ? sb[subfam]+0 : 0)
      ratio=(self>0 ? bs/self : 0)
      printf "%s\t%.1f\t%.1f\t%.4f\n", seq, bs, self, ratio
    }
  ' "$REAL_SELF_BITS" "$tmpdir/assigned_map.tsv" "$tmpdir/per_copy_bs.tsv" \
    > "$SIM_SCORES"

  log "  Similarity scores computed for $(wc -l < "$SIM_SCORES") assigned copies"
else
  log "Keeping existing $SIM_SCORES"
fi

###############################################################################
# 1) Build base ALL table from EXTRACTED fasta headers + ASSIGN
#    Output cols (12):
#      1 chr
#      2 start0
#      3 end
#      4 subfam_from_extracted (for trace)
#      5 assigned_subfam (from ASSIGN if present else extracted)
#      6 strand
#      7 best_bitscore (prefer ASSIGN col6 else extracted header 3rd token)
#      8 leak_flag ('.' for now)
#      9 conflict_flag ('.' for now)
#     10 note (status/votes/thr + any runner* if present)
#     11 sim_bitscore (single ssearch36 vs assigned consensus, or '.')
#     12 sim_ratio   (sim_bitscore / self_bits_real, or '.')
###############################################################################
log "Building $ALL (base join)"
all_tmp="$tmpdir/all.base.tsv"

# Ensure SIM_SCORES exists for the awk join (may be empty if ssearch36 produced no hits)
[[ -f "$SIM_SCORES" ]] || : > "$SIM_SCORES"

awk '
  BEGIN{FS="[ \t]+"; OFS="\t"}

  # FILE 1: read SIM_SCORES (seqID \t sim_bs \t self_bs \t sim_ratio)
  FILENAME==ARGV[1]{
    sim_bs[$1]=$2
    sim_ratio[$1]=$4
    next
  }

  # FILE 2: read ASSIGN (tab or space); keep whole line tail if it contains runner tags
  FILENAME==ARGV[2]{
    id=$1
    asf[id]=$2
    athr[id]=$3
    avotes[id]=$4
    ast[id]=$5
    abs[id]=($6==""? "" : $6)

    # if assignment line itself contains runner tags (rare), capture them from $0
    r=""; rbs=""; rr=""
    if(match($0,/runner=[^; \t]+/)){
      tmp=substr($0,RSTART,RLENGTH); sub(/^runner=/,"",tmp); r=tmp
    }
    if(match($0,/runner_bs=[0-9.]+/)){
      tmp=substr($0,RSTART,RLENGTH); sub(/^runner_bs=/,"",tmp); rbs=tmp
    }
    if(match($0,/runner_ratio=[0-9.]+/)){
      tmp=substr($0,RSTART,RLENGTH); sub(/^runner_ratio=/,"",tmp); rr=tmp
    }
    if(r!="" || rbs!="" || rr!=""){
      arunner[id]=r
      arunner_bs[id]=rbs
      arunner_ratio[id]=rr
    }
    next
  }

  # FILE 3: parse EXTRACTED fasta headers
  # >ctg:start-end(strand)|Subfam|Bits
  /^>/{
    h=substr($0,2)

    n=split(h,p,"|")
    loc=p[1]             # ctg:start-end(strand)
    sf_from_hdr=(n>=2?p[2]:"")
    bs_from_hdr=(n>=3?p[3]:"")

    # contig
    ctg=loc
    sub(/:.*/,"",ctg)

    # range part
    rng=loc
    sub(/^[^:]+:/,"",rng)        # start-end(strand)
    sub(/\([+-]\).*/,"",rng)     # strip strand suffix if present
    split(rng,xy,"-")
    s=xy[1]+0
    e=xy[2]+0

    # strand
    strand="."
    if(match(loc,/\([+-]\)/)){
      tmp=substr(loc,RSTART,RLENGTH)
      gsub(/[()]/,"",tmp)
      strand=tmp
    }

    # bed start0/end
    s0=s-1
    if(s0<0) s0=0

    assigned=(loc in asf ? asf[loc] : sf_from_hdr)
    thr=(loc in athr ? athr[loc] : "")
    votes=(loc in avotes ? avotes[loc] : "")
    st=(loc in ast ? ast[loc] : "NA")

    # best bitscore: prefer ASSIGN col6, else header 3rd token
    bs=""
    if(loc in abs && abs[loc]!="") bs=abs[loc]
    else bs=bs_from_hdr

    note="status="st";votes="votes";thr="thr

    # carry runner tags if present in ASSIGN line
    if(loc in arunner)      note=note";runner="arunner[loc]
    if(loc in arunner_bs)   note=note";runner_bs="arunner_bs[loc]
    if(loc in arunner_ratio)note=note";runner_ratio="arunner_ratio[loc]

    # sim cols: single ssearch36 bitscore and ratio
    sbs=(loc in sim_bs ? sim_bs[loc] : ".")
    srt=(loc in sim_ratio ? sim_ratio[loc] : ".")

    print ctg,s0,e,sf_from_hdr,assigned,strand,bs,".",".",note,sbs,srt
  }
' "$SIM_SCORES" "$ASSIGN" "$EXTRACTED" > "$all_tmp"

# sanity: avoid clobbering output with empty
[[ -s "$all_tmp" ]] || die "base build produced empty: $all_tmp"
mv -f "$all_tmp" "$ALL"

###############################################################################
# 2) LEAK annotation from runner_ratio>=0.90 in col10
#    (marks col8=LEAK when threshold met)
###############################################################################
log "Annotating LEAK (runner_ratio>=0.90) from col10 runner_ratio=..."
leak_tmp="$tmpdir/all.leak.tsv"

awk -F'\t' 'BEGIN{OFS="\t"}
  {
    rr=""
    if(match($10,/runner_ratio=[0-9.]+/)){
      tmp=substr($10,RSTART,RLENGTH); sub(/^runner_ratio=/,"",tmp); rr=tmp+0
    }
    if(rr!="" && rr>=0.90) $8="LEAK"
    print
  }' "$ALL" > "$leak_tmp"

[[ -s "$leak_tmp" ]] || die "LEAK step produced empty: $leak_tmp"
mv -f "$leak_tmp" "$ALL"

###############################################################################
# 3) CONFLICT annotation using regions.by_subfam.bed
#    Only rows whose intersected B4 contains comma are considered "conflict".
#    IMPORTANT: B has 4 cols => in -wa -wb output, B4 is field 8.
###############################################################################
log "Annotating CONFLICT using $CONFLICT_BED"

# stable row id = NR
loci_bed="$tmpdir/loci.nr.bed"
awk -F'\t' 'BEGIN{OFS="\t"} NF>=3{print $1,$2,$3,NR}' "$ALL" > "$loci_bed"
[[ -s "$loci_bed" ]] || die "loci.nr.bed is empty (ALL may be empty?)"

# make conflict map: row_id -> comma-subfam-list
conf_map="$tmpdir/loci.conflict.map.tsv"
bedtools intersect -wa -wb -a "$loci_bed" -b "$CONFLICT_BED" \
| awk 'BEGIN{FS=OFS="\t"}
       {
         id=$4
         sublist=$8   # B4 is field 8 (A4 + B4)
         if(index(sublist,",")>0) conf[id]=sublist
       }
       END{for(id in conf) print id,conf[id]}' > "$conf_map"

# apply conflict annotations (idempotent: won.t append conf_subfams twice)
conf_tmp="$tmpdir/all.conf.tsv"
if [[ -s "$conf_map" ]]; then
  awk -F'\t' 'BEGIN{OFS="\t"}
    NR==FNR{conf[$1]=$2; next}
    {
      id=FNR
      if(id in conf){
        sublist=conf[id]
        $9="CONFLICT"

        # conf_alt=yes if any token != assigned subfam ($5)
        calt="no"
        n=split(sublist,a,",")
        for(i=1;i<=n;i++) if(a[i]!=$5){calt="yes"; break}

        if($10 !~ /(^|;)conf_subfams=/){
          $10=$10 ";conf_subfams=" sublist ";conf_alt=" calt
        }
      }
      print
    }' "$conf_map" "$ALL" > "$conf_tmp"
else
  # No conflicts: just copy ALL as is
  cp "$ALL" "$conf_tmp"
fi

[[ -s "$conf_tmp" ]] || die "CONFLICT step produced empty: $conf_tmp"
mv -f "$conf_tmp" "$ALL"

###############################################################################
# 3.5) UNASSIGNED.TSV: soft-assign unassigned sequences via original search query
#    For each sequence NOT status=assigned, look up all_hits.labeled.bed to find
#    the original search query/queries that produced this hit. If multiple,
#    pick the one with the highest search score (col5 in labeled bed).
###############################################################################
UNASSIGNED_TSV="$OUT/unassigned.tsv"
ALL_VOTES="$OUT/all_votes.tsv"

log "Building $UNASSIGNED_TSV (soft assignment via original search query)"

# 1) Extract unassigned loci from ALL table as BED (not status=assigned)
unassigned_bed="$tmpdir/unassigned.bed"
awk -F'\t' 'BEGIN{OFS="\t"}
  $10 !~ /(^|;)status=assigned(;|$)/{
    print $1, $2, $3, NR, $6
  }' "$ALL" > "$unassigned_bed"

unassigned_n=$(wc -l < "$unassigned_bed")
log "  Unassigned loci to soft-assign: $unassigned_n"

if [[ $unassigned_n -gt 0 && -s "$HITS_LABELED" ]]; then
  # 2) Intersect unassigned loci with all_hits.labeled.bed to find original queries
  #    HITS_LABELED is 6-col: chr start end subfam score strand
  #    We use -wa -wb to get both sides
  soft_intersect="$tmpdir/unassigned.intersect.tsv"
  bedtools intersect -wa -wb -a "$unassigned_bed" -b "$HITS_LABELED" \
    > "$soft_intersect" 2>/dev/null || true

  # 3) For each unassigned sequence, pick the query with highest score
  #    A cols: chr(1) start(2) end(3) rowNR(4) strand(5)
  #    B cols: chr(6) start(7) end(8) subfam(9) score(10) strand(11)
  soft_best="$tmpdir/unassigned.soft_best.tsv"
  awk 'BEGIN{FS=OFS="\t"}{
    id = $4   # row NR from ALL
    subfam = $9
    score = $10 + 0
    if(!(id in best_score) || score > best_score[id]){
      best_score[id] = score
      best_subfam[id] = subfam
    }
    # also collect all queries for this locus
    if(id in all_subfams){
      # check if already listed
      if(index(all_subfams[id], subfam) == 0)
        all_subfams[id] = all_subfams[id] "," subfam
    } else {
      all_subfams[id] = subfam
    }
  }
  END{
    for(id in best_subfam){
      print id, best_subfam[id], best_score[id], all_subfams[id]
    }
  }' "$soft_intersect" | sort -k1,1n > "$soft_best"

  # 4) Build unassigned.tsv by joining ALL table, soft_best, and all_votes.tsv
  {
    printf "SeqID\tChr\tStart\tEnd\tStrand\tSoft_Subfamily\tSearch_Score\tAll_Queries\tReason\tssearch_Best_Subfamily\tssearch_Votes\tssearch_Bitscore\n"

    # Read soft_best into a lookup, then all_votes into a lookup, then process ALL
    # Use separate passes to avoid NR==FNR issues with empty files
    awk -F'\t' -v softfile="$soft_best" -v votefile="$ALL_VOTES" '
    BEGIN{
      OFS="\t"
      # Read soft_best
      while((getline line < softfile) > 0){
        n=split(line,a,"\t")
        soft_sf[a[1]] = a[2]
        soft_score[a[1]] = a[3]
        soft_all[a[1]] = a[4]
      }
      close(softfile)
      # Read all_votes.tsv
      if(votefile != ""){
        while((getline line < votefile) > 0){
          n=split(line,a,"\t")
          vote_cons[a[1]] = a[2]
          vote_n[a[1]] = a[3]
          vote_bs[a[1]] = a[4]
        }
        close(votefile)
      }
    }
    {
      nr = FNR
      if($10 !~ /(^|;)status=assigned(;|$)/){
        # Determine reason
        reason = "unknown"
        if($10 ~ /(^|;)status=no_unanimous/)       reason = "no_unanimous_votes"
        else if($10 ~ /(^|;)status=rejected_low_bitscore/) reason = "rejected_low_bitscore"
        else if($10 ~ /(^|;)status=NA/)             reason = "no_ssearch_data"

        # Sequence ID as chr:start-end(strand)
        s1 = $2 + 1  # convert bed start0 back to 1-based
        strand = ($6 == "." ? "+" : $6)
        seqid = $1 ":" s1 "-" $3 "(" strand ")"

        # Soft assignment from original search
        ssf = (nr in soft_sf ? soft_sf[nr] : ".")
        ssc = (nr in soft_score ? soft_score[nr] : ".")
        sal = (nr in soft_all ? soft_all[nr] : ".")

        # ssearch36 vote data
        vc = (seqid in vote_cons ? vote_cons[seqid] : ".")
        vn = (seqid in vote_n ? vote_n[seqid] : ".")
        vb = (seqid in vote_bs ? vote_bs[seqid] : ".")

        printf "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
          seqid, $1, s1, $3, strand, ssf, ssc, sal, reason, vc, vn, vb
      }
    }' "$ALL"
  } > "$UNASSIGNED_TSV"

  log "  Wrote $(( $(wc -l < "$UNASSIGNED_TSV") - 1 )) entries to $UNASSIGNED_TSV"

else
  # No unassigned or no labeled hits: create empty with header
  printf "SeqID\tChr\tStart\tEnd\tStrand\tSoft_Subfamily\tSearch_Score\tAll_Queries\tReason\tssearch_Best_Subfamily\tssearch_Votes\tssearch_Bitscore\n" \
    > "$UNASSIGNED_TSV"
  log "  No unassigned sequences to process"
fi

###############################################################################
# 4) summary.by_subfam.tsv
#    assigned_total: status=assigned
#    leak_n: among assigned, col8==LEAK
#    conf_alt_n: among assigned, col9==CONFLICT AND conf_alt=yes
#    sim_ratio: mean/median from col12 for assigned copies
###############################################################################
log "Building $SUMMARY"
sum_tmp="$tmpdir/summary.tmp.tsv"

# Count total extracted sequences
total_extracted=$(awk '/^>/{c++} END{print c+0}' "$EXTRACTED")

awk -F'\t' -v OFS="\t" -v utfile="$UNASSIGNED_TSV" -v total_ext="$total_extracted" '
  BEGIN{
    soft_total=0
    # Read unassigned.tsv to count soft assignments per subfamily
    while((getline line < utfile) > 0){
      n=split(line,a,"\t")
      if(a[1]=="SeqID") continue
      sf = a[6]
      if(sf != ".") {
        soft[sf]++
        soft_total++
      }
    }
    close(utfile)
  }
  {
    if($10 ~ /(^|;)status=assigned(;|$)/){
      tot[$5]++
      if($8=="LEAK") leak[$5]++
      if($9=="CONFLICT" && $10 ~ /(^|;)conf_alt=yes(;|$)/) conf[$5]++

      # Collect sim_ratio values (col12) for assigned copies
      if(NF>=12 && $12!="." && $12+0>0){
        sim_count[$5]++
        sim_sum[$5] += $12+0
        # store values for median
        sim_vals[$5, sim_count[$5]] = $12+0
        if(!($5 in sim_min) || $12+0 < sim_min[$5]) sim_min[$5] = $12+0
        if(!($5 in sim_max) || $12+0 > sim_max[$5]) sim_max[$5] = $12+0
      }
    }
    total_all++
  }
  END{
    for(k in tot) names[k]=1
    for(k in soft) names[k]=1

    print "subfam","firm_assigned","soft_assigned","total_assigned","leak_n","conf_alt_n","firm_pct","total_pct","leak_pct","conf_alt_pct","sim_mean","sim_median"
    for(k in names){
      f = (k in tot ? tot[k] : 0)
      s = (k in soft ? soft[k] : 0)
      t = f + s
      l = (k in leak ? leak[k] : 0)
      c = (k in conf ? conf[k] : 0)
      fp = (total_ext > 0 ? 100*f/total_ext : 0)
      tp = (total_ext > 0 ? 100*t/total_ext : 0)
      lp = (f > 0 ? 100*l/f : 0)
      cp = (f > 0 ? 100*c/f : 0)

      # similarity ratio stats
      sm = "."
      smed = "."
      if(k in sim_count && sim_count[k]>0){
        sm = sprintf("%.4f", sim_sum[k]/sim_count[k])
        # sort values for median (simple insertion sort)
        nc = sim_count[k]
        for(i=1; i<=nc; i++) sorted[i] = sim_vals[k, i]
        for(i=2; i<=nc; i++){
          v = sorted[i]
          j = i-1
          while(j>=1 && sorted[j]>v){ sorted[j+1]=sorted[j]; j-- }
          sorted[j+1] = v
        }
        if(nc%2==1) smed = sprintf("%.4f", sorted[int(nc/2)+1])
        else        smed = sprintf("%.4f", (sorted[nc/2]+sorted[nc/2+1])/2)
      }

      printf "%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%s\n", \
        k,f,s,t,l,c,fp,tp,lp,cp,sm,smed
    }
  }' "$ALL" > "$sum_tmp"

( head -n1 "$sum_tmp"; tail -n +2 "$sum_tmp" | sort -k2,2nr ) > "$SUMMARY"
rm -f "$sum_tmp"

###############################################################################
# 5) sanity
###############################################################################
runner_ratio_tags=$(awk -F'\t' '$10 ~ /(^|;)runner_ratio=/{c++} END{print c+0}' "$ALL")
leak_rows=$(awk -F'\t' '$8=="LEAK"{c++} END{print c+0}' "$ALL")
conf_rows=$(awk -F'\t' '$9=="CONFLICT"{c++} END{print c+0}' "$ALL")
conf_alt_yes_any=$(awk -F'\t' '$10 ~ /(^|;)conf_alt=yes(;|$)/{c++} END{print c+0}' "$ALL")
sim_scored=$(awk -F'\t' 'NF>=12 && $12!="." && $12+0>0{c++} END{print c+0}' "$ALL")

log "Done."
log "Outputs:"
log "  $SELF_BITS"
log "  $REAL_SELF_BITS"
log "  $SIM_SCORES"
log "  $ALL"
log "  $UNASSIGNED_TSV"
log "  $SUMMARY"
log "[SANITY] runner_ratio_tags=$runner_ratio_tags"
log "[SANITY] LEAK_rows=$leak_rows"
log "[SANITY] CONFLICT_rows=$conf_rows"
log "[SANITY] conf_alt_yes_any=$conf_alt_yes_any"
log "[SANITY] sim_ratio_scored=$sim_scored"
log "[SANITY] unassigned_soft_assigned=$(awk -F'\t' 'NR>1 && $6!="."{c++} END{print c+0}' "$UNASSIGNED_TSV")"
