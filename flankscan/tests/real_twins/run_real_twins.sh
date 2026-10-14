#!/usr/bin/env bash
# Real-data test of the fs9 twin rules on two multiplying squamate SINEs (2026-10-14):
#   Squam1  in Anolis carolinensis, UCSC anoCar2 (soft-masked, 36 % lower case): 79,225 full-length copies
#   Squam3C in Python bivittatus, NCBI RefSeq GCF_000186305.1 (soft-masked, 24 %): 6,713 full-length copies
# Consensuses: Squam1 and Squam3C from Vassetzky, Kosushkin & Ryskov 2023, Mobile DNA 14:21, Additional files 9 and 2.
#   run_real_twins.sh WORKDIR [THREADS]      needs curl, blastn/makeblastdb, bedtools, samtools, jellyfish, ssearch36, mafft, gawk
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); FS=$HERE/../..; W=${1:?WORKDIR}; T=${2:-4}; mkdir -p "$W/dl"; cd "$W"
[ -s dl/anoCar2.fa.gz ] || curl -sSf -o dl/anoCar2.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/anoCar2/bigZips/anoCar2.fa.gz
U=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/186/305/GCF_000186305.1_Python_molurus_bivittatus-5.0.2
[ -s dl/pyt.fna.gz ] || curl -sSf -o dl/pyt.fna.gz $U/GCF_000186305.1_Python_molurus_bivittatus-5.0.2_genomic.fna.gz
[ -s aca.fa ] || zcat dl/anoCar2.fa.gz > aca.fa
[ -s pyt.fa ] || zcat dl/pyt.fna.gz | sed 's/^>\([^ ]*\).*/>\1/' > pyt.fa
for g in aca pyt; do [ -s $g.fa.fai ] || samtools faidx $g.fa; [ -s db_$g.nsq ] || [ -s db_$g.00.nsq ] || makeblastdb -in $g.fa -dbtype nucl -out db_$g > /dev/null; done
run() {   # genome tag consensus length
    local g=$1 t=$2 c=$3 L=$4
    blastn -task blastn -query "$c" -db db_$g -evalue 1e-10 -word_size 11 -max_target_seqs 1000000 -max_hsps 1000000 -num_threads "$T" \
        -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' | sed 's/\tref|\([^|]*\)|/\t\1/' > $t.hits.tsv
    # full-length copies: >= 80 % of the consensus, >= 65 % identity; overlapping hits of one strand merged
    gawk -F'\t' -v OFS='\t' -v L=$L '{ if (($6-$5+1)/L < 0.8 || $3 < 65) next; s = ($7<$8) ? $7 : $8; e = ($7<$8) ? $8 : $7
        print $2, s-1, e, "h", $10, ($7<$8) ? "+" : "-" }' $t.hits.tsv | sort -k1,1 -k2,2n | bedtools merge -s -i - -c 5,6 -o max,distinct |
        gawk -v OFS='\t' -v T=$t '{ print $1, $2, $3, T "_" NR, 0, $5 }' > $t.copies.bed
    # genome-wide counts of the flank 20-mers only (jellyfish --if): the fs9 repeat filter needs nothing else
    bedtools flank -i $t.copies.bed -g $g.fa.fai -b 100 -s | bedtools getfasta -fi $g.fa -bed - > $t.flanks.fa
    jellyfish count -m 20 -C -s 50M -t "$T" --if $t.flanks.fa -o $t.jf $g.fa; jellyfish dump -c -L 21 $t.jf > $t.k20.txt
    for m in "1:85:70" "0:85:70" "0:95:95"; do IFS=: read sm i1 i2 <<< "$m"
        SOFTMASK=$sm ID1=$i1 ID2=$i2 KCOUNT=$PWD/$t.k20.txt bash "$FS/fs9_twins.sh" $g.fa $t.copies.bed fs9_${t}_sm${sm}_$i1 > fs9_${t}_sm${sm}_$i1.log 2>&1
        python3 -I "$HERE/scripts/pair_body_id.py" $g.fa $t.copies.bed fs9_${t}_sm${sm}_$i1/flank_groups.tsv body_${t}_sm${sm}_$i1.tsv
    done
    python3 -I "$HERE/scripts/random_pairs.py" $t.copies.bed 300 1 > rand_$t.tsv
    python3 -I "$HERE/scripts/pair_body_id.py" $g.fa $t.copies.bed rand_$t.tsv body_${t}_rand.tsv
    python3 "$HERE/scripts/summ.py" body_${t}_rand.tsv strict body_${t}_sm1_85.tsv "only:masked95" body_${t}_sm0_95.tsv "only:SOFTMASK=0" body_${t}_sm0_85.tsv | tee summary_$t.txt
}
run aca sq1 "$HERE/Squam1.fa" 329
run pyt sq3 "$HERE/Squam3C.fa" 247
