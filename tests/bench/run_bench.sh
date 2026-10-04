#!/usr/bin/env bash
# run_bench.sh BENCHDIR - time MAFFT modes on each case (base / +600 flanks); quality = Q vs the current
# mode (L-INS-i, 1 thread) and the continuation.py call on the result. Output: BENCHDIR/results.tsv
source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
B=$1; T=${THREADS:-16}; DISC=~/SINE_discriminator/site
COMMON="--nuc --reorder --preservecase --adjustdirection --quiet"
declare -A M=(
  [linsi1]="--localpair --maxiterate 1000 --ep 0.123 $COMMON"
  [linsi16]="--localpair --maxiterate 1000 --ep 0.123 --thread $T $COMMON"
  [linsi_it2]="--localpair --maxiterate 2 --ep 0.123 --thread $T $COMMON"
  [fftnsi2]="--retree 2 --maxiterate 2 --thread $T $COMMON"
  [fftns2]="--retree 2 --thread $T $COMMON"
)
ORDER="linsi1 linsi16 linsi_it2 fftnsi2 fftns2"
echo -e "case\tflanks\tmethod\tseconds\tQ_vs_linsi1\tcont5\tcont3" > $B/results.tsv
for d in $B/*/; do c=$(basename $d)
  for tag in base long; do
    for m in $ORDER; do
      s=$(date +%s.%N)
      mafft ${M[$m]} $d/$tag.fa > $d/$tag.$m.aln 2>/dev/null
      e=$(date +%s.%N); sec=$(echo "$e - $s" | bc)
      q=$(python3 ~/tmp/bench/qscore.py $d/$tag.linsi1.aln $d/$tag.$m.aln)
      ct=$(python3 - "$d/$tag.$m.aln" <<'PY'
import sys; sys.path.insert(0, "/home/toki/SINE_discriminator/site")
from fix_alignments import read_fa
import continuation as C
n, s = read_fa(sys.argv[1])
ci = next((i for i, x in enumerate(n) if "CONSENSUS" in x), 0)
e = C.extents(n, s, ci)
print("%s:%s\t%s:%s" % (e["5"]["status"], e["5"]["bp"], e["3"]["status"], e["3"]["bp"]))
PY
)
      printf "%s\t%s\t%s\t%.1f\t%s\t%s\n" $c $tag $m $sec $q "$ct" >> $B/results.tsv
    done
  done
done
echo DONE >> $B/results.tsv
