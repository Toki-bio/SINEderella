#!/usr/bin/env bash
# Oma orient repair: orient bank -> step1 sear -> step2-4 -> publish alignments.
#
# Hits inherit query direction from step1; this is the authoritative strand fix.
# Usage: bash launch_oma_rerun_orient_step1.sh [RUN_ROOT]
set -euo pipefail

RUN="${1:-/staging/tmp/scorpions/oma/run_oma}"
SINEDERELLA="${SINEDERELLA:-/staging/tmp/SINEderella}"
DISC="${DISC:-/staging/tmp/sinedisc}"
THREADS="${THREADS:-32}"
CHUNK_BP="${CHUNK_BP:-30000}"
FLANK="${FLANK:-50}"
LOG="${RUN}/rerun_orient_step1.log"
STAMP="$(date +%Y%m%d_%H%M)"

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG"; }
die(){ log "ERROR: $*"; exit 1; }

[[ -d "$RUN" ]] || die "RUN not found: $RUN"
[[ -d "$SINEDERELLA" ]] || die "SINEDERELLA not found: $SINEDERELLA"
[[ -s "$RUN/genome.clean.fa" ]] || die "missing genome.clean.fa"

exec >>"$LOG" 2>&1
log "=== oma orient + step1 rerun ==="
log "RUN=$RUN SINEDERELLA=$SINEDERELLA THREADS=$THREADS CHUNK_BP=$CHUNK_BP FLANK=$FLANK"

export PATH="/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin:$PATH"

cd "$SINEDERELLA"
git fetch origin main
git reset --hard origin/main
log "SINEderella at $(git rev-parse --short HEAD)"
chmod +x step7_boundary_refine.sh step8a_extract_alignments.sh run_publish_alignments.sh \
  publish/align_for_publish.sh docs/launch_oma_publish_post_rerun.sh 2>/dev/null || true

cd "$RUN"
[[ -s consensuses.clean.fa ]] || die "missing consensuses.clean.fa"

cp -f consensuses.clean.fa "consensuses.clean.fa.pre_orient.${STAMP}.bak"
log "Backed up consensuses.clean.fa"

for f in step1_search_extract.sh step2_asSINEment.sh step3_postprocess.sh step4_plots.sh \
         plot_subfamily.py canonicalize_consensus_bank.py consensus_bank_lib.py \
         rebuild_consensus_bank.py; do
  cp -f "$SINEDERELLA/$f" "./$f"
done
mkdir -p tools
cp -f "$SINEDERELLA/tools/orient_consensus_bank.py" "$SINEDERELLA/tools/orient_publish_aln.py" \
  "$SINEDERELLA/tools/audit_consensus_bank.py" tools/
chmod +x step1_search_extract.sh step2_asSINEment.sh step3_postprocess.sh step4_plots.sh \
  2>/dev/null || true

log "canonicalize (RC merge + tail orient)..."
python3 canonicalize_consensus_bank.py consensuses.clean.fa \
  -o consensuses.clean.fa.canon.tmp \
  --aliases "consensuses.aliases.orient_${STAMP}.tsv" \
  --min-id "${CANON_MIN_ID:-80}"
mv -f consensuses.clean.fa.canon.tmp consensuses.clean.fa

log "orient bank (idempotent pass)..."
python3 tools/orient_consensus_bank.py consensuses.clean.fa \
  --report "consensuses.orient_${STAMP}.tsv"

python3 <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, str(Path(".")))
from consensus_bank_lib import orient_by_simple_repeat_tail, read_fa
c = read_fa(Path("consensuses.clean.fa"))
for name in ("oma_SINE10", "oma_SINE27", "oma_grp080"):
    s = c.get(name, "")
    if not s:
        print(f"MISSING {name}")
        continue
    r = orient_by_simple_repeat_tail(s)
    print(f"{name}: orient={r.action} fwd={r.fwd_score:.1f} rev={r.rev_score:.1f} "
          f"start={r.seq[:24]} end={r.seq[-24:]}")
    if name == "oma_SINE10" and r.action != "kept":
        raise SystemExit(f"SINE10 not kept after orient ({r.action}) — abort")
PY

if [[ -d genome.clean_step1 ]]; then
  mv genome.clean_step1 "genome.clean_step1.pre_orient.${STAMP}"
  log "Moved genome.clean_step1 aside"
fi

log "step1 starting (sear with oriented bank)..."
set +e
bash ./step1_search_extract.sh ./consensuses.clean.fa ./genome.clean.fa "$CHUNK_BP" "$FLANK" \
  > >(tee "step1.orient_${STAMP}.stdout.log") \
  2> >(tee "step1.orient_${STAMP}.stderr.log" >&2)
step1_rc=$?
set -e
if (( step1_rc != 0 )); then
  tail -n 80 "step1.orient_${STAMP}.stderr.log" >&2 || true
  die "step1 failed rc=$step1_rc"
fi
if [[ ! -d genome.clean_step1 ]]; then
  mapfile -t _DIRS < <(find . -maxdepth 1 -mindepth 1 -type d -name '*_step1' -printf '%f\n' | sort)
  [[ ${#_DIRS[@]} -gt 0 ]] && mv -- "${_DIRS[0]}" genome.clean_step1
fi
[[ -s genome.clean_step1/extracted.fasta ]] || die "step1 missing extracted.fasta"

log "step2 starting..."
mkdir -p step2
(
  cd step2
  ../step2_asSINEment.sh ../consensuses.clean.fa ../genome.clean_step1/extracted.fasta \
    "$THREADS" "step2_output_rerun_orient_${STAMP}"
)
OUT="$(ls -dt step2/step2_output_rerun_orient_* 2>/dev/null | head -1)"
[[ -n "$OUT" && -d "$OUT" ]] || die "step2 output missing"
ln -sfn "$(basename "$OUT")" step2/step2_output_latest 2>/dev/null || true
log "step2 done: $OUT"

log "step3..."
./step3_postprocess.sh "$(pwd -P)" "$THREADS"

log "rebuild consensus bank..."
python3 rebuild_consensus_bank.py "$(pwd -P)"

log "step4..."
./step4_plots.sh "$(pwd -P)" "$THREADS"

log "publish alignments (orient bank + orient MSA)..."
export DISC
export PEEL_FLAGS="${PEEL_FLAGS:-$SINEDERELLA/oma_peel_border_flags.tsv}"
export PEEL_ALN_DIR="${PEEL_ALN_DIR:-/staging/tmp/scorpions/oma/sd/rebuild}"
export TRIM_DISPLAY_MODE="${TRIM_DISPLAY_MODE:-occupancy}"
bash "$SINEDERELLA/docs/launch_oma_publish_post_rerun.sh" "$RUN"

log "=== SINE10 check ==="
python3 <<'PY'
import re
from pathlib import Path
sys_path = Path(".")
aln = Path("alignments/oma_SINE10_top100.aln.fa")
if not aln.is_file():
    aln = Path("results/alignments/oma_SINE10_top100.aln.fa")
if not aln.is_file():
    print("MISSING SINE10 top100 aln"); raise SystemExit(0)
seqs, cur, buf = [], None, []
for line in aln.open():
    line = line.rstrip("\n\r")
    if line.startswith(">"):
        if cur is not None:
            seqs.append("".join(buf))
        cur = line[1:].split()[0]
        buf = []
    else:
        buf.append(line)
if cur is not None:
    seqs.append("".join(buf))
elem = re.sub(r"[^ACGT]", "", "".join(c for c in seqs[0] if c.isupper()))
print(f"SINE10 element len={len(elem)} start={elem[:30]} end={elem[-30:]}")
PY

log "=== COMPLETE ==="
log "Next: sync alignments/ to site repo and browser-verify SINE10 strand"
