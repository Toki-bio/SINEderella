#!/usr/bin/env bash
# Oma consensus repair re-run: canonicalize seeds -> step2 -> step3 -> rebuild -> step4 -> step6
set -euo pipefail

RUN="${1:-/staging/tmp/scorpions/oma/run_oma}"
SINEDERELLA="${SINEDERELLA:-/staging/tmp/SINEderella}"
THREADS="${THREADS:-32}"
LOG="${RUN}/rerun_consensus_repair.log"

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG"; }
die(){ log "ERROR: $*"; exit 1; }

[[ -d "$RUN" ]] || die "RUN not found: $RUN"
[[ -d "$SINEDERELLA" ]] || die "SINEDERELLA not found: $SINEDERELLA"

exec >>"$LOG" 2>&1
log "=== oma consensus repair re-run ==="
log "RUN=$RUN SINEDERELLA=$SINEDERELLA THREADS=$THREADS"

export PATH="/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin:$PATH"

# ── sync repo ───────────────────────────────────────────────────────────────
cd "$SINEDERELLA"
git stash push -m "pre-rerun $(date +%F)" 2>/dev/null || true
git fetch origin main
git reset --hard origin/main
log "SINEderella at $(git rev-parse --short HEAD)"

# ── copy fresh scripts into run ─────────────────────────────────────────────
cd "$RUN"
for f in step2_asSINEment.sh step3_postprocess.sh step4_plots.sh plot_subfamily.py \
         canonicalize_consensus_bank.py consensus_bank_lib.py rebuild_consensus_bank.py; do
  cp -f "$SINEDERELLA/$f" "./$f"
done
chmod +x step2_asSINEment.sh step3_postprocess.sh step4_plots.sh 2>/dev/null || true

[[ -s consensuses.clean.fa ]] || die "missing consensuses.clean.fa"
[[ -s genome.clean_step1/extracted.fasta ]] || die "missing extracted.fasta"

# ── backup + restore seeds + canonicalize ───────────────────────────────────
BAK="$(ls -t consensuses.clean.fa.pre_repair.*.bak 2>/dev/null | head -1 || true)"
if [[ -n "$BAK" && -s "$BAK" ]]; then
  cp -f "$BAK" consensuses.clean.fa
  log "Restored seeds from $BAK"
else
  cp -f consensuses.clean.fa "consensuses.clean.fa.pre_repair.$(date +%Y%m%d).bak"
  log "Backed up consensuses.clean.fa"
fi

python3 canonicalize_consensus_bank.py consensuses.clean.fa \
  -o consensuses.clean.fa.canon.tmp \
  --aliases consensuses.aliases.tsv \
  --min-id "${CANON_MIN_ID:-80}"
mv -f consensuses.clean.fa.canon.tmp consensuses.clean.fa
log "Canonicalized bank: $(grep -c '^>' consensuses.clean.fa) consensuses"
cat consensuses.aliases.tsv || true

# Drop corrupted oma_sub515 (not RC of grp080; 31N seed, n=49)
if grep -q '^>oma_sub515' consensuses.clean.fa 2>/dev/null; then
  awk '/^>/{hdr=$0; sub(/^>/,"",hdr); sub(/[ \t].*/,"",hdr); skip=(hdr=="oma_sub515")} !skip' \
    consensuses.clean.fa > consensuses.clean.fa.nosub515
  mv consensuses.clean.fa.nosub515 consensuses.clean.fa
  echo -e "oma_sub515\toma_grp080\tdropped corrupt seed (31N, not RC)\tmanual" >> consensuses.aliases.tsv
  log "Dropped oma_sub515 from bank (-> grp080 alias noted)"
fi

# ── step2 ───────────────────────────────────────────────────────────────────
log "step2 starting..."
mkdir -p step2
(
  cd step2
  ../step2_asSINEment.sh ../consensuses.clean.fa ../genome.clean_step1/extracted.fasta \
    "$THREADS" step2_output_rerun_"$(date +%Y%m%d)"
)
OUT="$(ls -dt step2/step2_output_rerun_* 2>/dev/null | head -1)"
[[ -n "$OUT" && -d "$OUT" ]] || die "step2 output missing"
log "step2 done: $OUT"

# symlink as primary output for downstream tools expecting step2_output*
ln -sfn "$(basename "$OUT")" step2/step2_output_latest 2>/dev/null || true

# ── step3 ───────────────────────────────────────────────────────────────────
log "step3 starting..."
./step3_postprocess.sh "$(pwd -P)" "$THREADS"
log "step3 done"

# ── rebuild + step4 ─────────────────────────────────────────────────────────
log "rebuild consensus bank from copies..."
python3 rebuild_consensus_bank.py "$(pwd -P)"

log "step4 starting..."
./step4_plots.sh "$(pwd -P)" "$THREADS"
log "step4 done -> $OUT/plots/"

# ── step6 report (stock, no publish chain) ──────────────────────────────────
if [[ -f "$SINEDERELLA/step6_report.py" ]]; then
  cp -f "$SINEDERELLA/step6_report.py" .
  mkdir -p results
  log "step6 report..."
  python3 step6_report.py "$(pwd -P)" --out results/report.html \
    --no-sineplot --species-code oma \
    || log "WARNING: step6 failed"
fi

# ── summary ─────────────────────────────────────────────────────────────────
log "=== audit focus subfamilies ==="
python3 <<'PY'
import glob, re, statistics
from pathlib import Path
RC=str.maketrans('ACGT','TGCA')
RUN=Path(".")
s2=sorted(glob.glob("step2/step2_output_rerun_*"))[-1]
cons=RUN/"consensuses.clean.fa"
focus=["oma_SINE10","oma_SINE27","oma_grp080"]
for f in focus:
    t=Path(s2)/"plots"/f"{f}_pctid.tsv"
    if not t.exists():
        print(f, "no pctid"); continue
    div=sorted(100-float(l.split()[1]) for l in open(t) if l.strip())
    print(f, "median_div", div[len(div)//2], "n", len(div))
aliases=RUN/"consensuses.aliases.tsv"
if aliases.exists():
    print("aliases:", aliases.read_text().strip())
print("bank_n", sum(1 for _ in open(cons) if _.startswith(">")))
PY

log "=== COMPLETE ==="
