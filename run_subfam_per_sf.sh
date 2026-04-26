#!/bin/bash
# run_subfam_per_sf.sh
# Run SubFam per subfamily for SAQ and CCR.
# Each subfamily: sample up to N_SAMPLE copies from assigned.fasta,
# skip if total < MIN_COPIES, then run SubFam + final mafft alignment.
# Output: {out_dir}/{sf}.al  (FASTA alignment: chunk-consensuses + sf-consensus)

set -euo pipefail

SAQ_RUN="/data/W/toki/Genomes/Mammalia/Eulipotyphla/saq/run_20260425_182219"
CCR_RUN="/data/W/toki/Genomes/Mammalia/Eulipotyphla/ccr/run_20260425_100809"
SUBFAM_BIN="/data/V/toki/bin/SubFam"
BINSIZE=50
N_SAMPLE=10000
MIN_COPIES=400
THREADS=$(nproc)

# ---- helper: run SubFam for one subfamily ----
run_one_sf() {
    local sf="$1"
    local fa="$2"       # input fasta (sampled copies)
    local cons="$3"     # single-sf consensus fasta
    local out_al="$4"   # final .al file path

    echo "[$(date '+%H:%M:%S')] START $sf ($(grep -c '^>' "$fa") copies)"

    local work
    work="$(dirname "$out_al")/_work_${sf}"
    rm -rf "$work"
    mkdir -p "$work"

    # SubFam needs to run from its work dir; input must be named input.fasta
    cp "$fa" "$work/input.fasta"
    (
        cd "$work"
        "$SUBFAM_BIN" input.fasta "$BINSIZE" > subfam.log 2>&1
    )

    if [[ ! -f "$work/input.clw" ]]; then
        echo "[$(date '+%H:%M:%S')] WARN  $sf: SubFam produced no input.clw (see $work/subfam.log)"
        return
    fi

    # Convert input.clw to FASTA if it looks like Clustal
    if head -1 "$work/input.clw" | grep -q "^CLUSTAL"; then
        awk '
          /^CLUSTAL/ || /^$/ { next }
          NF >= 2 {
            name=$1; seq=$2
            gsub(/[-.]/, "", seq)
            if (!(name in s)) { ord[++n]=name }
            s[name] = s[name] seq
          }
          END { for (i=1;i<=n;i++) { print ">" ord[i]; print s[ord[i]] } }
        ' "$work/input.clw" > "$work/input_reps.fasta"
    else
        cp "$work/input.clw" "$work/input_reps.fasta"
    fi

    # Final alignment: chunk-consensuses + sf consensus → .al
    cat "$work/input_reps.fasta" "$cons" > "$work/combined.fasta"
    mafft --thread "$THREADS" \
          --localpair \
          --maxiterate 1000 \
          --ep 0.123 \
          --nuc \
          --reorder \
          --preservecase \
          --quiet \
          "$work/combined.fasta" > "$out_al"

    if [[ -s "$out_al" ]]; then
        local n_seqs
        n_seqs=$(grep -c '^>' "$out_al")
        echo "[$(date '+%H:%M:%S')] OK    $sf → $(basename "$out_al") ($n_seqs seqs)"
    else
        echo "[$(date '+%H:%M:%S')] WARN  $sf: final mafft produced empty output"
    fi

    rm -rf "$work"
}

# ---- helper: process one species run ----
process_species() {
    local species="$1"
    local run_root="$2"

    local assigned="$run_root/step2/step2_output/assigned.fasta"
    local cons_fa="$run_root/consensuses.clean.fa"
    local out_dir="$run_root/results/subfam_alignments"

    mkdir -p "$out_dir"

    echo ""
    echo "========== $species =========="

    # Extract per-sf samples using python3
    /usr/local/bin/python3.12 - "$assigned" "$out_dir" "$N_SAMPLE" "$MIN_COPIES" <<'PYEOF'
import sys, random, os

assigned, out_dir, n_sample, min_copies = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

by_sf = {}
name = None
seq_lines = []
with open(assigned) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                parts = name.split('|')
                sf = parts[1] if len(parts) >= 2 else 'unknown'
                by_sf.setdefault(sf, []).append((name, ''.join(seq_lines)))
            name = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
if name:
    parts = name.split('|')
    sf = parts[1] if len(parts) >= 2 else 'unknown'
    by_sf.setdefault(sf, []).append((name, ''.join(seq_lines)))

rng = random.Random(42)
for sf, copies in sorted(by_sf.items()):
    n = len(copies)
    if n < min_copies:
        print('SKIP {}: {} copies < {}'.format(sf, n, min_copies), flush=True)
        continue
    take = min(n, n_sample)
    sampled = rng.sample(copies, take)
    print('{}: {} of {} copies'.format(sf, take, n), flush=True)
    fa = os.path.join(out_dir, sf + '.in.fasta')
    with open(fa, 'w') as fh:
        for hdr, seq in sampled:
            fh.write('>{}\n{}\n'.format(hdr, seq))
PYEOF

    # Extract individual sf consensus sequences
    /usr/local/bin/python3.12 - "$cons_fa" "$out_dir" <<'PYEOF'
import sys
from pathlib import Path

cons_fa, out_dir = sys.argv[1], sys.argv[2]
records = {}
name = None
seq_lines = []
with open(cons_fa) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                records[name] = ''.join(seq_lines)
            name = line[1:].split()[0]
            seq_lines = []
        else:
            seq_lines.append(line)
if name:
    records[name] = ''.join(seq_lines)

for sf_name, seq in records.items():
    out = Path(out_dir) / (sf_name + '.cons.fasta')
    out.write_text('>{}\n{}\n'.format(sf_name + '_CONSENSUS', seq))
    print('wrote consensus: {}'.format(out.name), flush=True)
PYEOF

    # Run SubFam for each extracted sf
    for fa in "$out_dir"/*.in.fasta; do
        sf=$(basename "${fa%.in.fasta}")
        cons="$out_dir/${sf}.cons.fasta"
        out_al="$out_dir/${sf}.al"

        if [[ ! -f "$cons" ]]; then
            echo "[$(date '+%H:%M:%S')] WARN  $sf: no consensus fasta found, skipping"
            continue
        fi

        run_one_sf "$sf" "$fa" "$cons" "$out_al"
    done

    echo "========== $species done =========="
}

# ---- main ----
process_species saq "$SAQ_RUN"
process_species ccr "$CCR_RUN"

echo ""
echo "All done."
