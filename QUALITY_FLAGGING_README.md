# SINE Consensus Quality Flagging System

## Problem

Some SINE families contain mixed data (real sequences + false positive hits + degraded fragments). These cause:
- Slow convergence (>30 iterations, sometimes hitting max 50)
- High iteration-to-iteration variance (bootstrap samples wildly inconsistent)
- Low-confidence consensuses not suitable for downstream analysis

## Solution: Two-Part System

### 1. Post-Run Analysis (`analyze_convergence.sh`)

After a batch completes, analyze all logs to identify problematic families:

```bash
bash analyze_convergence.sh /data/V/toki/eryx_jayakari_consensus
```

**Output:**
- CSV with metrics per SINE (iterations, status, avg_change)
- Lists of slow convergers (>30 iters)
- Lists of non-converged families (hit max 50 iters)
- Flagged for manual review / data cleaning

**Interpretation:**
- `ITERS > 30` → likely mixed/FP data
- `NOT_CONVERGED` → definite problem
- `AVG_CHANGE > 0.01` → high variance, noisy family

### 2. Smart Convergence (`sine_consensus_smart.sh`)

Improved version of `sine_consensus.sh` that:
- Tracks variance over last N iterations
- Detects high-variance families automatically
- Implements early stopping (stops trying to converge bad data)
- Flags output with quality status

**Usage:**
```bash
bash sine_consensus_smart.sh [options] input.fasta [basename]

# Options (same as original):
#  -n SIZE    Subsample size (default 100)
#  -m ITERS   Max iterations (default 50)
#  -t PCT     Consensus threshold (default 50)
#  -c PCT     Min coverage (default 30)
#  -k         Keep trace + log
```

**Output includes:**
- `QUALITY: LOW_CONFIDENCE` in log if high variance detected
- Stops early instead of wasting time on unconvergent families
- Consensus still output (usable but flagged as suspicious)

**Quality thresholds (tunable):**
```bash
VARIANCE_WINDOW=3        # examine last 3 iterations
VARIANCE_LIMIT=0.15      # if avg change > 15%, flag as bad
```

Lower VARIANCE_LIMIT = stricter (more families flagged).
Raise it if false positive rate too high.

## Workflow for Data Refinement

### Step 1: Run current batch (normal)
```bash
stdbuf -oL -eL bash run_all.sh 2>&1 | tee run_all.log
```

### Step 2: Post-analysis
```bash
bash analyze_convergence.sh /data/V/toki/eryx_jayakari_consensus
```

### Step 3: Manual review of flagged families
- Download log + trace alignment
- Inspect with MSA viewer
- Identify FP clusters to remove

### Step 4: Re-filter input FASTA
```bash
# Remove identified FP sequences from original SINE FASTA
# Rerun consensus on cleaned version
```

### Step 5: Future batches with smart script
```bash
# Deploy sine_consensus_smart.sh to server
# Use in run_all.sh

sed -i 's|sine_consensus.sh|sine_consensus_smart.sh|' run_all.sh
stdbuf -oL -eL bash run_all.sh 2>&1 | tee run_all_smart.log
```

## Files

- [analyze_convergence.sh](analyze_convergence.sh) — Post-run analysis
- [sine_consensus_smart.sh](sine_consensus_smart.sh) — Smart convergence with early stopping
- Original [sine_consensus.sh](sine_consensus.sh) — Standard (no variance detection)

## Next Steps

1. **Current batch completes** → run `analyze_convergence.sh`
2. **Identify problem families** → examine logs/alignments
3. **Implement data cleaning** or accept low-confidence consensuses
4. **Test smart script** on next batch or on flagged families
5. **Tune variance thresholds** based on observed FP/TN rates
