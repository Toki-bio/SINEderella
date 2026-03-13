# Alignment Deployment Checklist

## CRITICAL: Use extract_alignments.sh, NOT step5_align_subfamilies.sh

**The correct script is `extract_alignments.sh`.**  
`step5_align_subfamilies.sh` produces `.aln.fa` files WITHOUT flanking sequences.  
`extract_alignments.sh` produces `.core.fa`, `.best50.fa`, `.subfam.fa` WITH 30bp left + 70bp right flanks.  
**Flanks are mandatory for SINE boundary inspection on the MSA viewer.**

---

## What each tier contains

| File | Content | Use |
|------|---------|-----|
| `{sf}.core.fa` | 50 random copies, flanked 30L+70R, MAFFT aligned | Conservation canvas + boundary inspection |
| `{sf}.best50.fa` | Top 50 by bitscore, flanked 30L+70R, MAFFT aligned | Best-quality boundary inspection |
| `{sf}.subfam.fa` | SubFam consensuses (only generated if ≥400 copies) | Subfamily structure, may be absent |
| `{sf}.evidence.txt` | Zero-hit evidence note | Absent subfamilies |

**SubFam link (Tier C) is only present for high-copy subfamilies. Do not add it as a required link.**

---

## How to run for a new species set

```bash
# On server:
nohup extract_alignments.sh \
  /path/to/consensus_bank.fa \
  /path/to/multi_run_dir/ \
  /path/to/multi_run_dir/alignments/ \
  > /tmp/extract_alignments.log 2>&1 &
echo "PID: $!"
```

Then download:
```bash
# Single SCP -r (one connection):
scp -r -i ~/.ssh/id_ed25519 \
  toki@85.89.102.78:/path/to/multi_run_dir/alignments/ \
  C:\work\SINEdb-deploy\alignments_tmp\
```

Then merge into `alignments/`, git add, git commit, git push.

---

## Verification (do this before claiming deployment is complete)

1. Open the live page (combined.html or squamata.html)
2. Click a cell for a **non-snake** species
3. Click "Core" → MSA viewer must show sequences WITH lowercase flanking DNA before/after the SINE
4. Click "Best50" → same check
5. Click "SubFam" → if 404, that is expected for low-copy subfamilies — remove link or suppress
6. Repeat for 3 different taxa before claiming it works

**Do NOT claim deployment is complete without completing these checks.**

---

## Key server paths

- Non-snake multi-run: `/data/W/toki/Genomes/Reptiles/Squamata/multi_20260218_201305/`
- Snake multi-run: `/data/W/toki/Genomes/Reptiles/Squamata/Snakes/multi_20260213_041600/`
- Consensus bank (non-snake): found per-species at `{species}/run_{date}/consensuses.clean.fa`
- Deploy repo: `C:\work\SINEdb-deploy\`
- Alignment target: `C:\work\SINEdb-deploy\alignments\{Species}\`
