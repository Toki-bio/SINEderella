# Oma consensus-bank repair plan

Three workflow bugs block clean SINE10 discussion. Fixes land in SINEderella
`main` (canonicalize + rebuild + step4 `-3`). **Oma run on DRAGEN is not
repaired until the steps below are executed there.**

## Root cause (confirmed)

| Bug | Symptom | Fix |
|---|---|---|
| ± seeds kept as separate subfamilies | 6 names, 3 biology | `canonicalize_consensus_bank.py` before step1 |
| step4 pctid without `-3` | Divergence curves disagree for ± pairs | `step4_plots.sh` now uses `-g -3` (matches step2) |
| step4 uses AnnoSINE seeds | Ns in consensus; pctid vs wrong reference | `rebuild_consensus_bank.py` before step4 → `consensuses.rebuilt.fa` |

## Phase A — Audit current oma run (DRAGEN, read-only)

```bash
RUN=/staging/tmp/scorpions/oma/run_oma
cd /path/to/SINEderella   # sync repo first
python3 tools/audit_consensus_bank.py "$RUN" --min-id 90 \
  --focus oma_SINE27 oma_group34 oma_SINE10 oma_big76 oma_sub515 oma_grp080
```

Expected: three RC clusters among the six focus names; high N on seeds;
divergence medians diverge between pair members.

## Phase B — Repair pipeline (new code on DRAGEN)

1. Sync SINEderella repo to DRAGEN (`/staging/tmp/SINEderella` or run dir copy).
2. **Dry-run merge** on seeds only:
   ```bash
   python3 canonicalize_consensus_bank.py \
     "$RUN/consensuses.clean.fa" \
     -o /tmp/oma_cons.canon.fa \
     --aliases /tmp/oma_cons.aliases.tsv \
     --min-id 90 --dry-run
   ```
3. **Full re-run from step1** with canonical seeds (required — step1 search is per-consensus):
   ```bash
   cp /tmp/oma_cons.canon.fa /path/to/oma_seeds_canonical.fa
   # New run OR replace consensuses.clean.fa and re-run step1–4 with SKIP where safe
   export SKIP_CANONICALIZE=0   # default: merge RC
   export SKIP_REBUILD_CONS=0   # default: copy-majority before plots
   SINEderella genome.fa oma_seeds_canonical.fa
   ```
   **Do not** only patch `consensuses.clean.fa` on the existing run without
   re-running step1 — search hits and copy pools are per seed name.

## Phase C — After repair: SINE10 focus

Once one `oma_SINE10` subfamily exists (big76 merged or aliased):

- Re-run step4 plots → check single tight divergence curve
- Re-run publish border loop on `oma_SINE10` top100 only
- Then discriminator / manual review on one alignment geometry

## Env flags

| Variable | Default | Meaning |
|---|---|---|
| `SKIP_CANONICALIZE` | 0 | Set 1 to keep all AnnoSINE seed names |
| `CANON_MIN_ID` | 90 | RC/direct merge threshold (%) |
| `SKIP_REBUILD_CONS` | 0 | Set 1 to keep seed consensuses for step4 |

## Not in scope yet

- Post-hoc merge of an **existing** step2 assignment without step1 re-run
  (`merge_subfamilies_from_aliases.py` — TODO if full re-run is too expensive)
- Re-orienting gapped MSA consensus rows (your manual alignment stays valid;
  search/plot layer uses ungapped rebuilt consensuses)
