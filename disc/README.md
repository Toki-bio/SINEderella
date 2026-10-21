# disc/: the SINE-discriminator modules SINEderella runs

SINEderella's publish chain (`publish_run.sh`, `publish/align_for_publish.sh`, `publish/needs_border_loop.py`,
`step8a_extract_alignments.sh`) calls tools of SINE_discriminator: flank justification, the rebuilt consensus row, display
trimming, orientation and edge correction, the continuation rounds, the verdict and the page's "Subfamily alignments"
table. Until 2026-10-21 they ran from unversioned copies on each machine (therioserver `/home/toki/SINE_discriminator/site`,
DRAGEN `/staging/tmp/sinedisc`, `C:\work\SINE_discriminator\site`), so a run could not be reproduced from the repository.
This directory holds them, and the publish scripts use it unless `DISC` or `SINE_DISC` names another tree.

## What is here

The tools the publish chain calls, and every local module they import (closure by `import` statements; no dynamic imports,
no data files):

| called by SINEderella | imported by those |
|---|---|
| `boundary_justify.py`, `rebuild_consensus_row.py`, `trim_display_flanks.py`, `correct_published_aln.py`, `continuation.py`, `verdict.py`, `inject_disc_report.py`, `boundary.py` (needs_border_loop) | `edge_profile.py`, `extend_unique_flank.py`, `fix_alignments.py`, `flank_uniqueness.py`, `inject_oma_aln_section.py`, `measure_c.py`, `overall.py`, `parse_oma_copy_counts.py`, `profiles.py`, `test2_props.py`, `trim_flanks.py` |

## Where they came from

Copied 2026-10-21 from therioserver `/home/toki/SINE_discriminator/site`, the copy the rsi runs of 2026-10-16 to
2026-10-21 used. The local `C:\work\SINE_discriminator\site` copy differed in three files (`measure_c.py`, `profiles.py`,
`trim_flanks.py`); the server versions are the ones here. `inject_disc_report.py` includes the 2026-10-21 table change
(consensus length, assigned / SINE / grey / not SINE columns, `--bank`, `--families`).

Change these files here, in SINEderella, from now on; the old trees are no longer read by a default run.
