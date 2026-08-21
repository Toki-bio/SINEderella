# Research directions: two competing models of repeat divergence

SINEderella's current pipeline (consensus voting, boundary refinement,
pairwise-identity flank testing) assumes a **gradual, substitution-driven**
model of repeat evolution: every copy descends from one ancestral consensus
by accumulating point mutations (and occasional indels), so relatedness is
a matter of degree, measured as pairwise or consensus identity. This is the
model behind every distance-based step in the pipeline today — subfamily
assignment via ssearch36 voting, the boundary-refinement fraction-of-pairs
test, SubFam's chunk-consensus clustering.

There is a second, structurally different model worth treating as a first-
class research pathway rather than noise to be filtered out of the first:

## The modular / "repeat pangenome" (panconsensus) model

Some interspersed repeats — candidate case: `eri`'s larger SINE-adjacent
elements — may not be single evolving sequences at all, but **composites
assembled from a small library of reusable, independently-mobile modules**
(shared 3' tails, internal domains, TSD-adjacent motifs, or fragments
recruited from unrelated repeat families). Under this model, "how related
are two copies" is not a single identity number — it is a **presence/
absence and arrangement pattern over a shared module library**, the same
way a bacterial pangenome describes a species not as one reference genome
but as core genes (present in ~all strains) plus a shared pool of
accessory/variable genes that individual strains draw from.

The name fits directly: a "panconsensus" for a repeat family would be a
**module library plus a per-copy presence/absence/order matrix**, not a
single linear consensus sequence. Divergence between two copies is then
partly substitution-driven (within a shared module) and partly
combinatorial (which modules are present, and in what order/orientation).

### Motivating evidence already in hand

The `eri` e2-3 boundary-extension work (2026-08-21, see
[Tal-repo eri/LOG.md](https://github.com/Toki-bio/Tal/blob/main/eri/LOG.md))
found a population that is **bimodal**, not gradually diverging, at the
downstream flank: roughly half of sampled copies still share an exact
motif (`gcaggcaccgag...cccagcaataa`) even 100-120bp past the called SINE
end, while the other half look like independent background DNA at the same
offset. The gradual/pairwise-identity framing treated this as "boundary
under-calling" (a single true boundary, just measured wrong) and required
a fix (mean → fraction-of-pairs-above-threshold) to even detect the
bimodality rather than average it away.

The modular framing offers a different, competing explanation for the same
data: this is not one boundary being under-called for everyone — it may be
a **module that is present in some copies and absent in others**, e.g. a
recruited tail fragment that hitchhikes with some insertion events but not
others (independent of divergence time), or a nested/adjacent repeat
family whose own instances are getting pulled into the same window. These
predict different downstream signatures: gradual under-calling should
correlate boundary distance with subfamily age/divergence; a
present/absent module should instead produce a clean two-cluster split
with no such gradient, plus (testably) the "with-module" cluster's tail
sequence should itself match a well-defined, low-diversity motif rather
than a diffuse, aging-consistent one — which the raw sequences sampled so
far are at least consistent with, but this has not been tested directly.

## Why these are genuinely different research pathways, not two views of the same thing

| | Gradual / SNP-indel model | Modular / panconsensus model |
|---|---|---|
| Unit of comparison | Pairwise or consensus identity (one number) | Module presence/absence/order (a matrix) |
| "Subfamily" means | A cluster around one consensus, tightening with more data | A profile over a shared module library — two "subfamilies" can share some modules and not others |
| Divergence signature | Smooth, roughly clock-like decay | Discrete jumps at module gain/loss events; identity within a module can still be clock-like |
| Boundary "fuzziness" | Measurement/calling error around one true edge | Possibly a real biological on/off state, not error |
| Detection method | Alignment + identity threshold (what SINEderella already does) | Segment-and-cluster: break each copy into candidate blocks, cluster blocks across the whole population, ask which copies share which blocks |
| Existing infrastructure fit | Everything current (step2 assignment, step7 boundary refinement, SubFam) | None yet — needs a new comparison primitive |

## Proposed approach: additive, not a rewrite

Do not replace the gradual model — it is cheap, well-understood, and
already load-bearing for every existing SINEderella step. Instead, add a
**new, independent analysis step** that consumes SINEderella's normal
output (assigned loci + their genomic flanks) and asks a structurally
different question:

1. **Segment**: for a candidate subfamily's flanking/boundary region
   (exactly the zone the boundary-refinement step already extracts),
   break each copy's extended window into fixed- or breakpoint-detected
   sub-windows rather than treating it as one unit.
2. **Cluster blocks, not copies**: cluster those sub-windows across the
   whole sampled population (e.g. by k-mer profile or local alignment),
   independent of which parent copy they came from, to find recurring
   discrete motifs/blocks — this is the actual test for "is there a
   reusable module here" as opposed to "is there a shared ancestor."
3. **Build the presence matrix**: for each original copy, record which
   discovered blocks it does/doesn't contain, in what order. This matrix
   *is* the panconsensus — analogous to a pangenome's gene
   presence/absence matrix.
4. **Compare the two models' fit**: check whether block presence/absence
   correlates with anything the gradual model already measures
   (subfamily age proxies, tandem-repeat overlap, genomic location) —
   a real modular signal should show discrete clustering that a smooth
   divergence model cannot produce by construction.

This is deliberately scoped as **a new, standalone exploratory step**
(matching the project's existing modular-step philosophy — see
`step7_boundary_refine.sh`, designed to be runnable independently of the
main orchestrator) rather than a modification to any existing step. It
should be tried first against the `eri` e2-3 bimodal case, since that is
the one dataset already known to show the signature this model predicts.

## Status

Idea only, as of 2026-08-21 — not designed in implementation detail, not
scoped into a task spec, no code written. Recorded here so it is not lost
and so the two-pathway framing (gradual vs. modular) is explicit before
either gets more infrastructure built on top of it.
