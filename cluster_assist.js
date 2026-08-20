// cluster_assist.js — First-pass subfamily clustering ASSIST for SubFam output
// Usage: node cluster_assist.js input.aln.fasta output_dir [--recurse] [--split-fraction=0.15]
//
// Runs MSA-viewer's SINEClusterer (vendored in subfam_cluster_lib.js) against an
// already-aligned, fixed-width FASTA (e.g. SubFam's input.clw after the degap+mafft
// realignment described in MANUAL.md §6.1.1 — do not feed it the raw discordant
// input.clw directly).
//
// VALIDATED CAVEAT (2026-08-21, against Tal/saq ground truth — 9 known manual
// subfamilies s1..s9 summing to 600 chunk-consensuses):
//   - Default settings only reliably resolve the coarsest structure. On a simple
//     2-subfamily case (hedgehog/eri) it recovered the correct 2-way split.
//     On saq's real 9-subfamily case it collapsed everything into 2 giant blobs
//     (293 + 307 seqs) plus one 6-seq cluster — no fine structure at all.
//   - Recursive re-clustering of oversized clusters (--recurse) does NOT rescue
//     this: re-running the algorithm on an already-homogeneous sub-pool with the
//     same default thresholds finds nothing further to split (tested: both the
//     293- and 307-seq saq sub-pools returned zero additional clusters even with
//     the relaxed-upper-bound retry).
//   - Conclusion: use this as a fast first pass to find the 1-2 most obvious
//     dominant groups and cut down what needs manual review — NOT as a
//     replacement for manual clustering, especially on datasets you expect to
//     contain more than 2-3 real subfamilies. Always inspect cluster sizes
//     against your prior expectation before trusting the output.
//
// Output (in output_dir):
//   summary.tsv               — cluster_id, size, n_perfect_features, n_occurrences
//   cluster_N.fasta           — aligned members of cluster N (gapped, from input)
//   cluster_N_consensus.fasta — single per-column majority consensus (ungapped)
//   unassigned.fasta          — sequences that didn't join any cluster ("noise")

const fs = require('fs');
const path = require('path');
const { SINEClusterer } = require('./subfam_cluster_lib.js');

const args = process.argv.slice(2);
const positional = args.filter(a => !a.startsWith('--'));
const flags = args.filter(a => a.startsWith('--'));
if (positional.length < 2) {
  console.error('Usage: node cluster_assist.js input.aln.fasta output_dir [--recurse] [--split-fraction=0.15]');
  process.exit(1);
}
const [inputPath, outDir] = positional;
const doRecurse = flags.includes('--recurse');
const splitFractionArg = flags.find(f => f.startsWith('--split-fraction='));
const SPLIT_FRACTION = splitFractionArg ? parseFloat(splitFractionArg.split('=')[1]) : 0.15;
const MIN_SIZE = 3;

if (!fs.existsSync(outDir)) fs.mkdirSync(outDir, { recursive: true });

function parseFasta(text) {
  return text.split(/(?=^>)/m).filter(Boolean).map(r => {
    const lines = r.split('\n');
    return { id: lines[0].slice(1).trim(), seq: lines.slice(1).join('').trim() };
  });
}

function consensus(members) {
  const len = members[0].seq.length;
  let out = '';
  for (let p = 0; p < len; p++) {
    const counts = {};
    for (const m of members) {
      const ch = m.seq[p];
      if (ch === '-' || ch === '.') continue;
      counts[ch] = (counts[ch] || 0) + 1;
    }
    let best = '-', bc = 0;
    for (const [ch, c] of Object.entries(counts)) if (c > bc) { bc = c; best = ch; }
    out += best;
  }
  return out.replace(/-/g, '');
}

const text = fs.readFileSync(inputPath, 'utf8');
const recs = parseFasta(text);
const lens = new Set(recs.map(r => r.seq.length));
if (lens.size !== 1) {
  console.error(`Input is not a fixed-width alignment (${lens.size} distinct lengths found). ` +
    'Degap and realign first (see MANUAL.md 6.1.1) before running this tool.');
  process.exit(1);
}
console.log(`loaded ${recs.length} seqs, aln length ${recs[0].seq.length}`);

const TOTAL = recs.length;
const finalClusters = [];
let noise = [];

function recurse(members, depth) {
  if (members.length < MIN_SIZE) { noise.push(...members); return; }
  const clusterer = new SINEClusterer(members);
  const result = clusterer.cluster({});
  for (const c of result.clusters) {
    const seqs = c.sequences.map(s => ({ id: s.id, seq: s.seq }));
    if (doRecurse && seqs.length > TOTAL * SPLIT_FRACTION && seqs.length < members.length) {
      const clustersBefore = finalClusters.length;
      const noiseBefore = noise.length;
      recurse(seqs, depth + 1);
      if (finalClusters.length === clustersBefore) {
        // recursion found nothing further — undo whatever it dumped into noise
        // and keep the original group intact as a single cluster instead
        noise.length = noiseBefore;
        finalClusters.push(seqs);
      }
    } else {
      finalClusters.push(seqs);
    }
  }
  noise.push(...result.unassigned.map(s => ({ id: s.id, seq: s.seq })));
}

recurse(recs, 0);
finalClusters.sort((a, b) => b.length - a.length);

const summaryLines = ['cluster_id\tsize'];
finalClusters.forEach((members, i) => {
  const id = `cluster_${i + 1}`;
  summaryLines.push(`${id}\t${members.length}`);
  fs.writeFileSync(path.join(outDir, `${id}.fasta`),
    members.map(m => `>${m.id}\n${m.seq}\n`).join(''));
  fs.writeFileSync(path.join(outDir, `${id}_consensus.fasta`),
    `>${id}_consensus_${members.length}seqs\n${consensus(members)}\n`);
});
fs.writeFileSync(path.join(outDir, 'summary.tsv'), summaryLines.join('\n') + '\n');
fs.writeFileSync(path.join(outDir, 'unassigned.fasta'),
  noise.map(m => `>${m.id}\n${m.seq}\n`).join(''));

console.log(`\n${finalClusters.length} clusters, sizes: ${finalClusters.map(c => c.length).join(', ')}`);
console.log(`${noise.length} sequences unassigned/noise`);
console.log(`Wrote output to ${outDir}/`);
console.log('\nASSIST TOOL ONLY — review cluster sizes against your expectation before trusting. ' +
  'See header comment for the validated caveat on multi-subfamily datasets.');
