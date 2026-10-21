#!/usr/bin/env python3
"""Per-copy edge profiling and classification -- the general, reusable rule
this was built to be, not a one-off patch for oma_SINE16.

Motivating case (Sergei, 2026-09-10): a family-level flank-uniqueness call
(one number, or one island_scan over the whole flank window) collapses three
genuinely different situations into the same verdict:

  SHARP        nearly all copies transition to background at ~the same,
               short offset -> one clean shared boundary, extend once and done.
  GRADED_DECAY transition offsets spread smoothly with no real gaps -> a
               LINE-like decaying tail, not a sharp SINE edge; there is no
               single "the" boundary to extend to.
  MULTI_GROUP  transition offsets cluster into 2+ separated groups -> these
               are probably structurally different populations (a real
               subfamily split, or some copies sitting in a different,
               longer shared context) and should be peeled apart, not
               averaged together.
  UNRESOLVED   a copy never reaches background within the explored length
               (max_off) -- must be reported as exactly that, "still shared
               at Nbp, not further explored", never silently treated as
               either independent or as a settled shared boundary.

All four are real, valid answers. The point of this module is to always
report which one actually holds, with the evidence (group sizes, edge
ranges, unresolved count), rather than forcing every case into a single
"independent y/n" call.

Cross-check: the SAME positional pattern should show up as elevated
mosaicism (report_profile.py's rank-1-residual mosaic track) in the
disputed region for MULTI_GROUP and GRADED_DECAY cases with real subset
structure -- a family-level match-rate split necessarily produces a
rank-1-inconsistent (mosaic) signal there. Verified on oma_SINE16: the
mosaic track dips to its cleanest point right where a shared short motif
sits (~310bp) and rises on both sides of it, consistent with the per-copy
edge groups found independently here.
"""
import sys
import numpy as np

GAP = 4
CODE = {"a": 0, "c": 1, "g": 2, "t": 3, "A": 0, "C": 1, "G": 2, "T": 3}
WIN = 15          # smoothing window for the per-copy match-rate curve
BG_MARGIN = 0.12  # match rate must fall within this of background to call "edge"

# cluster_edges thresholds
MIN_GAP_RATIO = 2.5
MIN_GAP_ABS = 15
MIN_GROUP = 3

# classification thresholds
SHARP_MAX_MEDIAN = 20      # a single group with median edge <= this reads SHARP
SHARP_MIN_FRAC = 0.80      # ... and covering at least this fraction of resolved copies
UNRESOLVED_NOTABLE_FRAC = 0.05  # unresolved copies above this fraction is reported


def read_aln(path):
    names, seqs, cur = [], [], None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            names.append(line[1:]); seqs.append([]); cur = seqs[-1]
        elif cur is not None:
            cur.append(line)
    seqs = ["".join(s) for s in seqs]
    L = max(len(s) for s in seqs)
    A = np.full((len(seqs), L), GAP, dtype=np.int8)
    for i, s in enumerate(seqs):
        for j, ch in enumerate(s):
            A[i, j] = CODE.get(ch, GAP)
    return names, A


def consensus_index(names):
    for i, h in enumerate(names):
        if "CONSENSUS" in h.upper():
            return i
    return 0


def per_copy_edges(path, side, max_off=1000):
    """side: 'L' or 'R'. Returns (results, bg_rate, names) where results is
    a list of (copy_index_in_others, edge_offset_or_None, curve). edge is
    None when the copy never reaches background within max_off -- reported
    as unresolved, never defaulted to 0 or to max_off."""
    names, A = read_aln(path)
    ci = consensus_index(names)
    cons = A[ci]
    nz = np.where(cons != GAP)[0]
    lo, hi = int(nz[0]), int(nz[-1])
    C = np.delete(A, ci, axis=0)
    block = C[:, :lo] if side == "L" else C[:, hi + 1:]
    seqs = []
    for r in block:
        x = r[r != GAP]
        seqs.append(x[::-1] if side == "L" else x)

    L_off = min(max_off, max((len(s) for s in seqs), default=0))
    maj = np.full(L_off, GAP, dtype=np.int8)
    for off in range(L_off):
        col = np.array([s[off] if len(s) > off else GAP for s in seqs], dtype=np.int8)
        b = col[col != GAP]
        if len(b) >= 8:
            maj[off] = int(np.bincount(b, minlength=4)[:4].argmax())

    bg_matches, bg_n = 0, 0
    for off in range(L_off):
        if maj[off] == GAP:
            continue
        col = np.array([s[off] if len(s) > off else GAP for s in seqs], dtype=np.int8)
        b = col[col != GAP]
        bg_matches += int((b == maj[off]).sum())
        bg_n += len(b)
    bg_rate = bg_matches / float(bg_n) if bg_n else 0.25

    results = []
    for i, s in enumerate(seqs):
        n = min(len(s), L_off)
        if n < WIN * 2:
            results.append((i, None, None))
            continue
        match = np.array([1.0 if (s[j] != GAP and maj[j] != GAP and s[j] == maj[j])
                           else (np.nan if (s[j] == GAP or maj[j] == GAP) else 0.0)
                           for j in range(n)])
        curve = np.full(n, np.nan)
        for j in range(n):
            a, b = max(0, j - WIN // 2), min(n, j + WIN // 2 + 1)
            w = match[a:b]
            w = w[np.isfinite(w)]
            if len(w) >= 3:
                curve[j] = float(np.mean(w))
        edge = None
        run = 0
        for j in range(n):
            if np.isfinite(curve[j]) and curve[j] <= bg_rate + BG_MARGIN:
                run += 1
                if run >= WIN:
                    edge = j - WIN + 1
                    break
            else:
                run = 0
        results.append((i, edge, curve))
    return results, bg_rate, names


def cluster_edges(edges, min_gap_ratio=MIN_GAP_RATIO, min_gap_abs=MIN_GAP_ABS, min_group=MIN_GROUP):
    """Split sorted edge offsets into ALL plausible groups, not just the
    single largest gap. A gap counts as a real split when it is both
    absolutely large (>= min_gap_abs, matching the smoothing window's own
    resolution) AND large relative to the median gap elsewhere, so one
    sparse stretch doesn't get chopped into spurious one-copy groups."""
    s = sorted(edges)
    n = len(s)
    if n < 2 * min_group:
        return [s] if s else []
    gaps = [s[i + 1] - s[i] for i in range(n - 1)]
    med_gap = sorted(gaps)[len(gaps) // 2] if gaps else 0
    cuts = [i + 1 for i, g in enumerate(gaps)
            if g >= min_gap_abs and g >= min_gap_ratio * max(1, med_gap)]
    groups, start = [], 0
    for c in cuts:
        groups.append(s[start:c]); start = c
    groups.append(s[start:])
    merged = []
    for g in groups:
        if merged and len(g) < min_group:
            merged[-1].extend(g)
        else:
            merged.append(list(g))
    return merged


def classify(path, side, max_off=1000):
    """The operational call: SHARP | GRADED_DECAY | MULTI_GROUP, plus
    unresolved-fraction reporting. Returns a dict, not a single label --
    the groups/evidence must travel with the call so a verdict can act on
    them (e.g. which copies to peel into a subfamily)."""
    results, bg_rate, names = per_copy_edges(path, side, max_off)
    n_total = len(results)
    edges = [e for _, e, _ in results if e is not None]
    unresolved_idx = [i for i, e, _ in results if e is None]
    groups = cluster_edges(edges)
    unresolved_frac = len(unresolved_idx) / float(n_total) if n_total else 0.0

    if not groups:
        pattern = "UNRESOLVED" if unresolved_idx else "NO_DATA"
    elif len(groups) == 1:
        g = groups[0]
        med = sorted(g)[len(g) // 2]
        frac = len(g) / float(len(edges)) if edges else 0.0
        pattern = "SHARP" if (med <= SHARP_MAX_MEDIAN and frac >= SHARP_MIN_FRAC) else "GRADED_DECAY"
    else:
        pattern = "MULTI_GROUP"

    group_summaries = []
    for gi, g in enumerate(groups):
        gset = set(g)
        members = [names[i + 1] for i, e, _ in results if e is not None and e in gset]
        group_summaries.append({
            "group": gi, "n": len(g), "edge_min": min(g), "edge_max": max(g),
            "edge_median": sorted(g)[len(g) // 2], "members": members,
        })

    return {
        "side": side, "bg_rate": round(bg_rate, 3), "n_total": n_total,
        "n_resolved": len(edges), "n_unresolved": len(unresolved_idx),
        "unresolved_frac": round(unresolved_frac, 3),
        "unresolved_members": [names[i + 1] for i in unresolved_idx],
        "pattern": pattern, "groups": group_summaries,
    }


def recommend(classifications):
    """Turn one or more classify() results (typically one per side) into the
    FINAL-vs-NEEDS_ADDITIONAL_WORK status + concrete next steps this was
    asked to produce -- so a verdict never silently reports a number without
    saying whether the analysis is actually done."""
    steps = []
    for c in classifications:
        side = c["side"]
        if c["pattern"] == "GRADED_DECAY":
            steps.append("%s side: graded decay (median edge %s over %d copies) -- "
                         "likely a LINE-like tail, not a sharp SINE boundary; "
                         "re-examine whether this family is SINE at all on this side"
                         % (side, c["groups"][0]["edge_median"] if c["groups"] else "?",
                            c["n_resolved"]))
        elif c["pattern"] == "MULTI_GROUP":
            sizes = ", ".join("%d@%d-%dbp" % (g["n"], g["edge_min"], g["edge_max"])
                               for g in c["groups"])
            steps.append("%s side: %d distinct edge groups (%s) -- "
                         "peel into separate subgroups before re-scoring, "
                         "do not average them into one boundary"
                         % (side, len(c["groups"]), sizes))
        if c["unresolved_frac"] > UNRESOLVED_NOTABLE_FRAC:
            steps.append("%s side: %d copies (%.0f%%) still shared at the max explored "
                         "length -- extend further or state explicitly as unresolved, "
                         "never treat as independent" % (side, c["n_unresolved"],
                                                          100 * c["unresolved_frac"]))
    status = "FINAL" if not steps else "NEEDS_ADDITIONAL_WORK"
    return status, steps


def main():
    path, side = sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "R"
    c = classify(path, side)
    print("side=%s pattern=%s bg_rate=%.3f  n_total=%d n_unresolved=%d (%.0f%%)"
          % (c["side"], c["pattern"], c["bg_rate"], c["n_total"],
             c["n_unresolved"], 100 * c["unresolved_frac"]))
    for g in c["groups"]:
        print("  group %d: n=%d edge=[%d,%d] median=%d  e.g. %s"
              % (g["group"], g["n"], g["edge_min"], g["edge_max"], g["edge_median"],
                 g["members"][0] if g["members"] else "-"))
    status, steps = recommend([c])
    print("status:", status)
    for s in steps:
        print(" -", s)


if __name__ == "__main__":
    main()
