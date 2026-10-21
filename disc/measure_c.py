#!/usr/bin/env python3
"""MEASURE(), rebuilt on the consensus anchor.

Differences from measure.py, all following from the two corrections:
  * the element is delimited by the KNOWN subfamily consensus, added to the
    alignment as a row - no threshold, no smoothing, no free parameter;
  * flanks are 70 bp, not 300, so MAFFT is not forced to align unrelated DNA;
  * flank identity is measured on UNGAPPED sequence walking outward from each
    copy's own element edge, never on aligned columns - aligning unrelated
    flanks is what inflated the background to 0.50 against a true 0.25.
"""
import os, sys, json, glob
import numpy as np
from fix_alignments import consensus_index

GAP = 4
CODE = {"a": 0, "c": 1, "g": 2, "t": 3, "A": 0, "C": 1, "G": 2, "T": 3}
W, STEP = 20, 4


def read_aln(path):
    names, seqs, cur = [], [], None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            names.append(line[1:])
            seqs.append([])
            cur = seqs[-1]
        elif cur is not None:
            cur.append(line)
    seqs = ["".join(s) for s in seqs]
    if not seqs:
        return names, None
    L = max(len(s) for s in seqs)
    A = np.full((len(seqs), L), GAP, dtype=np.int8)
    for i, s in enumerate(seqs):
        for j, ch in enumerate(s):
            A[i, j] = CODE.get(ch, GAP)
    return names, A


def pair_identity_ungapped(strings, rng, n_pairs=300, maxlen=70):
    """Mean pairwise identity of raw sequences compared position by position
    from a shared anchor. No aligner, so no similarity is manufactured."""
    ss = [s for s in strings if len(s) >= 25]
    if len(ss) < 4:
        return float("nan")
    vals = []
    for _ in range(n_pairs):
        a, b = rng.integers(0, len(ss), 2)
        if a == b:
            continue
        x, y = ss[a][:maxlen], ss[b][:maxlen]
        m = min(len(x), len(y))
        if m >= 25:
            vals.append(float(np.mean(x[:m] == y[:m])))
    return float(np.mean(vals)) if vals else float("nan")


def find_tsd(left, right, lo=6, hi=20):
    """Exact direct repeat flanking the element. The consensus anchor gives an
    exact edge in every copy, which is what makes this testable at all."""
    if len(left) < lo or len(right) < lo:
        return 0
    ls = "".join("ACGT"[c] for c in left[-(hi + 6):])
    rs = "".join("ACGT"[c] for c in right[:hi + 6])
    for k in range(min(hi, len(ls), len(rs)), lo - 1, -1):
        for i in range(len(ls) - k + 1):
            if ls[i:i + k] in rs:
                return k
    return 0


def measure(path, seed=0):
    rng = np.random.default_rng(seed)
    names, A = read_aln(path)
    v = {"set": os.path.basename(path).replace(".aln.fa", "")}
    if A is None or A.shape[0] < 9:
        v["error"] = "too_few_sequences"
        return v

    k = consensus_index(names)
    cons = A[k]
    nz = np.where(cons != GAP)[0]
    if len(nz) < 60:
        v["error"] = "consensus_absent"
        return v
    lo, hi = int(nz[0]), int(nz[-1])
    C = np.delete(A, k, axis=0)
    n = C.shape[0]
    v["n_copies"] = int(n)
    v["aln_len"] = int(A.shape[1])
    v["cons_bp"] = int(len(nz))
    v["cons_span"] = int(hi - lo + 1)
    v["cons_stretch"] = float(v["cons_span"] / v["cons_bp"])

    # ---- element length per copy, delimited by the consensus --------------
    inside = C[:, lo:hi + 1]
    elem_len = (inside != GAP).sum(axis=1).astype(float)
    v["elem_len_med"] = float(np.median(elem_len))
    v["elem_len_cv"] = float(np.std(elem_len) / max(1e-9, np.mean(elem_len)))
    v["elem_len_iqr"] = float(np.subtract(*np.percentile(elem_len, [75, 25])))
    v["frac_full"] = float(np.mean(elem_len > 0.9 * v["cons_bp"]))

    # ---- identity to the known consensus, per copy ------------------------
    cin = cons[lo:hi + 1]
    present = inside != GAP
    agree = (inside == cin[None, :]) & present
    with np.errstate(invalid="ignore"):
        ident = agree.sum(axis=1) / np.maximum(present.sum(axis=1), 1)
    ident = ident[present.sum(axis=1) >= 40]
    if len(ident) < 6:
        v["error"] = "no_supported_copies"
        return v
    v["cons_identity_med"] = float(np.median(ident))
    v["cons_identity_iqr"] = float(np.subtract(*np.percentile(ident, [75, 25])))
    v["_ident"] = ident            # consumed below, once the background is known

    # ---- flank: ungapped, walking outward from each copy's own edge -------
    lefts, rights = [], []
    for i in range(n):
        row = C[i]
        l = row[:lo]
        l = l[l != GAP]
        r = row[hi + 1:]
        r = r[r != GAP]
        lefts.append(l[::-1])          # reversed: index 0 = adjacent to element
        rights.append(r)
    v["flank_bp_med"] = float(np.median([len(x) + len(y) for x, y in zip(lefts, rights)]))
    v["flank_id_L"] = pair_identity_ungapped(lefts, rng)
    v["flank_id_R"] = pair_identity_ungapped(rights, rng)
    fl = np.nanmean([v["flank_id_L"], v["flank_id_R"]])
    v["flank_id"] = float(fl)
    v["cliff"] = float(v["cons_identity_med"] - fl)

    # Support is judged against THIS family's own contrast, not an absolute cut.
    # A constant 0.55 marks an old but clean family as contaminated, because
    # genuine old copies fall below any fixed threshold.
    _id = v.pop("_ident")
    _q75 = float(np.percentile(_id, 75))
    _contrast = _q75 - fl
    if _contrast < 0.15:
        v["frac_supported"] = 0.0          # no element: nothing supports it
        v["support_threshold"] = None
    else:
        _thr = max(fl + 0.10, fl + 0.45 * _contrast)
        v["frac_supported"] = float(np.mean(_id >= _thr))
        v["support_threshold"] = float(_thr)

    # ---- per-copy identity landscape and rank ----------------------------
    starts = np.arange(0, max(1, inside.shape[1] - W + 1), STEP)
    M = np.full((n, len(starts)), np.nan, dtype=np.float32)
    ca = np.concatenate([np.zeros((n, 1)), agree.astype(float).cumsum(axis=1)], axis=1)
    cp = np.concatenate([np.zeros((n, 1)), present.astype(float).cumsum(axis=1)], axis=1)
    for j, s in enumerate(starts):
        e = min(s + W, inside.shape[1])
        den = cp[:, e] - cp[:, s]
        with np.errstate(invalid="ignore"):
            M[:, j] = np.where(den >= 5, (ca[:, e] - ca[:, s]) / np.maximum(den, 1), np.nan)
    # Dropping every row with a NaN window leaves almost nothing once copies
    # carry indels, which silently disabled the whole rank statistic. Impute
    # the column median instead and record how much had to be filled.
    colmed = np.nanmedian(M, axis=0)
    fillfrac = float(np.mean(~np.isfinite(M)))
    Mc = np.where(np.isfinite(M), M, colmed[None, :])
    Mc = Mc[np.isfinite(Mc).all(axis=1)]
    v["rank_fill_frac"] = fillfrac
    if Mc.shape[0] >= 8 and Mc.shape[1] >= 4:
        S = np.linalg.svd(Mc, compute_uv=False)
        tot = float((S ** 2).sum())
        v["rank1_frac"] = float(S[0] ** 2 / tot)
        v["rank2_frac"] = float(S[1] ** 2 / tot)
        nulls = []
        for _ in range(40):
            P = Mc.copy()
            for j in range(P.shape[1]):
                rng.shuffle(P[:, j])
            Sn = np.linalg.svd(P, compute_uv=False)
            nulls.append(Sn[0] ** 2 / float((Sn ** 2).sum()))
        v["rank1_null"] = float(np.mean(nulls))
        v["rank1_excess"] = float(v["rank1_frac"] - v["rank1_null"])
        v["n_rank_rows"] = int(Mc.shape[0])

    # ---- per-copy edge deviation from the consensus edge ------------------
    resL, resR = [], []
    for i in range(n):
        p = np.where(inside[i] != GAP)[0]
        if len(p) < 40:
            continue
        resL.append(int(p[0]))
        resR.append(int(inside.shape[1] - 1 - p[-1]))
    if len(resL) >= 6:
        resL, resR = np.array(resL, float), np.array(resR, float)
        v["resL_med"] = float(np.median(resL))
        v["resR_med"] = float(np.median(resR))
        v["resL_iqr"] = float(np.subtract(*np.percentile(resL, [75, 25])))
        v["resR_iqr"] = float(np.subtract(*np.percentile(resR, [75, 25])))
        v["res_asymmetry"] = float(np.log((v["resL_iqr"] + 1) / (v["resR_iqr"] + 1)))

    # ---- one-sided bonuses, now that edges are exact ----------------------
    tsd = [find_tsd(lefts[i][::-1], rights[i]) for i in range(min(n, 120))]
    v["tsd_frac"] = float(np.mean([t > 0 for t in tsd]))
    v["tsd_len_med"] = float(np.median([t for t in tsd if t > 0])) if any(tsd) else 0.0
    _u = [l[:20] for l in lefts if len(l) >= 20]
    up = np.concatenate(_u) if _u else np.array([])
    v["arich_score"] = float(np.mean(np.isin(up, [0, 3]))) if len(up) > 100 else float("nan")
    _d = [r[:25] for r in rights if len(r) >= 25]
    dn = np.concatenate(_d) if _d else np.array([])
    v["polyA_score"] = float(np.mean(dn == 0)) if len(dn) > 100 else float("nan")
    return v


def _one(f):
    try:
        return measure(f)
    except Exception as exc:
        return {"set": os.path.basename(f).replace(".aln.fa", ""), "error": "exception: %s" % exc}


def main():
    import multiprocessing as mp
    d = sys.argv[1] if len(sys.argv) > 1 else "aln_c"
    out = sys.argv[2] if len(sys.argv) > 2 else "features_c.jsonl"
    files = sorted(glob.glob(os.path.join(d, "*.aln.fa")))
    print("measuring %d" % len(files), flush=True)
    with mp.Pool(max(1, mp.cpu_count() - 2)) as pool, open(out, "w") as fh:
        for i, v in enumerate(pool.imap_unordered(_one, files, chunksize=4)):
            fh.write(json.dumps(v) + "\n")
            if i % 100 == 0:
                print("  %d" % i, flush=True)
    print("wrote", out)


if __name__ == "__main__":
    main()
