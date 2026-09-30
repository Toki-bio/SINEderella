#!/usr/bin/env python3
"""Positional alignment-composition profile, ported into SINEderella itself
(step6_report.py embeds this) from SINE_discriminator's profiles.py +
measure_c.py, 2026-09-09 -- so every species' standard HTML report gets the
per-position diagram, not only the discriminator's oma-specific site.

Merged from two files with no behavior change other than dropping the
`fix_alignments` cross-repo import (consensus_index copied in verbatim,
3 lines, see below) and renaming the module. If the discriminator's originals
change, re-diff before re-porting -- this is a copy, not a shared dependency,
by deliberate choice (SINEderella must not runtime-depend on
SINE_discriminator).

Positional profiles: every statistic that varies along the element, plotted
against nucleotide position rather than collapsed to one number.

The x axis has three zones stitched together, so the Tier-2 cliff is visible as
a curve rather than inferred from a scalar:

    -70 .. -1     left flank, offset outward from the element edge
      0 .. L-1    consensus positions of the element itself
     +1 .. +70    right flank, offset outward from the element edge

Flank positions are indexed from each copy's OWN edge on ungapped sequence, so
no aligner is involved out there and the background is the true one.

Tracks
  pair_id   mean pairwise identity between copies - the one metric defined
            identically in all three zones, so the cliff is directly readable
  cover     fraction of copies present - shows truncation as a ramp
  cons_id   identity to the known consensus (element only)
  mosaic    per-window residual after a rank-1 fit of the copy x window
            identity matrix. Rank 1 means "every copy is just older or younger";
            what does not fit that is positional structure specific to a subset
            of copies, i.e. the spec's mosaicism, localised.
  at        A+T fraction - the A-rich insertion signature and the poly-A tail

Structural calls (measure(), from measure_c.py) -- A box, B box, self-similarity
to the tRNA-derived head, TSD, flank identity/cliff, mosaicism rank-1 excess --
are all designed to come back as None/0/nan when the feature genuinely is not
there. The report renderer must show "not detected", never hide the row, for
families that are not tRNA-derived (7SL/Alu/FLAM/unclassified): the same
uniform structural-feature panel applies to every subfamily, tRNA-derived or
not, so absence displays honestly instead of a track disappearing.
"""
import os
import sys
import json
import glob
import numpy as np

GAP = 4
CODE = {"a": 0, "c": 1, "g": 2, "t": 3, "A": 0, "C": 1, "G": 2, "T": 3}
W, STEP = 20, 4
FL = 70          # flank offsets profiled either side
WIN, PSTEP = 20, 4


def consensus_index(names):
    """Copied verbatim from SINE_discriminator's fix_alignments.py -- trivial
    (3 lines), not worth a cross-repo import for. Keep in sync by inspection,
    not by dependency."""
    for i, h in enumerate(names):
        if "CONSENSUS" in h.upper():
            return i
    return 0


def read_fa(p):
    """Raw string sequences (not the numeric-coded matrix read_aln() below
    builds) -- copied verbatim from SINE_discriminator's fix_alignments.py,
    needed by report_verdict.py and report_flank_uniqueness.py."""
    names, seqs, cur, buf = [], [], None, []
    for line in open(p, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = line[1:]
            names.append(cur)
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


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
    """Per-alignment structural summary stats (the chip row above the chart).
    Every field is designed to come back None/0.0/nan on genuine absence --
    the renderer must display that as "not detected", not omit the row."""
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

    inside = C[:, lo:hi + 1]
    elem_len = (inside != GAP).sum(axis=1).astype(float)
    v["elem_len_med"] = float(np.median(elem_len))
    v["elem_len_cv"] = float(np.std(elem_len) / max(1e-9, np.mean(elem_len)))
    v["elem_len_iqr"] = float(np.subtract(*np.percentile(elem_len, [75, 25])))
    v["frac_full"] = float(np.mean(elem_len > 0.9 * v["cons_bp"]))

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
    _ident = ident

    lefts, rights = [], []
    for i in range(n):
        row = C[i]
        l = row[:lo]
        l = l[l != GAP]
        r = row[hi + 1:]
        r = r[r != GAP]
        lefts.append(l[::-1])
        rights.append(r)
    v["flank_bp_med"] = float(np.median([len(x) + len(y) for x, y in zip(lefts, rights)]))
    v["flank_id_L"] = pair_identity_ungapped(lefts, rng)
    v["flank_id_R"] = pair_identity_ungapped(rights, rng)
    fl = np.nanmean([v["flank_id_L"], v["flank_id_R"]])
    v["flank_id"] = float(fl)
    v["cliff"] = float(v["cons_identity_med"] - fl)

    _q75 = float(np.percentile(_ident, 75))
    _contrast = _q75 - fl
    if _contrast < 0.15:
        v["frac_supported"] = 0.0
        v["support_threshold"] = None
    else:
        _thr = max(fl + 0.10, fl + 0.45 * _contrast)
        v["frac_supported"] = float(np.mean(_ident >= _thr))
        v["support_threshold"] = float(_thr)

    starts = np.arange(0, max(1, inside.shape[1] - W + 1), STEP)
    M = np.full((n, len(starts)), np.nan, dtype=np.float32)
    ca = np.concatenate([np.zeros((n, 1)), agree.astype(float).cumsum(axis=1)], axis=1)
    cp = np.concatenate([np.zeros((n, 1)), present.astype(float).cumsum(axis=1)], axis=1)
    for j, s in enumerate(starts):
        e = min(s + W, inside.shape[1])
        den = cp[:, e] - cp[:, s]
        with np.errstate(invalid="ignore"):
            M[:, j] = np.where(den >= 5, (ca[:, e] - ca[:, s]) / np.maximum(den, 1), np.nan)
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


def pair_identity_col(col):
    """Mean pairwise identity within one column of symbols (gaps excluded)."""
    b = col[col != GAP]
    n = len(b)
    if n < 4:
        return np.nan, n
    cnt = np.bincount(b, minlength=4)[:4].astype(float)
    same = float((cnt * (cnt - 1)).sum())
    return same / (n * (n - 1)), n


def profile(path):
    """Per-position track data for the interactive diagram, plus the
    A-box/B-box/self-similarity motif tracks. Returns None if no consensus
    row can be found (feeds the "not detected" panel, not an error)."""
    names, A = read_aln(path)
    if A is None:
        return None
    # SINEderella's own step8a output names the consensus row just
    # `>{species}_{subfamily}` (no "CONSENSUS" substring) and relies on
    # rebuild_consensus_row.py's reordering to put it at index 0 -- so use
    # consensus_index()'s 0-fallback here too (measure() already did; this
    # function's stricter "return None if no CONSENSUS-named row" check was
    # wrong for real SINEderella alignments, confirmed 2026-09-09 against
    # oma_SINE10_top100.aln.fa, whose header is plain `>oma_SINE10`).
    k = consensus_index(names)
    cons = A[k]
    nzc = np.where(cons != GAP)[0]
    if len(nzc) < 1:
        return None
    lo, hi = int(nzc[0]), int(nzc[-1])
    C = np.delete(A, k, axis=0)
    n = C.shape[0]
    L = len(nzc)

    xs, pair, cover, consid, at = [], [], [], [], []

    lefts = []
    for i in range(n):
        r = C[i, :lo]
        r = r[r != GAP]
        lefts.append(r[::-1])
    for off in range(FL, 0, -1):
        col = np.array([l[off - 1] if len(l) >= off else GAP for l in lefts], dtype=np.int8)
        p, m = pair_identity_col(col)
        xs.append(-off)
        pair.append(p)
        cover.append(m / n)
        consid.append(np.nan)
        b = col[col != GAP]
        at.append(float(np.mean(np.isin(b, [0, 3]))) if len(b) > 3 else np.nan)

    for p_i, c in enumerate(nzc):
        col = C[:, c]
        p, m = pair_identity_col(col)
        xs.append(p_i)
        pair.append(p)
        cover.append(m / n)
        present = col != GAP
        consid.append(float(np.mean(col[present] == cons[c])) if present.sum() >= 4 else np.nan)
        b = col[present]
        at.append(float(np.mean(np.isin(b, [0, 3]))) if len(b) > 3 else np.nan)

    rights = []
    for i in range(n):
        r = C[i, hi + 1:]
        rights.append(r[r != GAP])
    for off in range(1, FL + 1):
        col = np.array([r[off - 1] if len(r) >= off else GAP for r in rights], dtype=np.int8)
        p, m = pair_identity_col(col)
        xs.append(L - 1 + off)
        pair.append(p)
        cover.append(m / n)
        consid.append(np.nan)
        b = col[col != GAP]
        at.append(float(np.mean(np.isin(b, [0, 3]))) if len(b) > 3 else np.nan)

    El = C[:, nzc]
    pres = El != GAP
    agree = (El == cons[nzc][None, :]) & pres
    starts = np.arange(0, max(1, L - WIN + 1), PSTEP)
    Mw = np.full((n, len(starts)), np.nan)
    for j, s in enumerate(starts):
        e = min(s + WIN, L)
        den = pres[:, s:e].sum(axis=1)
        with np.errstate(invalid="ignore"):
            Mw[:, j] = np.where(den >= 6, agree[:, s:e].sum(axis=1) / np.maximum(den, 1), np.nan)
    keep = np.isfinite(Mw).mean(axis=1) > 0.8
    Mk = Mw[keep]
    colmed = np.nanmedian(Mk, axis=0)
    Mk = np.where(np.isfinite(Mk), Mk, colmed[None, :])
    mosaic = np.full(len(starts), np.nan)
    if Mk.shape[0] >= 8:
        U, S, Vt = np.linalg.svd(Mk - Mk.mean(), full_matrices=False)
        fit = np.outer(U[:, 0] * S[0], Vt[0])
        resid = (Mk - Mk.mean()) - fit
        denom = np.linalg.norm(Mk - Mk.mean(), axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mosaic = np.where(denom > 1e-9, np.linalg.norm(resid, axis=0) / denom, np.nan)
    mos_x = (starts + WIN // 2).tolist()

    import re as _re
    IU = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
          "W": "[AT]", "S": "[GC]", "K": "[GT]", "M": "[AC]", "N": "[ACGT]"}
    cseq = "".join("ACGT"[c] for c in cons[nzc])

    def raw_scan(seq, motif):
        pat = [IU[c] for c in motif]
        kk = len(motif)
        return [sum(1 for j in range(kk) if _re.match(pat[j], seq[i + j]))
                for i in range(len(seq) - kk + 1)]

    _rng = np.random.default_rng(0)

    def motif_track(motif, n_null=120):
        """z-score against shuffled-consensus null, else the whole element
        sits in a noisy 0.33-0.56 band and the real hit is not readable."""
        kk = len(motif)
        if L < kk + 5:
            return [None] * L
        obs = raw_scan(cseq, motif)
        arr = np.frombuffer(cseq.encode(), dtype="S1")
        null = []
        for _ in range(n_null):
            sh = arr.copy()
            _rng.shuffle(sh)
            null.extend(raw_scan(sh.tobytes().decode(), motif))
        mu, sd = float(np.mean(null)), float(np.std(null)) or 1.0
        z = [(v - mu) / sd for v in obs]
        return [round(float(x), 3) for x in z] + [None] * (L - len(z))

    abox_p = motif_track("TRGCNNARYGG")   # A box -- tRNA-Pol-III internal promoter
    bbox_p = motif_track("GWTCRANNC")     # B box -- ditto; both come back all-None
                                            # (never absent from the JSON) on a
                                            # non-tRNA-derived family: motif_track
                                            # still runs, the z-scores are just flat.

    head_lo, head_hi = 5, min(80, L)
    head = cseq[head_lo:head_hi]
    hl = len(head)
    raw_self = []
    for i in range(L):
        if i + hl > L:
            raw_self.append(None)
            continue
        seg = cseq[i:i + hl]
        raw_self.append(sum(1 for a, b in zip(seg, head) if a == b) / hl)
    _v = [x for x in raw_self if x is not None]
    _bg = float(np.median(_v)) if _v else 0.25
    _sd = float(np.std([x for x in _v if x < _bg + 0.15])) or 0.05
    self_p = []
    for i, x in enumerate(raw_self):
        if x is None or (head_lo - hl // 2 <= i <= head_hi):
            self_p.append(None)
        else:
            self_p.append(round((x - _bg) / _sd, 3))

    f = lambda a: [None if (x is None or not np.isfinite(x)) else round(float(x), 4) for x in a]
    return {"x": xs, "elem_len": L, "n": int(n),
            "pair_id": f(pair), "cover": f(cover), "cons_id": f(consid), "at": f(at),
            "mosaic_x": mos_x, "mosaic": f(mosaic),
            "motif_x": list(range(L)),
            "abox_p": abox_p, "bbox_p": bbox_p, "selfsim_p": self_p}


def profile_and_measure(path):
    """Convenience for step6_report.py: one file read, both halves."""
    p = profile(path)
    m = measure(path)
    return p, m
