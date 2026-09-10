#!/usr/bin/env python3
# Ported into SINEderella from SINE_discriminator's verdict.py (2026-09-09),
# import lines repointed to report_profile.py (this repo's port of
# measure_c.py / fix_alignments.consensus_index) -- logic below is otherwise
# unmodified. See that file's own docstring (right below) for the design
# rationale, and read this comment before touching anything for the caveat.
#
# Side-car JSON files this module optionally reads (flankdecay.json,
# flank_islands.json, microsat.json, region_part.json) come from separate
# SINE_discriminator scripts (flankdecay.py, islands_corpus.py, microsat.py,
# region_partition.py) that run on a 400bp-flank corpus vanilla SINEderella
# does not build. Their absence is handled already (try/except -> {}): the
# verdict still computes, just without the flags/groups those side-cars would
# add (FLANK_ISLANDS, MICROSATELLITE_*, CONSENSUS_UNDEREXTENDED, and g_uniq
# stays neutral/unmeasured). Same for the optional test2_props import used
# for HETEROGENEOUS_SELECTION. Porting those side-cars natively is future
# work, not done here.
"""A weighted verdict with named sub-cases, instead of a binary call.

Two things the earlier builds got wrong, both raised by Sergei:

1. A single accept/reject throws away most of what was measured. The verdict is
   now a weighted score over the evidence groups, with every component visible so
   a disagreement can be traced to the variable that caused it.

2. "Negative" was one bucket. It is not one thing. A set can fail because there
   is no element at all, or because there IS an element and something else is
   also true - copies drawn from several structurally different subfamilies, or
   copies sitting inside other repeats. Those are different pieces of future
   work, so they get different names.

**A core of real-looking copies is always reported, however contaminated the set
is.** Twenty or thirty genuine copies inside a heavily polluted candidate are
worth knowing about - they are where the next clean family comes from.

The weights are stated, not fitted. Fitting now would fit the synthetic
negatives, which is the thing this project has been careful not to do.
"""
import json, glob, os, sys
import numpy as np
import report_profile as M
from report_profile import consensus_index

W = {                       # evidence group -> weight
    "element": 0.45,        # is there an element supported by the copies
    "homogeneity": 0.25,    # is it ONE element, uniform across copies
    "uniqueness": 0.20,     # are the copies in unique genomic locations
    "insertion": 0.10,      # one-sided bonuses: TSDs and the insertion site
}
SUPPORT = 0.55              # legacy constant, kept for reference only


def turbulence(el, pres):
    """How unsettled is the element along its length?

    Sergei on hedgehog e2-3: "can be SINE, but very turbulent at many places,
    should be regarded as SINE but with caution and need manual reinspection."
    Measured as the fraction of element columns whose identity falls more than
    0.25 below the element's own median. A clean family sits at 0.06-0.08;
    e2-3 is 0.18.

    This does not change the SINE verdict. It marks a set as one a person should
    look at, which is a different thing from a set being wrong.
    """
    idc = []
    for c in range(el.shape[1]):
        b = el[:, c][pres[:, c]]
        if len(b) >= 8:
            cnt = np.bincount(b, minlength=4)[:4].astype(float)
            idc.append((cnt * (cnt - 1)).sum() / (len(b) * (len(b) - 1)))
        else:
            idc.append(np.nan)
    idc = np.array(idc, float)
    if not np.isfinite(idc).sum() > 40:
        return None
    med = float(np.nanmedian(idc))
    return float(np.nanmean(idc < med - 0.25))


def contamination_split(ident, bg, min_gap=0.10, min_frac=0.05):
    """Contamination is a SEPARATE MODE, not a low tail.

    Asking "how many copies fall below a line" penalises any diverged family:
    a clean family at 30 % divergence has copies at identity 0.47-0.77 and every
    absolute or proportional cut clips its lower half. But that distribution is
    unimodal - largest internal gap 0.036 - while a genuinely contaminated set is
    bimodal, largest gap 0.247. The shape is divergence-independent; the location
    is not.

    Returns the identity cut separating the two modes, or None if the
    distribution is unimodal, i.e. no contamination.
    """
    d = np.sort(ident)
    n = len(d)
    if n < 20:
        return None
    lo, hi = max(1, int(n * 0.03)), min(n - 1, int(n * 0.97))
    if hi - lo < 3:
        return None
    gaps = np.diff(d[lo:hi + 1])
    j = int(np.argmax(gaps))
    gap = float(gaps[j])
    cut = float((d[lo + j] + d[lo + j + 1]) / 2)
    below = int(np.sum(ident < cut))
    above = n - below
    # BOTH modes must be real. Guarding only the lower one lets a single
    # high-identity outlier stand as an entire "clean family": on Timema
    # SINE_42 (32 copies) the largest gap sat just under a lone copy at
    # identity 1.0, so one copy became the family and the other 31 - which
    # 24 of, by the support threshold, genuinely support the consensus -
    # were called contamination. That collapsed n_supported to 1 and the
    # set scored 0.0 NO_ELEMENT despite a cliff of 0.32.
    if gap < min_gap or min(below, above) < max(3, min_frac * n):
        return None                       # unimodal: a diverged family, not a mixture
    if np.median(d[d >= cut]) - bg < 0.20:
        return None                       # the upper mode is not an element either
    return cut


def support_threshold(ident, bg):
    """A copy supports the consensus if its identity is well above background
    FOR THIS FAMILY. A constant 0.55 flags an old but clean family as
    contaminated - SIMCLEAN at age 0.30 scored 80 with CONTAMINATED - because
    genuine old copies fall below any absolute cut. Scale it instead between the
    measured flank background and the family's own upper quartile.

    The guard matters: in a set with no element q75 is barely above background,
    and a proportional threshold would then pass almost everything. Below a
    minimum contrast there is nothing to support, so nothing is supported."""
    q75 = float(np.percentile(ident, 75))
    contrast = q75 - bg
    if contrast < 0.15:
        return None                      # no element: support is undefined
    return max(bg + 0.10, bg + 0.45 * contrast)
FLANK_SHARE = 0.55          # ungapped flank identity above which two copies share a flank
BG = 0.25                   # unrelated DNA


# Flank-decay readings, computed on 400 bp flanks by flankdecay.py. When a set
# has one, it REPLACES the within-set flank-sharing test: that test compared
# copies over 70 bp and called adjacent similarity "nesting", which produced
# false NESTED_COPIES on hedgehog e2-3 and e2-4. Decay distance distinguishes
# 50 bp of shared context from a satellite that never ends.
try:
    _DECAY = json.load(open("flankdecay.json"))
except Exception:
    _DECAY = {}

# Flank similarity ISLANDS, computed on the SAME 400 bp flanks by
# islands_corpus.py. Sergei saw these by eye in aca_SINE_0 and said they can sit
# far out on either side. They are not what flank_bg measures: flank_bg is one
# average over the whole flank, so a few hundred island columns diluted by a
# thousand ordinary ones vanish into background. aca_SINE_0 reads 0.271 there -
# perfectly normal - while 17 % of its flank is inside an island.
#
# The bands come from his own calls, not from a fit. Of 55 judged sets:
#   <= 0.067  all 31 he called a plain SINE
#   0.19-0.37 all 3 he asked for flank checks on (SINE_NEEDS_GENOME_CHECK,
#             SINE_WITH_CAUTION, GREY_OR_BADLY_PREPARED)
#   >= 0.50   all 7 he rejected or called unusable
# The gaps either side of the middle band are wide, so 0.10 and 0.45 sit in
# empty space rather than on top of any judged set.
try:
    _ISLANDS = json.load(open("flank_islands.json"))
except Exception:
    _ISLANDS = {}

ISLAND_NOTE = 0.10          # above this, say so
ISLAND_STRONG = 0.45        # above this, the flanks are shared outright

# Microsatellite content, from microsat.py. Sergei asked for this by name on
# hyd_SINE_7: "looks like combination of microsatellites to me, not sine - not
# caught by your filters! - needs new criteria on microsatellite content?" and
# on hyd_SINE_17: "surrounded by long 2-nt microsatellites".
#
# It is a genuinely new axis: over the whole 673-set labelled corpus EVERY
# class has a median element-microsatellite content of 0.000, POS and NEGSAT
# and NEGLINEORF alike, with the highest 90th percentile anywhere at 0.076. So
# no existing class exercises it and no threshold can be fitted from them - the
# bands come from his hydra calls: every candidate he accepted is at or below
# 0.065, hyd_SINE_17 is 0.285 and hyd_SINE_7 is 0.468.
#
# The flank threshold is set at 0.20 rather than 0.10 because the synthetic SIM*
# sets carry 0.12-0.17 in their generated flanks - a property of the simulator,
# not of any real locus, and not something to fire on.
try:
    _MSAT = json.load(open("microsat.json"))
except Exception:
    _MSAT = {}

# Region consistency, from region_partition.py: when the copies split into two
# groups, are they the SAME copies in every third of the element?
#
# Sergei separates two things that look alike in any single number. Two real
# subfamilies split the same way along the whole element, and he is explicit
# that this does not matter - "subfamily ambiguity doesn't matter, we answer
# sine/not-sine". A mosaic splits differently in different places: a copy sides
# with one group at the 5' end and another in the middle.
#
# On the corpus the measure runs 0.56 for SIMSUBFAM and up to 0.87 for
# MIXSUBFAM against 0.010-0.012 for NEGMOSAIC and SIMMOSAIC. It is NOT used to
# reject anything: real subfamily mixtures have a median of only 0.064, so any
# threshold that caught hyd_SINE_5 would also catch sets he considers fine. It
# is reported inside the subfamily note so the split can be read for what it is.
try:
    _REGPART = json.load(open("region_part.json"))
except Exception:
    _REGPART = {}

MSAT_ELEM = 0.15
MSAT_FLANK = 0.20


def sat(x, lo, hi):
    return float(min(1.0, max(0.0, (x - lo) / (hi - lo))))


def core_window(el, pres, cons_el, bg, min_len=50, min_gain=0.15, max_frac=0.80):
    """Is the consensus longer than the part the copies actually support?

    Sergei on MIR1_Amn: "weak old sine, consensus needs refinement and
    shortening, otherwise difficult but very real SINE". It scored 0.0
    NO_ELEMENT. Per column its identity runs 0.31 at the 5' end, 0.68-0.71
    across positions 80-140, and 0.29 at the 3' end, while COVERAGE stays
    0.92-1.00 throughout - the copies are present and aligned at the ends, they
    simply do not match there. That is an over-extended consensus, not an absent
    element, and rejecting the set outright is a false negative on a real family.

    The test is NOT "some window scores better than the mean" - in any profile
    with structure, some window always does, and a first version of this flagged
    clean AluJr, AluSc and AluY. The signature is that the EXCLUDED part carries
    essentially no signal: its identity sits at genomic background, meaning the
    consensus runs off the end of the element into unrelated sequence.

    Returns (gain, start, end, core_identity, tail_identity) or None.
    """
    if el.shape[1] < min_len * 2:
        return None
    # el holds base codes, not agreement - compare against the consensus row.
    # A first version summed el directly and produced identities like 48.6
    # instead of fractions, which silently made every result meaningless.
    # Score only columns where the CONSENSUS carries a base. The span lo..hi
    # also contains columns where the consensus itself is gapped - those are
    # insertions in some copies, not positions the copies fail to match, and
    # counting them drags the measured identity far below the truth (0.22 vs
    # 0.51 on MIR1_Amn).
    keep = cons_el != M.GAP
    if keep.sum() < min_len * 2:
        return None
    el = el[:, keep]
    pres = pres[:, keep]
    cons_el = cons_el[keep]
    agree = (el == cons_el[None, :]) & pres
    ok = pres.sum(axis=0)
    with np.errstate(invalid="ignore"):
        colid = np.where(ok > 0, agree.sum(axis=0) / np.maximum(ok, 1), np.nan)
    if not np.isfinite(colid).any():
        return None
    L = el.shape[1]
    whole = float(np.nanmean(colid))
    filled = np.nan_to_num(colid, nan=whole)
    cs = np.concatenate([[0.0], np.cumsum(filled)])
    total = cs[-1]
    best = None
    for wl in range(min_len, int(L * max_frac) + 1, 5):
        for a in range(0, L - wl + 1, 5):
            inside = (cs[a + wl] - cs[a]) / wl
            n_out = L - wl
            if n_out < 20:
                continue
            outside = (total - (cs[a + wl] - cs[a])) / n_out
            if best is None or inside - outside > best[0]:
                best = (inside - outside, a, a + wl, inside, outside)
    if best is None:
        return None
    sep, a, b, inside, outside = best
    # the core must be clearly better AND the discarded part must be at
    # background - otherwise this is ordinary internal variation, not an
    # over-extended boundary.
    if (inside - whole) < min_gain:
        return None
    # The discarded part must carry little signal - but "little" is not "none".
    # On MIR1_Amn the tails sit at 0.395 against a 0.25 background: weak
    # residual similarity, which is what an ancient element's diverged ends
    # actually look like. A 0.12 bound rejected it; 0.18 admits it while
    # staying far below any clean family's outside value (~0.85 for the Alus,
    # ~0.50 for MIR3/MIRc, both of which Sergei confirmed are clean).
    if outside > bg + 0.18:
        return None
    return (round(float(inside - whole), 3), int(a), int(b),
            round(float(inside), 3), round(float(outside), 3))


def parts(path):
    names, A = M.read_aln(path)
    if A is None:
        return None
    k = consensus_index(names)
    cons = A[k]
    nz = np.where(cons != M.GAP)[0]
    if len(nz) < 60:
        return None
    lo, hi = int(nz[0]), int(nz[-1])
    idx = [i for i in range(len(names)) if i != k]
    C = A[idx]
    el = C[:, lo:hi + 1]
    pres = el != M.GAP
    agree = (el == cons[lo:hi + 1][None, :]) & pres
    with np.errstate(invalid="ignore"):
        ident = agree.sum(axis=1) / np.maximum(pres.sum(axis=1), 1)
    ident[pres.sum(axis=1) < 30] = 0.0
    lefts, rights = [], []
    for r in C:
        l = r[:lo]
        lefts.append(l[l != M.GAP][::-1])
        rr = r[hi + 1:]
        rights.append(rr[rr != M.GAP])
    return {"names": [names[i] for i in idx], "el": el, "pres": pres, "ident": ident,
            "cons_el": cons[lo:hi + 1],
            "lefts": lefts, "rights": rights, "cons_bp": len(nz), "n": len(idx)}


def flank_sharing(seqs, minlen=35):
    """Copy x copy flank identity, ungapped, anchored at each copy's own edge.
    Above background means two copies sit in the same genomic context - nested in
    a host repeat, inside a duplication, or in a satellite array."""
    n = len(seqs)
    L = min(70, max((len(s) for s in seqs), default=0))
    if L < minlen or n < 4:
        return np.zeros((n, n)), np.zeros(n, bool)
    Mx = np.full((n, L), 4, dtype=np.int8)
    for i, s in enumerate(seqs):
        m = min(L, len(s))
        Mx[i, :m] = s[:m]
    ok = Mx != 4
    same = np.zeros((n, n))
    for b in range(4):
        m = (Mx == b).astype(np.float32)
        same += m @ m.T
    both = ok.astype(np.float32) @ ok.astype(np.float32).T
    with np.errstate(invalid="ignore", divide="ignore"):
        F = np.where(both >= minlen, same / np.maximum(both, 1), 0.0)
    np.fill_diagonal(F, 0.0)
    return F, (F > FLANK_SHARE).any(axis=1)


def subfamily_split(el, pres, ident, min_grp=12, thr=None):
    """Do the supported copies fall into structurally different groups?

    Copy x copy identity inside the element, split on the leading eigenvector.
    A real family varies smoothly with age; two subfamilies with different
    structure separate into blocks whose between-group identity is clearly below
    their within-group identity."""
    sel = np.where(ident >= (SUPPORT if thr is None else thr))[0]
    if len(sel) < 2 * min_grp:
        return None
    E = el[sel]
    P = pres[sel]
    ok = P.astype(np.float32)
    same = np.zeros((len(sel), len(sel)), np.float32)
    for b in range(4):
        m = ((E == b) & P).astype(np.float32)
        same += m @ m.T
    both = ok @ ok.T
    with np.errstate(invalid="ignore", divide="ignore"):
        S = np.where(both >= 40, same / np.maximum(both, 1), np.nan)
    if not np.isfinite(S).any():
        return None
    S = np.where(np.isfinite(S), S, np.nanmean(S))
    Sc = S - S.mean(axis=0, keepdims=True)
    try:
        w, v = np.linalg.eigh(Sc @ Sc.T)
    except np.linalg.LinAlgError:
        return None
    g = v[:, -1] > np.median(v[:, -1])
    if g.sum() < min_grp or (~g).sum() < min_grp:
        return None
    within = np.nanmean([S[np.ix_(g, g)].mean(), S[np.ix_(~g, ~g)].mean()])
    between = S[np.ix_(g, ~g)].mean()
    gap = float(within - between)

    # Null: split on a RANDOM vector instead of the leading eigenvector. Any
    # split of a homogeneous family yields a positive gap whose size grows with
    # divergence and shrinks with n - which is why the raw gap tracked age at
    # +0.99 and copy number at -0.97 in families with no subfamily structure at
    # all. Subtracting the null removes both dependencies at their source.
    rng = np.random.default_rng(0)
    null = []
    for _ in range(24):
        r = rng.normal(size=len(sel))
        gr = r > np.median(r)
        if gr.sum() < 4 or (~gr).sum() < 4:
            continue
        w = np.nanmean([S[np.ix_(gr, gr)].mean(), S[np.ix_(~gr, ~gr)].mean()])
        null.append(w - S[np.ix_(gr, ~gr)].mean())
    nullmean = float(np.mean(null)) if null else 0.0
    return {"gap": gap, "gap_excess": float(gap - nullmean), "gap_null": nullmean,
            "sizes": [int(g.sum()), int((~g).sum())],
            "within": float(within), "between": float(between),
            "members": [int(sel[i]) for i in np.where(g)[0]]}


def verdict(path):
    p = parts(path)
    if p is None:
        return None
    n, ident = p["n"], p["ident"]
    dec = _DECAY.get(os.path.basename(path).replace(".aln.fa", ""))
    isl = _ISLANDS.get(os.path.basename(path).replace(".aln.fa", ""))
    ms = _MSAT.get(os.path.basename(path).replace(".aln.fa", ""))
    rp = _REGPART.get(os.path.basename(path).replace(".aln.fa", ""))

    # background first, because the support threshold is defined against it
    # A side only contributes if it actually HAS flank sequence to measure. The
    # curated Timema `subfam` alignments are element-only - median flank length
    # is 0 bases and 193 of 200 copies have no left flank at all - but a handful
    # of copies do carry a few bases, those few happen to match, and the
    # background came out as 1.000. With background above element identity there
    # is no cliff by definition, so t3-1 and t3-2 scored 0.0 NO_ELEMENT while
    # top100 and rand100 of the SAME subfamilies scored 100.0.
    #
    # Same principle as the uniqueness fallback and the over-extended consensus:
    # a measurement that could not be made must not be scored as evidence.
    MIN_FLANK_PAIRS = 20
    MIN_FLANK_BP = 15
    fl_pre = []
    flank_measured = True
    for seqs in (p["lefts"], p["rights"]):
        usable = [q for q in seqs if len(q) >= MIN_FLANK_BP]
        if len(usable) < MIN_FLANK_PAIRS:
            continue                      # not enough flank on this side to say anything
        F0, _ = flank_sharing(seqs)
        vv = F0[np.triu_indices(len(seqs), 1)]
        vv = vv[vv > 0]
        if len(vv) >= MIN_FLANK_PAIRS:
            fl_pre.append(float(np.median(vv)))
    if fl_pre:
        bg_meas = float(np.mean(fl_pre))
    else:
        bg_meas = BG                      # unrelated-DNA constant; nothing measurable
        flank_measured = False

    # ...but when the consensus stops short of the element, the sequence just
    # outside it IS the element, so bg_meas is measuring the element against
    # itself. hyd_SINE_0 came out at 0.614 that way, above its own element
    # identity of 0.560, and every copy then failed to reach the support
    # threshold: 0 of 100 supported, NO_ELEMENT, on a family Sergei reads as a
    # good SINE and RepBase calls SINE2-2B_HM.
    #
    # The decay profile walks outward in 25 bp steps, so its far end is
    # unrelated DNA no matter where the consensus ends. Prefer it.
    flank_bg_far = None
    if dec:
        far = []
        for side in ("prof_L", "prof_R"):
            prof = dec.get(side) or []
            tail = [x for x in prof[-4:] if x is not None]
            if tail:
                far.append(float(np.mean(tail)))
        if far:
            flank_bg_far = float(np.mean(far))
    bg_true = bg_meas if flank_bg_far is None else max(BG, flank_bg_far)
    # Two separate questions, previously conflated into one threshold:
    #   is there an element at all      -> support_threshold
    #   is this set a MIXTURE           -> contamination_split
    thr = support_threshold(ident, bg_true)
    cut = contamination_split(ident, bg_true)
    if thr is None:
        supported = np.zeros(n, bool)     # no element: nothing supports it
    elif cut is not None:
        supported = ident >= cut          # a real mixture: split at the mode boundary
    else:
        supported = np.ones(n, bool)      # unimodal: every copy belongs, however diverged
    n_sup = int(supported.sum())

    cw = core_window(p["el"], p["pres"], p["cons_el"], bg_meas)
    if cw:
        # Support was being judged against a boundary this very test has just
        # shown to be wrong, which is circular: MIR1_Amn had 0 of 100 copies
        # "supported" and scored 0.0 NO_ELEMENT, on a family Sergei calls
        # "difficult but very real". Re-measure per-copy identity over the core
        # window only, and judge support there.
        _keep = p["cons_el"] != M.GAP
        _el = p["el"][:, _keep][:, cw[1]:cw[2]]
        _pr = p["pres"][:, _keep][:, cw[1]:cw[2]]
        _ce = p["cons_el"][_keep][cw[1]:cw[2]]
        _ag = (_el == _ce[None, :]) & _pr
        with np.errstate(invalid="ignore"):
            _id = _ag.sum(axis=1) / np.maximum(_pr.sum(axis=1), 1)
        _id[_pr.sum(axis=1) < 20] = 0.0
        if float(np.mean(_id)) > float(np.mean(ident)):
            ident = _id
            thr = support_threshold(ident, bg_meas)
            cut = contamination_split(ident, bg_meas)
            supported = ident >= (cut if cut is not None else (thr if thr is not None else SUPPORT))
            n_sup = int(supported.sum())
    v_cut = cut

    FL, sharedL = flank_sharing(p["lefts"])
    FR, sharedR = flank_sharing(p["rights"])
    shared = sharedL | sharedR
    n_shared = int(shared.sum())
    global_elev = float(np.median(np.concatenate([FL[FL > 0], FR[FR > 0]]))) \
        if (FL > 0).any() or (FR > 0).any() else 0.0

    # a "core" copy: supports the consensus AND sits in unique genomic sequence
    core = supported & ~shared
    n_core = int(core.sum())

    # Structure is tested on the SUPPORTED copies only. Contamination masks a
    # scramble - seg_diag against scrambling went from -0.61 on the joint grid
    # to -0.92 once contamination was low - so pruning has to come first.
    sub = subfamily_split(p["el"], p["pres"], ident, thr=thr)

    # Is this one family with subfamily structure, or two different things that
    # should be separated and re-analysed before any verdict is given?
    #
    # Sergei on AluYk11: "truly a mixture, top proper about half longer
    # similarity sequences need separate analysis as a good SINE candidate, then
    # bottom more discordant shorter ones need additional fresh re-run of whole
    # analysis. This non-homogenous selection problem is different from bad SINE
    # set and needs additional work before verdict."
    #
    # The discriminator is the LENGTH of the two groups, not their sequence
    # divergence. Measured over every set carrying SUBFAMILY_NOTE, relative
    # median length difference put his three flagged sets at the top - AluYk11
    # 0.309, CAS 0.133, AluYh9 0.122 - with a clear drop to 0.096 and below for
    # every set he did not flag. Ordinary subfamilies differ in sequence at
    # roughly equal length; a mixture differs in length.
    # Which copies form the two groups matters as much as the measure. The
    # subfamily split is driven by sequence structure; on CAS and AluYh9 - both
    # of which Sergei called mixtures - it groups them so that the two halves
    # differ by only 0.013 and 0.012 in identity, the LOWEST of all nine cases.
    # Splitting instead on each copy's own identity, length and coverage
    # recovers 0.209 and 0.128 and matches his reading on all nine.
    het = None
    _rg = None
    try:
        import test2_props as _T2
        _rg = _T2.row_groups(path)
    except Exception:
        _rg = None
    if _rg and _rg.get("rowsplit_dident") is not None:
        _did = float(_rg["rowsplit_dident"])
        if _did >= 0.126 and sub and "sizes" in sub:
            het = {"d_ident": round(_did, 3),
                   "rel_len": round(float(_rg.get("rowsplit_dlen", 0)) /
                                    max(1.0, float(_rg.get("rowsplit_dlen", 1))), 3),
                   "med": [int(_rg.get("rowsplit_dlen", 0)), 0],
                   "sizes": sub["sizes"]}
    if False and sub and "members" in sub:
        _el = (p["el"] != M.GAP).sum(axis=1).astype(float)
        _m = np.zeros(len(_el), bool)
        _m[sub["members"]] = True
        if _m.sum() >= 5 and (~_m).sum() >= 5:
            _a, _b = float(np.median(_el[_m])), float(np.median(_el[~_m]))
            _rel = abs(_a - _b) / max(_a, _b, 1.0)
            # Group LENGTH difference was the original test. Measured against
            # nine of Sergei's own judgements it does not work: he called
            # s8_225seqs "not mixture" at 0.228, higher than three of the four
            # sets he DID call mixtures (0.122-0.138). The ranges overlap
            # completely, so no threshold exists. The measure is invalid, not
            # miscalibrated.
            #
            # What separates them is how much the two groups differ in
            # IDENTITY, which is also what he actually describes: "top proper
            # about half longer similarity sequences ... bottom more discordant".
            # On the same nine: his mixtures 0.128-0.215, his non-mixtures
            # 0.012-0.124. No overlap.
            _ia = float(np.mean(ident[_m])) if _m.any() else 0.0
            _ib = float(np.mean(ident[~_m])) if (~_m).any() else 0.0
            _did = abs(_ia - _ib)
            if _did >= 0.126:
                het = {"d_ident": round(_did, 3), "rel_len": round(_rel, 3),
                       "med": [round(_a), round(_b)], "sizes": sub["sizes"]}
    turb = turbulence(p["el"], p["pres"])

    elen = (p["el"] != M.GAP).sum(axis=1).astype(float)
    # length spread over the SUPPORTED copies, not the core: excluding copies for
    # sharing flanks would hide truncation that the set really has
    cv_core = float(np.std(elen[supported]) / max(1e-9, np.mean(elen[supported])))         if n_sup >= 8 else 1.0
    id_core = float(np.median(ident[core])) if n_core else 0.0

    # Same guard as bg_meas above - this is the value that actually feeds
    # g_elem, so leaving it unguarded reproduced the whole failure: an
    # element-only alignment yielded flank_bg 1.000, which makes the cliff
    # (id_all - flank_bg) negative and forces NO_ELEMENT on a perfectly good
    # family. Both computations must refuse to report a background they could
    # not measure.
    fl_id = []
    for seqs in (p["lefts"], p["rights"]):
        usable = [q for q in seqs if len(q) >= MIN_FLANK_BP]
        if len(usable) < MIN_FLANK_PAIRS:
            continue
        F, _ = flank_sharing(seqs)
        v = F[np.triu_indices(len(seqs), 1)]
        v = v[v > 0]
        if len(v) >= MIN_FLANK_PAIRS:
            fl_id.append(float(np.median(v)))
    flank_bg = float(np.mean(fl_id)) if fl_id else BG

    bg_cliff = bg_true

    tsd = [M.find_tsd(p["lefts"][i][::-1], p["rights"][i])
           for i in np.where(core)[0][:120]]
    tsd_frac = float(np.mean([t > 0 for t in tsd])) if tsd else 0.0

    # ---- evidence groups, each 0-1 --------------------------------------
    # The cliff is measured over ALL copies, not over the core. Measuring it on
    # the copies that were selected for matching the consensus is circular - it
    # is how random loci first scored 0.67 here, and the same trap that made the
    # gap-based pruning rule manufacture families out of noise.
    # mean, not median: with 30 % contamination the median is still a real copy,
    # so a median cliff cannot see contamination at all. The score must describe
    # the set as submitted; recoverability is carried by n_core and the flags.
    id_all = float(np.mean(ident))

    # An over-extended consensus depresses id_all across its whole length and
    # can drive a real family to NO_ELEMENT. When a contiguous window is much
    # better supported than the whole span, judge the element on that window and
    # report the coordinates so the consensus can be trimmed.
    # "Is there an element" is a BINARY question, so the evidence must saturate
    # once it is answered. A cliff of 0.28 means copies match the consensus at
    # 0.57 while unrelated DNA matches at 0.29 - already unambiguous. Scoring on
    # to 0.55 made the statistic report youth, not existence, and cost a clean
    # 30 %-divergent family a third of its element score.
    g_elem = (sat(id_all - bg_cliff, 0.10, 0.30) ** 0.5) *              (sat(n_sup / max(1, n), 0.15, 0.98) ** 0.5)
    # Subfamily structure is deliberately NOT scored. Sergei: "subfamily
    # ambiguity doesn't matter, we answer sine/not-sine, not this sine or that
    # sine." Two subfamilies in one set is still a set of SINEs. It stays as a
    # note so the copies can be split later, but it must not cost score.
    # Length spread is softened deliberately. Truncated copies are still copies
    # of that family - Sergei on the truncated sets: "more like a SINE", "looks
    # more like a LINE", never "not an element". Scored on a 0.10-0.50 range it
    # notes heterogeneity without destroying the verdict; at 0.05-0.25 it drove
    # the truncated class to a mean of 9.7.
    g_homog = 1 - sat(cv_core, 0.10, 0.50)
    uniq_measured = True
    if dec:
        # BOTH parts of the decay curve matter. Distance alone rates a LINE
        # fragment as isolated, because its shared flank runs only ~50-75 bp
        # before the copies truncate - yet identity right at the edge is 0.89,
        # which says the element does not end there at all.
        g_uniq = (1 - sat(dec["decay_max"], 75, 300)) *                  (1 - sat(dec.get("edge_max") or 0.0, 0.45, 0.85))
    else:
        # No 400 bp decay data. The within-set flank-sharing test that used to
        # stand in here is measured over ~50-70 bp and is documented above as
        # producing false nesting calls at that width. Measured across the whole
        # corpus it turned out to carry NO discriminating power: every genuine
        # negative (NEGRAND, NEGSAT, NEGSEGDUP, NEGLINEORF, 30 sets) scores 0.0
        # both with and without it - the element gate does all the rejecting.
        # Its only measurable effect was to suppress real SINEs, including three
        # Timema candidates an expert confirmed by eye (38.6/29.3/71.3 ->
        # 61.9/44.5/100.0). An absent measurement must not be scored as guilt,
        # so uniqueness is neutral here and the absence is reported instead.
        g_uniq = 1.0
        uniq_measured = False
    g_ins = sat(tsd_frac, 0.10, 0.60)

    # The element group GATES the score rather than voting in it. As a co-equal
    # geometric term, a set with no element (g_elem 0.31) still scored 56 because
    # its flanks were unique and its lengths uniform - both true and both
    # irrelevant when there is nothing there. Uniqueness and homogeneity only
    # mean something once an element exists.
    rest = (max(g_homog, 1e-6) ** 0.55) * (max(g_uniq, 1e-6) ** 0.45)
    score = 100 * min(1.0, g_elem * rest * (1 + W["insertion"] * g_ins))

    # ---- sub-cases, not one "negative" bucket ---------------------------
    flags = []
    if dec:
        call = dec.get("call")
        if call == "SATELLITE_OR_DUPLICATION":
            flags.append({"code": "NOT_ISOLATED", "n": dec["decay_max"],
                          "text": "Similarity between copies continues %d bp beyond the "
                                  "element and does not reach background — a satellite "
                                  "array or a duplicated region, not independent "
                                  "insertions." % dec["decay_max"]})
        elif call == "ELEMENT_CONTINUES":
            flags.append({"code": "FRAGMENT_OF_LONGER", "n": dec["decay_max"],
                          "text": "The sequence flanking these copies is itself shared "
                                  "(identity %.2f at the edge, reaching background only "
                                  "%d bp out). They look like fragments of a longer "
                                  "element — a LINE 3' end or similar — rather than a "
                                  "short element with its own boundaries."
                                  % (dec["edge_max"], dec["decay_max"])})
        elif call == "ADJACENT_SIMILARITY":
            flags.append({"code": "ADJACENT_SIMILARITY", "n": dec["decay_max"],
                          "text": "Copies share about %d bp of sequence just outside the "
                                  "element, then reach background. Worth noting, but they "
                                  "are in independent locations."
                                  % dec["decay_max"]})
    # Flank sharing is reported whatever the core count. Previously this sat
    # behind "n_core >= 8", so a set whose copies ALL share flanking sequence
    # scored near zero with no explanation at all - hedgehog e1-4 and e2-2 both
    # did exactly that, at flank identity 0.68 against 0.28 for a normal family.
    # ...but only when there were flanks to measure. On an element-only
    # alignment (the curated `subfam` variants carry no genomic flanks at all)
    # the few copies that do have a stub of sequence match each other trivially,
    # and this fired as if the loci were inside a satellite.
    if flank_measured and (not dec) and global_elev - BG > 0.15 and n >= 15:
        flags.append({"code": "SHARED_FLANKS", "n": n_shared,
                      "text": "Flanking sequence is %.2f identical between copies against "
                              "%.2f for unrelated DNA — these loci are not in independent "
                              "genomic contexts. A satellite array, a duplicated region, or "
                              "copies inside one host repeat." % (global_elev, BG)})
    # Reported, never subtracted from the score. All three sets this band was
    # calibrated on he still leaned towards SINE on - "requires post-processing
    # with proving uniqness of at least some flanks on whole-genome level" is a
    # next step, not a rejection. Penalising the score would have turned e2-4
    # (97.4, which he accepts with a caveat) into a reject.
    if isl and isl.get("frac", 0) >= ISLAND_NOTE:
        frac, ncol, nisl = isl["frac"], isl["cols"], isl["islands"]
        already = flank_bg >= 0.40
        if frac >= ISLAND_STRONG:
            txt = ("%.0f %% of the flank sits inside patches of raised similarity "
                   "(%d columns in %d patches). At this level the copies are not in "
                   "independent places at all." % (100 * frac, ncol, nisl))
        elif already:
            txt = ("%.0f %% of the flank sits inside patches of raised similarity "
                   "(%d columns in %d patches), which agrees with the raised flank "
                   "average of %.2f." % (100 * frac, ncol, nisl, flank_bg))
        else:
            txt = ("%.0f %% of the flank sits inside patches of raised similarity "
                   "(%d columns in %d patches), yet the flank average is an "
                   "ordinary %.2f. The similarity is localised, so averaging hides "
                   "it. Some copies are probably inside a shared larger repeat or "
                   "duplication: worth checking those flanks against the whole "
                   "genome before treating every copy as an independent insertion."
                   % (100 * frac, ncol, nisl, flank_bg))
        flags.append({"code": "FLANK_ISLANDS", "n": [ncol, nisl, frac], "text": txt})

    # The mirror of CONSENSUS_OVEREXTENDED, which had no counterpart until
    # hyd_SINE_0. When the similarity carries on past the consensus edge but
    # ENDS within a couple of hundred bases, the copies are not inside a
    # satellite and are not fragments of a LINE - the consensus simply stops
    # too early, and the missing piece is about as long as the decay distance.
    # A LINE fragment runs on for thousands of bases; this does not.
    # ...and the flank beyond that stretch must be ORDINARY. Without this guard
    # the reason also fired on four NEGLINEORF sets and a segmental duplication,
    # whose consensuses are 260-300 bp and whose decay also reads under 200 bp.
    # What separates them is what the rest of the flank looks like: an
    # under-extended consensus has a normal flank once past the missing piece
    # (hyd_SINE_0 has 6 % of its flank in islands), while a LINE ORF or a
    # duplication has a flank that is shared all the way out (37-82 %).
    if (dec and dec.get("call") == "ELEMENT_CONTINUES"
            and dec.get("decay_max", 0) <= 200
            and flank_bg_far is not None and flank_bg_far < 0.40
            and (isl is None or isl.get("frac", 0) < 0.15)):
        flags.append({"code": "CONSENSUS_UNDEREXTENDED",
                      "n": dec["decay_max"],
                      "text": "Similarity carries on about %d bp past the end of the "
                              "consensus and then reaches background (%.2f), so that "
                              "stretch is the element rather than its context. The "
                              "consensus is too short by roughly that much - extend it "
                              "and re-run. Until then the flank average is being "
                              "measured on element sequence, which understates "
                              "everything that depends on it."
                              % (dec["decay_max"], flank_bg_far)})

    if ms and ms.get("msat_elem", 0) >= MSAT_ELEM:
        flags.append({"code": "MICROSATELLITE_ELEMENT", "n": ms["msat_elem"],
                      "text": "%.0f %% of the element is simple tandem repeat in the "
                              "median copy (%.0f %% in the top tenth). Every family "
                              "Sergei accepted in hydra sits at or below 7 %%, and no "
                              "class in the labelled corpus reaches 1 %%. A consensus "
                              "built mostly out of microsatellite will collect loci by "
                              "their repeat content, not by common descent."
                              % (100 * ms["msat_elem"], 100 * ms.get("msat_elem_p90", 0))})

    if ms and (ms.get("msat_flank") or 0) >= MSAT_FLANK:
        flags.append({"code": "MICROSATELLITE_FLANK", "n": ms["msat_flank"],
                      "text": "%.0f %% of the flanking sequence is simple tandem repeat "
                              "in the median copy. The copies are sitting in "
                              "microsatellite tracts, so the flanks cannot show whether "
                              "the insertions are independent and short-read mapping "
                              "into such tracts is unreliable."
                              % (100 * ms["msat_flank"])})

    # hyd_SINE_6 scored 100.0 with 25 of 100 copies in the core, and Sergei read
    # it as "mosaic columns in flanks, few copies, not firm edges, unknown
    # middle part". n_core was already in the output but nothing said it out
    # loud, so a set held up by a quarter of its copies looked identical to one
    # held up by all of them.
    # n_core == 0 used to slip past this entirely (the old `n_core >= 1` guard
    # excluded exactly the worst case: a family with a small core still hit
    # the flag and its cap, but the family with NO core at all -- SINE16 (oma),
    # 0 of 100 -- did not). Zero is not a special case that deserves LESS
    # scrutiny than one; if anything it deserves more. Found 2026-09-10.
    if n >= 20 and n_core / float(n) < 0.35:
        flags.append({"code": "SMALL_CORE", "n": [n_core, n],
                      "text": "Only %d of %d copies form the core the consensus is "
                              "actually built on (%.0f %%). The rest are too divergent "
                              "or too truncated to support it, so the score describes a "
                              "minority of what was collected. Look at whether the other "
                              "%d belong here at all."
                              % (n_core, n, 100.0 * n_core / n, n - n_core)})

    if g_elem < 0.25:
        flags.append({"code": "NO_ELEMENT",
                      "text": "No element: too few copies support the consensus above background."})
    if n_core >= 8:
        if cut is not None and (n - n_sup) / max(1, n) > 0.03:
            flags.append({"code": "CONTAMINATED", "n": int(n - n_sup),
                          "text": "%d of %d copies do not support the consensus; prune them and "
                                  "re-score (see prune.py)." % (n - n_sup, n)})
        if sub and sub.get("gap_excess", sub["gap"]) > 0.03:
            flags.append({"code": "SUBFAMILY_NOTE", "n": sub["sizes"],
                          "text": "Supported copies split into structurally different groups of "
                                  "%d and %d (within-group identity %.2f, between %.2f, excess "
                                  "over a random split %.3f). Both are SINEs - this does not "
                                  "affect the verdict, but the copies can be split if the "
                                  "subfamilies are wanted separately.%s"
                                  % (sub["sizes"][0], sub["sizes"][1], sub["within"],
                                     sub["between"], sub.get("gap_excess", 0),
                                     "" if not rp else
                                     (" The same copies stay together along the whole element "
                                      "(consistency %.2f), so this is two subfamilies."
                                      % rp["consistency"] if rp["consistency"] >= 0.20 else
                                      " But group membership shuffles between the 5' end, the "
                                      "middle and the 3' end (consistency %.2f, against 0.56 for "
                                      "a real subfamily split): the copies are not two families "
                                      "so much as patchwork, which is worth looking at by eye."
                                      % rp["consistency"]))})
        if (not dec) and n_shared / max(1, n) > 0.15:
            kind = ("a satellite or segmental duplication" if global_elev - BG > 0.15
                    else "another repeat or a duplication")
            flags.append({"code": "FLANKS_UNMEASURED", "n": n_shared,
                          "text": "%d of %d copies appear to share flanking sequence, which "
                                  "would put them inside %s - but this set has no 400 bp "
                                  "flank-decay data, and over the ~50-70 bp actually present "
                                  "that reading is not reliable: adjacent context and true "
                                  "nesting look identical at this width. Uniqueness is scored "
                                  "as neutral rather than penalised, and this is reported as "
                                  "an unmade measurement, not as evidence against the family. "
                                  "Re-run with 400 bp flanks to settle it."
                                  % (n_shared, n, kind)})
    if cw:
        flags.append({"code": "CONSENSUS_OVEREXTENDED", "n": [cw[1], cw[2]],
                      "text": "The copies support consensus positions %d-%d (identity %.2f) "
                              "much better than the consensus as a whole (%.2f). Coverage is "
                              "even across the full length, so the copies are present at the "
                              "ends and simply do not match there - the consensus is longer "
                              "than the element it describes. Trim it to that window and "
                              "re-run; the family is real, the boundary is not."
                              % (cw[1], cw[2], cw[3], cw[3] - cw[0])})
    if n < 30:
        flags.append({"code": "INSUFFICIENT_COPIES", "n": n,
                      "text": "Only %d copies. The subfamily-split and flank-decay tests need "
                              "about 30 before they mean anything, so this score is an "
                              "impression rather than a measurement - it should not be read as "
                              "a confident accept OR a confident reject. Gather more copies "
                              "from the genome before deciding." % n})
    if het:
        flags.append({"code": "HETEROGENEOUS_SELECTION", "n": het["sizes"],
                      "text": "This looks like two different things collected together: the "
                              "%d and %d copies split into groups differing by %.2f in identity "
                              "to the consensus (his mixtures run 0.13-0.22, single families "
                              "0.01-0.12). Median lengths are %d and %d bp. Separate them and "
                              "re-run each part - the score is withheld, because scoring a "
                              "mixture answers neither question."
                              % (het["sizes"][0], het["sizes"][1], het["d_ident"],
                                 het["med"][0], het["med"][1])})
    if not flank_measured:
        flags.append({"code": "NO_FLANKS_PRESENT", "n": 0,
                      "text": "This alignment carries no genomic flanks, so nothing outside the "
                              "element could be measured: the background is assumed rather than "
                              "observed, and the flank-based tests (isolation, nesting, shared "
                              "context, TSD) are simply not available. The element itself is "
                              "still judged normally."})
    if cv_core > 0.18 and n_core >= 15:
        flags.append({"code": "TRUNCATED_COPIES", "n": round(cv_core, 3),
                      "text": "Copy length varies widely (CV %.2f against 0.05 for a clean "
                              "family). Many copies are 5'-truncated. Still the same family, "
                              "but the truncated copies may be worth setting aside."
                              % cv_core})
    # NEEDS_REVIEW (turbulence > 0.12) was REMOVED 2026-09-01.
    #
    # Turbulence counts COLUMNS dipping 0.25 below the element's own median
    # identity. Two things were wrong with it. First, the 0.25 is absolute while
    # the available range is not: at median identity 0.88 a column only has to
    # fall to 0.63, which ordinary variation does constantly, whereas at 0.55 it
    # would have to reach 0.30 - background - which essentially never happens.
    # So it was most sensitive exactly where families are most homogeneous and
    # blind where they are least: the old Alus scored 0.125-0.162 and every MIR
    # scored 0.000-0.008.
    #
    # Second and more fundamental: Sergei reviewed these by eye and every real
    # problem he identified was about a SUBSET OF COPIES - "11 lower sequences",
    # "top half vs bottom", "2 variants", "11 in the middle". Those are row-wise.
    # Turbulence is column-wise. It was answering a question nobody asked, and
    # the row-wise detectors (CONTAMINATED, SUBFAMILY_NOTE) already caught every
    # set he flagged.
    #
    # Measured cost of keeping it: 129 sets flagged, 91 with it as the ONLY
    # flag, all scoring 90.1-100, including 8 POS and 60 SQ known-good sets. It
    # never moved a score. It was noise on high-confidence positives, and he
    # named 10 clean families it had wrongly flagged.

    if n_core >= 20 and score < 55:
        flags.append({"code": "RECOVERABLE_CORE", "n": n_core,
                      "text": "Despite the problems above, %d copies look like genuine elements. "
                              "That is enough to build a clean family from — worth following up "
                              "even though the set as submitted scores poorly." % n_core})

    # ---- reasons that CAP the score --------------------------------------
    #
    # Until now a reason could fire and change nothing. hyd_SINE_7 is 52 % simple
    # repeat and scored 77; hyd_SINE_6 rests on 25 of 100 copies and scored 93;
    # hyd_SINE_2 has 8 copies and scored 95 while its own flag text said "this
    # score is an impression rather than a measurement". Sergei rejected all
    # three. A number that reads as acceptance while the text underneath says
    # the evidence is not there is worse than no number.
    #
    # These three do not reduce the score by some amount - there is no principled
    # amount. They cap it just below the acceptance line, which states exactly
    # what is true: on this evidence the question cannot be answered yes.
    # Anything the copies genuinely show is still in the groups and the flags.
    # INSUFFICIENT_COPIES was in this list and should never have been. It says
    # the set cannot be measured, not that the evidence points away from a SINE,
    # and capping on it turned HUM__hum__AluYb9 - a known Alu subfamily Sergei
    # calls "no problems good sine" - into a 45. That is the exact mistake
    # already written down as a rule: an absent measurement must never be scored
    # as guilt. It now sets a separate state instead of moving the score.
    capped_by = []
    for f in flags:
        if f["code"] in ("MICROSATELLITE_ELEMENT", "SMALL_CORE"):
            capped_by.append(f["code"])
    if capped_by and score > 45.0:
        score = 45.0

    # "cannot be assessed" is a third answer, not a low one. The score is still
    # reported, but it describes too little evidence to stand on its own.
    assessable = not any(f["code"] in ("INSUFFICIENT_COPIES", "NO_FLANKS_PRESENT")
                         for f in flags)

    return {"set": os.path.basename(path).replace(".aln.fa", ""),
            "score": round(float(score), 1),
            "deferred": bool(het),
            "flank_measured": bool(flank_measured),
            "uniqueness_measured": bool(uniq_measured),
            "groups": {"element": round(float(g_elem), 3),
                       "homogeneity": round(float(g_homog), 3),
                       "uniqueness": round(float(g_uniq), 3),
                       "insertion": round(float(g_ins), 3)},
            "weights": W,
            "n": n, "n_supported": n_sup, "n_core": n_core, "n_shared": n_shared,
            "core_identity": round(id_core, 3), "all_identity": round(float(np.median(ident)), 3), "flank_bg": round(flank_bg, 3),
            "core_len_cv": round(cv_core, 3), "tsd_frac": round(tsd_frac, 3),
            "flank_bg_far": None if flank_bg_far is None else round(flank_bg_far, 3),
            "msat_elem": None if not ms else ms.get("msat_elem"),
            "msat_flank": None if not ms else ms.get("msat_flank"),
            "core_frac": round(n_core / float(max(1, n)), 3),
            "region_consistency": None if not rp else rp.get("consistency"),
            "capped_by": capped_by,
            "assessable": assessable,
            "island_frac": None if not isl else isl.get("frac"),
            "island_cols": None if not isl else isl.get("cols"),
            "support_threshold": None if thr is None else round(float(thr), 3),
            "contamination_cut": None if v_cut is None else round(float(v_cut), 3),
            "turbulence": None if turb is None else round(turb, 3),
            "subfamily": sub and {k: v for k, v in sub.items() if k != "members"},
            "flags": flags}


def _one(f):
    try:
        return verdict(f)
    except Exception as exc:
        return {"set": os.path.basename(f).replace(".aln.fa", ""), "error": str(exc)}


def main():
    import multiprocessing as mp
    d = sys.argv[1] if len(sys.argv) > 1 else "aln_v2"
    out = sys.argv[2] if len(sys.argv) > 2 else "verdicts.json"
    files = sorted(glob.glob(os.path.join(d, "*.aln.fa")))
    files = [f for f in files if not os.path.basename(f).startswith("SQ__")]
    res = {}
    with mp.Pool(max(1, mp.cpu_count() - 2)) as pool:
        for v in pool.imap_unordered(_one, files, chunksize=4):
            if v:
                res[v["set"]] = v
    json.dump(res, open(out, "w"), separators=(",", ":"))

    from collections import defaultdict
    byc = defaultdict(list)
    for s, v in res.items():
        if "error" not in v:
            byc[s.split("__")[0]].append(v)
    print("%-12s %5s %7s   %-38s %s" % ("class", "sets", "score", "evidence groups", "n_core"))
    print("-" * 86)
    for c in ("POS", "MIXED10", "MIXED30", "NEGJITTER", "NEGSPLICE", "NEGTRUNC5", "NEGRAND"):
        v = byc.get(c, [])
        if not v:
            continue
        g = lambda k: np.mean([x["groups"][k] for x in v])
        print("%-12s %5d %7.1f   el %.2f  hom %.2f  uniq %.2f  ins %.2f   %5.1f"
              % (c, len(v), np.mean([x["score"] for x in v]),
                 g("element"), g("homogeneity"), g("uniqueness"), g("insertion"),
                 np.mean([x["n_core"] for x in v])))
    print("\nflag frequency by class")
    for c in ("POS", "MIXED10", "MIXED30", "NEGJITTER", "NEGSPLICE", "NEGTRUNC5", "NEGRAND"):
        v = byc.get(c, [])
        if not v:
            continue
        cnt = defaultdict(int)
        for x in v:
            for f in x["flags"]:
                cnt[f["code"]] += 1
        print("  %-12s %s" % (c, dict(cnt) or "clean"))


def full_verdict(path_top100, path_flank_L=None, path_flank_R=None):
    """verdict() plus the edge_profile.py status/recommendation layer
    (Sergei, 2026-09-10): every verdict must say not just a score but
    whether the analysis is DONE (status "FINAL") or needs more work
    (status "NEEDS_ADDITIONAL_WORK", with concrete steps) -- flank
    extension, subgrouping, mosaicism resolution, etc.

    path_flank_L/R are long-flank (up to edge_profile.MAX_OFF, typically
    1000bp) validation alignments anchored on the current consensus --
    NOT the short 50L/70R display alignment (verdict()'s own flank_sharing
    test already knows it can't see far enough; this is what CAN, when
    those files exist). Pass None for a side that hasn't been built yet;
    the function still returns the verdict() score, just without an edge
    classification for that side (never silently assumed clean).
    """
    v = verdict(path_top100)
    if v is None:
        return {"verdict": None, "status": "NEEDS_ADDITIONAL_WORK",
                "recommended_next_steps": ["no element found in top100 -- verdict() itself failed"]}

    try:
        import edge_profile as E
    except ImportError:
        return {"verdict": v, "edge": None, "status": "FINAL" if not v.get("flags") else "NEEDS_ADDITIONAL_WORK",
                "recommended_next_steps": ["edge_profile.py not available -- edge pattern not checked"]}

    edge_classifications = []
    for side, path in (("L", path_flank_L), ("R", path_flank_R)):
        if path:
            try:
                edge_classifications.append(E.classify(path, side))
            except Exception as exc:
                edge_classifications.append({"side": side, "pattern": "ERROR", "error": str(exc)})

    edge_status, edge_steps = E.recommend([c for c in edge_classifications if c.get("pattern") not in (None, "ERROR")])

    # combine with verdict()'s own flags into one status/steps list, plus a
    # human label describing WHAT this looks like when the picture is clear
    # enough to say (the "clean SINE" / "LINE-like rugged edge" / "needs
    # subgrouping" distinction asked for) -- never invented when ambiguous.
    steps = list(edge_steps)
    if v.get("flags"):
        for f in v["flags"]:
            steps.append("verdict flag %s: %s" % (f["code"], f.get("text", "")[:160]))

    patterns = {c.get("pattern") for c in edge_classifications if c.get("pattern")}
    score = v.get("score", 0)
    if not steps and score >= 90:
        label = "clean SINE (independent flanks, single sharp edge, high score)"
    elif "GRADED_DECAY" in patterns:
        label = "LINE-like: at least one edge decays gradually rather than a sharp SINE boundary"
    elif "MULTI_GROUP" in patterns:
        label = "mixed population: edge groups suggest real subgroup structure, not one family"
    elif score < 55:
        label = "not SINE-like on current evidence"
    else:
        label = "grey zone: score and edge pattern do not cleanly agree"

    status = "FINAL" if not steps else "NEEDS_ADDITIONAL_WORK"
    return {
        "verdict": v,
        "edge": {c["side"]: c for c in edge_classifications},
        "label": label,
        "status": status,
        "recommended_next_steps": steps,
    }


if __name__ == "__main__":
    main()
