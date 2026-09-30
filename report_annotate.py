#!/usr/bin/env python3
# Ported into SINEderella from SINE_discriminator's annotate.py (2026-09-09),
# import repointed to report_profile.py, plus the same consensus-row-detection
# fix as report_profile.profile() (SINEderella's own headers carry no
# "CONSENSUS_" substring -- use consensus_index()'s 0-fallback instead).
# Logic otherwise unmodified.
"""Structural feature annotation, to be drawn as a track under the profile.

Two kinds of feature, and the distinction matters when reading them:

  CONSENSUS features are properties of the query consensus, not of the copy set.
  Two sets searched with the same consensus get identical A/B boxes. They are
  descriptive, and they cannot discriminate between such sets - which is
  consistent with the spec treating them as one-sided evidence.

  COPY features are properties of the loci: TSDs, the poly-A actually present in
  each copy's 3' flank. These do vary between sets and do carry evidence.

Detected offline and reliably:
  abox, bbox      Pol III type-2 internal promoter, degenerate motif match
  trna_region     A-box start to B-box end, the tRNA-related head. This is an
                  operational definition from the boxes, NOT a database match to
                  a real tRNA - see requires-reference below.
  simple_repeat   homopolymer / di- / tri-nucleotide runs, in the consensus and
                  in the copies' 3' flank (the poly-A tail and any microsatellite)
  internal_dup    near-identical segment pairs within the consensus
  conserved_core  the most conserved block, empirically. NOT the CORE-SINE
                  domain, which needs its reference consensus to assert.
  tsd             from the per-copy detection in measure_c

Requires reference data, not attempted here, marked as such in the output:
  tRNA / 5S rRNA identity   needs a tRNA database (tRNAscan-SE, or ssearch36
                            against a tRNA set)
  LINE 3' end               needs LINE sequences; candidates exist on KIT at
                            .../teu/line/LINE_candidates.fa
  CORE domain               needs the CORE-SINE reference consensus
"""
import json, os, glob, sys, re
import numpy as np
import report_profile as M

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
         "W": "[AT]", "S": "[GC]", "K": "[GT]", "M": "[AC]", "N": "[ACGT]",
         "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]"}

# Pol III type-2 internal promoter, the standard degenerate descriptions
ABOX = "TRGCNNARYGG"
BBOX = "GWTCRANNC"


def to_re(m):
    return "".join(IUPAC[c] for c in m)


def best_motif(seq, motif, lo, hi, name):
    """Best partial match to a degenerate motif in seq[lo:hi]. Reports the
    number of mismatches, because an exact regex hit is rare in a diverged
    consensus and a silent miss would be worse than a scored near-match."""
    pat = [IUPAC[c] for c in motif]
    L = len(motif)
    best = None
    for i in range(max(0, lo), min(len(seq) - L + 1, hi)):
        mm = sum(0 if re.match(pat[j], seq[i + j]) else 1 for j in range(L))
        if best is None or mm < best[1]:
            best = (i, mm)
    if best is None:
        return None
    i, mm = best
    return {"type": name, "start": i, "end": i + L - 1, "mismatches": mm,
            "seq": seq[i:i + L], "max_mm": L}


def simple_repeats(seq, min_units=4, min_len=8, offset=0, tag="simple_repeat"):
    """Homopolymer, di- and tri-nucleotide runs."""
    out = []
    for period in (1, 2, 3):
        i = 0
        while i < len(seq) - period:
            unit = seq[i:i + period]
            if "N" in unit:
                i += 1
                continue
            j = i + period
            while j + period <= len(seq) and seq[j:j + period] == unit:
                j += period
            units = (j - i) // period
            if units >= min_units and (j - i) >= min_len:
                out.append({"type": tag, "start": i + offset, "end": j - 1 + offset,
                            "unit": unit, "units": units})
                i = j
            else:
                i += 1
    # drop a run fully contained in a longer one
    out.sort(key=lambda d: (d["start"], -(d["end"] - d["start"])))
    keep = []
    for d in out:
        if not any(k["start"] <= d["start"] and d["end"] <= k["end"] for k in keep):
            keep.append(d)
    return keep


def internal_dups(seq, k=12, min_len=16, max_mm_frac=0.15):
    """Near-identical segment pairs inside the consensus: seed on shared k-mers,
    extend along the diagonal, keep pairs that stay above the identity floor."""
    pos = {}
    for i in range(len(seq) - k + 1):
        pos.setdefault(seq[i:i + k], []).append(i)
    seen, out = set(), []
    for kmer, ps in pos.items():
        if len(ps) < 2:
            continue
        for a in range(len(ps)):
            for b in range(a + 1, len(ps)):
                i, j = ps[a], ps[b]
                d = j - i
                if d < k or (i, d) in seen:
                    continue
                s, e = i, i + k
                while e < len(seq) - d and seq[e] == seq[e + d]:
                    e += 1
                while s > 0 and seq[s - 1] == seq[s - 1 + d]:
                    s -= 1
                if e - s >= min_len:
                    for x in range(s, e):
                        seen.add((x, d))
                    out.append({"type": "internal_dup", "start": s, "end": e - 1,
                                "partner_start": s + d, "partner_end": e - 1 + d,
                                "len": e - s})
    out.sort(key=lambda d: -d["len"])
    return out[:6]


def terminators(seq, offset=0, tag="terminator"):
    """RNA Pol III terminates at a run of T on the non-template strand.

    strong    TTTT or longer - the canonical signal
    moderate  TCTTT / TGTTT / TATTT - a T-run interrupted by one non-T, which
              still terminates but less efficiently, so read-through is more
              likely and the element may acquire a 3' extension
    """
    out, i = [], 0
    while i < len(seq) - 3:
        if seq[i:i + 4] == "TTTT":
            j = i + 4
            while j < len(seq) and seq[j] == "T":
                j += 1
            out.append({"type": tag, "class": "strong", "start": i + offset,
                        "end": j - 1 + offset, "seq": seq[i:j], "run": j - i})
            i = j
            continue
        i += 1
    taken = set()
    for d in out:
        taken.update(range(d["start"] - offset, d["end"] - offset + 1))
    for i in range(len(seq) - 4):
        w = seq[i:i + 5]
        if w[0] == "T" and w[1] in "CGA" and w[2:] == "TTT" and not (set(range(i, i + 5)) & taken):
            out.append({"type": tag, "class": "moderate", "start": i + offset,
                        "end": i + 4 + offset, "seq": w})
            taken.update(range(i, i + 5))
    return sorted(out, key=lambda d: d["start"])


def annotate(path, prof=None):
    names, A = M.read_aln(path)
    k = M.consensus_index(names)  # see report_profile.profile() for why: this repo's
    # own headers are plain `>{species}_{subfamily}`, no CONSENSUS_ substring
    cons = A[k]
    nzc = np.where(cons != M.GAP)[0]
    cseq = "".join("ACGT"[c] for c in cons[nzc])
    L = len(cseq)
    C = np.delete(A, k, axis=0)
    n = C.shape[0]
    lo, hi = int(nzc[0]), int(nzc[-1])

    feats = []

    # --- consensus features -------------------------------------------------
    a = best_motif(cseq, ABOX, 0, min(60, L), "abox")
    b = best_motif(cseq, BBOX, 30, min(140, L), "bbox")
    for m in (a, b):
        if m and m["mismatches"] <= 4:
            m["kind"] = "consensus"
            feats.append(m)
    if a and b and a["mismatches"] <= 4 and b["mismatches"] <= 4 and b["start"] > a["end"]:
        feats.append({"type": "trna_region", "kind": "consensus",
                      "start": a["start"], "end": b["end"],
                      "note": "A-box to B-box; operational, not a tRNA database match"})

    for d in simple_repeats(cseq):
        d["kind"] = "consensus"
        feats.append(d)
    for d in terminators(cseq):
        d["kind"] = "consensus"
        feats.append(d)
    for d in internal_dups(cseq):
        d["kind"] = "consensus"
        feats.append(d)

    # most conserved block, empirical - not the CORE domain
    if prof and prof.get("cons_id"):
        xs = np.array(prof["x"], float)
        ci_ = np.array([np.nan if v is None else v for v in prof["cons_id"]], float)
        m = (xs >= 0) & (xs < L) & np.isfinite(ci_)
        if m.sum() > 40:
            y = ci_[m]
            w = 30
            ker = np.ones(w) / w
            sm = np.convolve(y, ker, mode="valid")
            s = int(np.argmax(sm))
            feats.append({"type": "conserved_core", "kind": "consensus",
                          "start": s, "end": s + w - 1,
                          "note": "most conserved 30 bp; empirical, not the CORE-SINE domain"})

    # --- copy features ------------------------------------------------------
    lefts, rights = [], []
    for i in range(n):
        l = C[i, :lo]
        lefts.append(l[l != M.GAP][::-1])
        r = C[i, hi + 1:]
        rights.append(r[r != M.GAP])

    # per-copy terminator: the functional one is the first T-run at or just past
    # the element 3' end, so search each copy's own downstream sequence
    strong = moderate = none_ = 0
    dists = []
    for r in rights:
        sq = "".join("ACGT"[c] for c in r[:40])
        t = terminators(sq)
        if not t:
            none_ += 1
            continue
        first = t[0]
        dists.append(first["start"])
        if first["class"] == "strong":
            strong += 1
        else:
            moderate += 1
    tot = max(1, len(rights))
    feats.append({"type": "terminator_copies", "kind": "copies",
                  "strong_frac": round(strong / tot, 3),
                  "moderate_frac": round(moderate / tot, 3),
                  "none_frac": round(none_ / tot, 3),
                  "dist_med": int(np.median(dists)) if dists else None,
                  "note": "first Pol III terminator within 40 bp downstream of each copy"})

    tsd = [M.find_tsd(lefts[i][::-1], rights[i]) for i in range(min(n, 120))]
    hits = [t for t in tsd if t > 0]
    if hits:
        feats.append({"type": "tsd", "kind": "copies",
                      "frac": round(len(hits) / len(tsd), 3),
                      "len_med": int(np.median(hits)),
                      "note": "direct repeat flanking both edges, per copy"})

    # poly-A / microsatellite actually present in the copies' 3' flank:
    # majority base per offset, then run-detect on that
    if rights:
        maj = []
        for off in range(1, 61):
            col = np.array([r[off - 1] if len(r) >= off else M.GAP for r in rights],
                           dtype=np.int8)
            bb = col[col != M.GAP]
            maj.append("ACGT"[int(np.bincount(bb, minlength=4)[:4].argmax())] if len(bb) >= 6 else "N")
        for d in simple_repeats("".join(maj), min_units=4, min_len=8,
                                offset=L, tag="tail_repeat"):
            d["kind"] = "copies"
            d["note"] = "in the copies' 3' flank, majority base per offset"
            feats.append(d)

    return {"set": os.path.basename(path).replace(".aln.fa", ""),
            "elem_len": L, "features": feats,
            "not_attempted": ["tRNA/5S database identity", "LINE 3' end",
                              "CORE-SINE domain"]}


def main():
    prof = {}
    if os.path.exists("profiles.json"):
        prof = json.load(open("profiles.json"))
    sets = sys.argv[1:] or list(prof.keys())
    out = {}
    for s in sets:
        p = os.path.join("aln_c", s + ".aln.fa")
        if not os.path.exists(p):
            continue
        a = annotate(p, prof.get(s))
        if a:
            out[s] = a
    json.dump(out, open("annotations.json", "w"), separators=(",", ":"))
    print("annotated %d sets" % len(out))
    for s in list(out)[:3]:
        print(" ", s, [f["type"] for f in out[s]["features"]])


if __name__ == "__main__":
    main()
