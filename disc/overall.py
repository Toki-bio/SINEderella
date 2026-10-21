#!/usr/bin/env python3
"""The Overall call on a published plate: the verdict score plus what the copies were seen to do.

verdict.py withholds its score below 30 copies (INSUFFICIENT_COPIES -> "Cannot assess"), so that a
small but clean family is not rejected for being small (HUM AluYb9). That rule stays. What it also
did was hide negative evidence the copies DO show: on rsi (2026-09-27) MEG-RL had 20 copies, no
supported element (NO_ELEMENT), a consensus the copies cut ~45 bp short at the 3' end, and left and
right flanks shared by 14-20 % of copies - and the report said only "Cannot assess". The flank
context column (flank_uniqueness.scan) was never part of the Overall call at all.

Rules:
  * assessable (>= 30 copies, flanks present): the verdict score as before.
  * not assessable, with negative evidence in the verdict flags: "Doubtful", reasons listed.
  * not assessable, nothing negative observed: "Cannot assess", as before.
Every call carries its reasons, shown in the tooltip.

The flank-context scan is listed ("also seen") but does not move the call. Tested 2026-09-27 on the
224 published plates scoring >= 75: 30 have rand100 "high" shared flanks, among them hum CAS, FLAM_C,
FRAM and Rhin-1 in six bats - random copies sharing a flank group is common in real families, so it
cannot demote a call on its own.
"""

NEGATIVE = {
    "NO_ELEMENT": "the copies do not support the consensus above background",
    "SMALL_CORE": "too few copies form the core the consensus is built on",
    "MICROSATELLITE_ELEMENT": "the element is largely simple repeat",
    "NOT_ISOLATED": "copies sit in satellites or duplications",
    "SHARED_FLANKS": "copies share flanking sequence",
    "FRAGMENT_OF_LONGER": "the element continues past the consensus (part of a longer repeat)",
    "ELEMENT_CONTINUES": "the element continues past the consensus (part of a longer repeat)",
    "CONSENSUS_OVEREXTENDED": "the consensus is longer than the element the copies support",
    "TANDEM_ARRAY": "many copies sit in tandem clusters: copies of a repeated unit, not insertions",
}


# Plain wording (his review 2026-09-28: "Flank context and Overall comments should be more standardized,
# more readable without jargon"). One fixed layout for every row; no internal codes on the page.
CAP_WORDS = {code: text for code, text in NEGATIVE.items()}


def side_line(side_name, side):
    """'left side: 12 of 97 copies share flanking DNA with another copy (largest group: 5 copies)'"""
    if not side or not side.get("measured"):
        why = (side or {}).get("reason", "no flank sequence")
        return "%s side: not measured (%s)" % (side_name, str(why).replace("_", " "))
    n = int(side.get("n_measured", side.get("n", 0)) or 0)
    k = int(round(side.get("shared_copy_frac", 0) * n))
    if k == 0:
        return "%s side: all %d copies have their own flanking DNA" % (side_name, n)
    return ("%s side: %d of %d copies share flanking DNA with another copy (largest group: %d copies)"
            % (side_name, k, n, side.get("largest_cluster", 0)))


def flank_context_text(top, rand):
    """(worst severity, tooltip lines) for the Flank context chip - same layout on every row."""
    lines = ["Are the copies independent insertions? Checked by comparing the DNA just outside each copy."]
    worst = None
    for label, r in (("Top 100", top), ("100 random", rand)):
        if not r:
            continue
        if r.get("error"):
            lines.append("%s: not measured (%s)." % (label, str(r["error"]).replace("_", " ")))
            continue
        lines.append("%s - %s; %s." % (label, side_line("left", r.get("left")), side_line("right", r.get("right"))))
        w = r.get("worst_flag")
        if w == "high" or (w == "medium" and worst is None):
            worst = w
    lines.append({
        None: "Reading: flanking DNA differs between copies, as expected for independent insertions.",
        "medium": "Reading: a small group of copies shares flanking DNA - common in real families (tandem "
                  "copies, duplicated regions); worth a look, it does not change the call.",
        "high": "Reading: many randomly chosen copies share flanking DNA, so some of these loci may not be "
                "independent insertions - look for tandem copies, duplicated regions, or copies inside a "
                "larger repeat.",
    }[worst])
    return worst, lines


def _flags(v):
    return {f["code"]: f for f in v.get("flags", [])}


def _context(top, rand):
    """(severity, short plain lines) for the flagged sides of the top100 / rand100 flank scans."""
    sev, out = None, []
    for label, r in (("100 random", rand), ("top 100", top)):
        if not r or r.get("error"):
            continue
        for f in r.get("flags", []):
            out.append("%s, %s" % (label, side_line(f["side"], r.get(f["side"]))))
        w = r.get("worst_flag")
        if w == "high" or (w == "medium" and sev is None):
            sev = w
    return sev, out


def overall_call(v, top_scan=None, rand_scan=None):
    """-> {"label", "kind" (ok/edge/warn/muted), "reasons": [str]}"""
    if not v or v.get("error"):
        return {"label": "n/a", "kind": "muted",
                "reasons": [v.get("error", "not scored") if v else "not scored"]}
    if v.get("deferred"):
        return {"label": "Deferred", "kind": "edge",
                "reasons": ["Mixture of two groups - split before a firm call."]}

    fl = _flags(v)
    n = v.get("n", 0)
    neg = [NEGATIVE[c] + (" (%d/%d copies support it)" % (v.get("n_supported", 0), n)
                          if c == "NO_ELEMENT" else "")
           for c in NEGATIVE if c in fl]
    # the two ELEMENT_CONTINUES/FRAGMENT texts are identical; keep one
    neg = list(dict.fromkeys(neg))
    sev, ctx = _context(top_scan, rand_scan)
    also = (["Also noted (does not change the call): some copies share flanking DNA -"] + ["- " + c for c in ctx]
            if ctx else [])

    if not v.get("assessable", True):
        why = []
        if "INSUFFICIENT_COPIES" in fl:
            why.append("only %d copies (a score needs about 30)" % n)
        if "NO_FLANKS_PRESENT" in fl:
            why.append("no flanking DNA to measure")
        head = "No score: not enough evidence - " + "; ".join(why) + "."
        if neg:
            return {"label": "Doubtful", "kind": "warn",
                    "reasons": [head, "Evidence against:"]
                    + ["- " + r for r in neg] + also}
        return {"label": "Cannot assess", "kind": "muted",
                "reasons": [head, "Nothing against it was seen."] + also}

    s = float(v.get("score", 0))
    if s >= 90:
        label, kind = "SINE", "ok"
    elif s >= 75:
        label, kind = "SINE (caveats)", "ok"
    elif s >= 55:
        label, kind = "Grey zone", "edge"
    else:
        label, kind = "Not SINE", "warn"
    reasons = ["Score %.0f of 100: %s." % (s, label)]
    if v.get("capped_by"):
        reasons.append("The score is held down because: %s."
                       % "; ".join(CAP_WORDS.get(c, c.lower().replace("_", " ")) for c in v["capped_by"]))
    shown = set(CAP_WORDS.get(c) for c in v.get("capped_by", []))
    rest = [r for r in neg if r.split(" (")[0] not in shown]
    if rest:
        reasons += ["Evidence against:"] + ["- " + r for r in rest]
    reasons += also
    return {"label": label, "kind": kind, "reasons": reasons}
