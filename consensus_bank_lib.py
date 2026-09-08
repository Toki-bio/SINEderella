"""Shared helpers for consensus bank audit, canonicalize, rebuild."""
from __future__ import annotations

import glob
import re
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

RC = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

# Defaults; override via env in callers if needed
TAIL_WINDOW = 60
TAIL_MAX_PERIOD = 6
TAIL_MIN_REPEATS = 3
ORIENT_MIN_MARGIN = 0.15
START_SLOP = 5  # allow tandem repeat to start a few bp in from 5'


def tandem_unit_weight(unit: str) -> float:
    """AT-rich units score as tail signal; G/C-only units are weak (body, not tail)."""
    if not unit:
        return 0.0
    at = sum(1 for c in unit if c in "AT") / len(unit)
    return 0.15 + 0.85 * at


def tandem_bp_score(bp: int, unit: str) -> float:
    return float(bp) * tandem_unit_weight(unit)


@dataclass
class OrientResult:
    seq: str
    action: str  # kept | flipped | undecided | empty
    fwd_score: float = 0.0
    rev_score: float = 0.0
    detail: dict[str, Any] = field(default_factory=dict)


def read_fa(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name, buf = None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(re.sub(r"\s", "", line))
    if name is not None:
        out[name] = "".join(buf).upper()
    return out


def write_fa(path: Path, seqs: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for name in sorted(seqs):
            seq = seqs[name]
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def rc(s: str) -> str:
    return s.translate(RC)[::-1]


def ungap(s: str) -> str:
    """ACGT-only sequence for identity / orientation (strip gaps and IUPAC)."""
    return re.sub(r"[^ACGT]", "", s.upper())


def identity(a: str, b: str) -> float:
    a, b = ungap(a), ungap(b)
    if not a or not b:
        return 0.0
    L = min(len(a), len(b))
    matches = sum(1 for i in range(L) if a[i] == b[i])
    return 100.0 * matches / L


def longest_tandem_at_end(
    s: str,
    max_period: int = TAIL_MAX_PERIOD,
    min_repeats: int = TAIL_MIN_REPEATS,
) -> dict[str, Any]:
    """Longest exact simple tandem repeat anchored at the 3' end of s."""
    s = ungap(s)
    best: dict[str, Any] = {"score": 0.0, "period": 0, "repeats": 0, "unit": "", "bp": 0}
    if len(s) < max_period * min_repeats:
        return best
    for period in range(1, max_period + 1):
        unit = s[-period:]
        count = 0
        pos = len(s)
        while pos >= period:
            if s[pos - period : pos] == unit:
                count += 1
                pos -= period
            else:
                break
        if count >= min_repeats:
            bp = count * period
            score = tandem_bp_score(bp, unit)
            if score > best["score"]:
                best = {"score": score, "period": period, "repeats": count, "unit": unit, "bp": bp}
    return best


def longest_tandem_near_start(
    s: str,
    max_period: int = TAIL_MAX_PERIOD,
    min_repeats: int = TAIL_MIN_REPEATS,
    max_slop: int = START_SLOP,
) -> dict[str, Any]:
    """Longest AT-weighted tandem repeat starting within max_slop bp of 5'."""
    s = ungap(s)
    best: dict[str, Any] = {"score": 0.0, "period": 0, "repeats": 0, "unit": "", "bp": 0}
    if len(s) < max_period * min_repeats:
        return best
    slop = min(max_slop + 1, len(s))
    for offset in range(slop):
        sub = s[offset:]
        for period in range(1, max_period + 1):
            if len(sub) < period * min_repeats:
                continue
            unit = sub[:period]
            count = 0
            pos = 0
            while pos + period <= len(sub):
                if sub[pos : pos + period] == unit:
                    count += 1
                    pos += period
                else:
                    break
            if count >= min_repeats:
                bp = count * period
                score = tandem_bp_score(bp, unit)
                if score > best["score"]:
                    best = {
                        "score": score,
                        "period": period,
                        "repeats": count,
                        "unit": unit,
                        "bp": bp,
                        "offset": offset,
                    }
    return best


def longest_tandem_at_start(
    s: str,
    max_period: int = TAIL_MAX_PERIOD,
    min_repeats: int = TAIL_MIN_REPEATS,
) -> dict[str, Any]:
    """Longest exact simple tandem repeat anchored at the 5' end of s."""
    s = ungap(s)
    best: dict[str, Any] = {"score": 0.0, "period": 0, "repeats": 0, "unit": "", "bp": 0}
    if len(s) < max_period * min_repeats:
        return best
    for period in range(1, max_period + 1):
        unit = s[:period]
        count = 0
        pos = 0
        while pos + period <= len(s):
            if s[pos : pos + period] == unit:
                count += 1
                pos += period
            else:
                break
        if count >= min_repeats:
            bp = count * period
            score = tandem_bp_score(bp, unit)
            if score > best["score"]:
                best = {"score": score, "period": period, "repeats": count, "unit": unit, "bp": bp}
    return best


def tail_strand_score(
    s: str,
    window: int = TAIL_WINDOW,
    max_period: int = TAIL_MAX_PERIOD,
    min_repeats: int = TAIL_MIN_REPEATS,
) -> tuple[float, dict[str, Any]]:
    """Higher when a simple tandem repeat sits at 3' rather than 5'."""
    u = ungap(s)
    if not u:
        return 0.0, {}
    tail3 = u[-window:] if len(u) > window else u
    tail5 = u[:window] if len(u) > window else u
    end3 = longest_tandem_at_end(tail3, max_period, min_repeats)
    end5 = longest_tandem_near_start(tail5, max_period, min_repeats)
    score = end3["score"] - 0.5 * end5["score"]
    return score, {"3prime": end3, "5prime": end5}


def orient_by_simple_repeat_tail(
    seq: str,
    window: int = TAIL_WINDOW,
    max_period: int = TAIL_MAX_PERIOD,
    min_repeats: int = TAIL_MIN_REPEATS,
    min_margin: float = ORIENT_MIN_MARGIN,
) -> OrientResult:
    u = ungap(seq)
    if not u:
        return OrientResult(seq="", action="empty")
    fwd, df = tail_strand_score(u, window, max_period, min_repeats)
    ru = rc(u)
    rev, dr = tail_strand_score(ru, window, max_period, min_repeats)
    scale = max(abs(fwd), abs(rev), 1.0)
    detail = {"fwd": df, "rev": dr, "fwd_score": fwd, "rev_score": rev}
    if abs(fwd - rev) < min_margin * scale:
        return OrientResult(seq=u, action="undecided", fwd_score=fwd, rev_score=rev, detail=detail)
    if rev > fwd:
        return OrientResult(seq=ru, action="flipped", fwd_score=fwd, rev_score=rev, detail=detail)
    return OrientResult(seq=u, action="kept", fwd_score=fwd, rev_score=rev, detail=detail)


def at_rich_3prime_score(s: str, tail: int = 30) -> float:
    u = ungap(s)
    if not u:
        return 0.0
    t = u[-tail:]
    return sum(1 for c in t if c in "AT") / len(t)


def pick_canonical(a: str, b: str) -> tuple[str, str, str]:
    sa, sb = ungap(a), ungap(b)
    id_direct = identity(sa, sb)
    id_rc = identity(sa, rc(sb))
    oa = orient_by_simple_repeat_tail(sa)
    ob = orient_by_simple_repeat_tail(sb)
    fa, _ = tail_strand_score(oa.seq)
    fb, _ = tail_strand_score(ob.seq)
    if id_rc > id_direct + 1:
        if fa >= fb:
            return oa.seq, "+", (
                f"RC pair; keep A (tail score {fa:.1f} vs {fb:.1f})"
            )
        return ob.seq, "+", (
            f"RC pair; keep B (tail score {fb:.1f} vs {fa:.1f})"
        )
    if fa >= fb:
        return oa.seq, "+", f"direct; keep A (tail score {fa:.1f})"
    return ob.seq, "+", f"direct; keep B (tail score {fb:.1f})"


def find_rc_clusters(cons: dict[str, str], min_id: float) -> list[list[str]]:
    names = sorted(cons)
    parent = {n: n for n in names}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i, a in enumerate(names):
        sa = cons[a]
        for b in names[i + 1 :]:
            sb = cons[b]
            if max(identity(sa, sb), identity(sa, rc(sb))) >= min_id:
                union(a, b)
    groups: dict[str, list[str]] = {}
    for n in names:
        groups.setdefault(find(n), []).append(n)
    return [sorted(g) for g in groups.values() if len(g) > 1]


def majority_consensus(seqs: list[str], tie_to_n: bool = False) -> str:
    """Column-wise majority on ungapped ACGT only; skip columns with no data."""
    seqs = [re.sub(r"[^ACGT]", "", s.upper()) for s in seqs if s]
    if not seqs:
        return ""
    L = max(len(s) for s in seqs)
    out = []
    for i in range(L):
        counts: dict[str, int] = {}
        for s in seqs:
            if i >= len(s):
                continue
            c = s[i]
            counts[c] = counts.get(c, 0) + 1
        if not counts:
            continue
        best_n = max(counts.values())
        winners = [b for b, n in counts.items() if n == best_n]
        if len(winners) == 1:
            out.append(winners[0])
        elif tie_to_n:
            out.append("N")
        else:
            out.append(winners[0])
    return "".join(out)


def find_step2_out(run_root: Path) -> Path | None:
    cands = sorted(glob.glob(str(run_root / "step2" / "step2_output*")))
    return Path(cands[-1]) if cands else None


def rebuild_from_subfams(run_root: Path, max_seqs: int = 500) -> dict[str, str]:
    s2 = find_step2_out(run_root)
    if not s2:
        return {}
    subdir = s2 / "subfamilies"
    if not subdir.is_dir():
        return {}
    rebuilt: dict[str, str] = {}
    for fa in sorted(subdir.glob("*.fasta")):
        seqs = list(read_fa(fa).values())[:max_seqs]
        if seqs:
            rebuilt[fa.stem] = majority_consensus(seqs, tie_to_n=False)
    return rebuilt


def pctid_stats(run_root: Path, sf: str) -> dict | None:
    s2 = find_step2_out(run_root)
    if not s2:
        return None
    tsv = s2 / "plots" / f"{sf}_pctid.tsv"
    if not tsv.is_file():
        return None
    vals = []
    with open(tsv) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                vals.append(100.0 - float(parts[1]))
    if not vals:
        return None
    vals.sort()
    return {
        "n": len(vals),
        "mean": statistics.mean(vals),
        "median": statistics.median(vals),
    }


def orient_at_rich_3prime(seq: str) -> str:
    """Canonical orientation for consensus bank (simple tandem repeat tail at 3')."""
    r = orient_by_simple_repeat_tail(seq)
    if r.action == "empty":
        return ungap(seq)
    if r.action == "undecided":
        u = ungap(seq)
        if at_rich_3prime_score(u) < at_rich_3prime_score(rc(u)):
            return rc(u)
        return u
    return r.seq
