"""Shared helpers for consensus bank audit, canonicalize, rebuild."""
from __future__ import annotations

import glob
import re
import statistics
from pathlib import Path

RC = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


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
    if id_rc > id_direct + 1:
        score_a = at_rich_3prime_score(sa)
        score_b = at_rich_3prime_score(sb)
        if score_a >= score_b:
            return sa, "+", f"RC pair; keep A (AT3'={score_a:.2f} vs {score_b:.2f})"
        return rc(sb), "-", f"RC pair; flip B (AT3'={score_b:.2f} vs {score_a:.2f})"
    if at_rich_3prime_score(sa) >= at_rich_3prime_score(sb):
        return sa, "+", "direct; keep A by AT3'"
    return sb, "+", "direct; keep B by AT3'"


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
            out.append(winners[0])  # deterministic: ACGT order via max key
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
    u = ungap(seq)
    if at_rich_3prime_score(u) < at_rich_3prime_score(rc(u)):
        return rc(u)
    return u
