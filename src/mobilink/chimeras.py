"""Chimera flags without a reference (P2 step 2.5, D50): overlap breakpoints and adapters inside reads.

Overlap breakpoints (minimap2 -x ava-ont among the carrier reads):
  - overlaps with an alignment block >= min_block (2,000 bp) and nmatch / block >= min_density (0.25);
  - an overlap end is internal when both reads continue past it by > term (200 bp): the two reads diverge there;
  - internal ends within cluster (200 bp) on a read form a candidate breakpoint; its support is the number of distinct
    partner reads;
  - a candidate is vetoed when an overlap spans it by >= span (500 bp) on each side;
  - it is flagged when support >= min_support (2) and it lies >= min_edge (500 bp) from both read ends.
  Mock communities (7 public runs): 0.6-2.2% of clean reads flagged, 73-100% of two-species chimeras, 0.7-7.5% of true
  links lost; with min_support 1 and min_edge 200, 12-38% of clean reads were flagged (most false breakpoints sit
  200-400 bp from a read end). Reads without overlaps cannot be judged ("assessable" False).
Adapters inside reads (optional; needs edlib): the ligation adapter core (all Kit 14 kits) and the rapid primer, plus
  native barcodes with their flanks if a barcode table is given, found in 500 bp windows anywhere >= 150 bp from both
  read ends. Specific (0.07-0.2% of clean mock reads) but current basecallers split most such reads already.
"""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd

LSK_CORE = "CCTGTACTTCGTTCAGTTACGTATTGC"
RAD_PRIMER = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
NB_FRONT, NB_REAR = "AAGGTTAA", "CAGCACCT"
_COMP = str.maketrans("ACGTN", "TGCAN")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def overlap_breakpoints(
    paf: pd.DataFrame,
    min_block: int = 2000,
    min_density: float = 0.25,
    term: int = 200,
    cluster: int = 200,
    span: int = 500,
):
    """paf: mobilink.loci.read_paf columns (q, qlen, qs, qe, strand, t, tlen, ts, te, nmatch, block).
    Returns (set of assessable reads, breakpoints DataFrame read_id, x, support, edge)."""
    o = paf[(paf.block >= min_block) & (paf.nmatch >= min_density * paf.block) & (paf.q != paf.t)]
    plus = (o.strand == "+").to_numpy()
    fw = pd.DataFrame({"r": o.q.to_numpy(), "a": o.qs.to_numpy(), "b": o.qe.to_numpy(), "rl": o.qlen.to_numpy(),
                       "p": o.t.to_numpy(),
                       "p_left": np.where(plus, o.ts > term, o.te < o.tlen - term),
                       "p_right": np.where(plus, o.te < o.tlen - term, o.ts > term)})
    rv = pd.DataFrame({"r": o.t.to_numpy(), "a": o.ts.to_numpy(), "b": o.te.to_numpy(), "rl": o.tlen.to_numpy(),
                       "p": o.q.to_numpy(),
                       "p_left": np.where(plus, o.qs > term, o.qe < o.qlen - term),
                       "p_right": np.where(plus, o.qe < o.qlen - term, o.qs > term)})
    d = pd.concat([fw, rv], ignore_index=True)
    assessed = set(d.r)
    cols = ["read_id", "x", "support", "edge"]
    left = d[(d.a > term) & d.p_left][["r", "a", "p", "rl"]].rename(columns={"a": "x"})
    right = d[(d.b < d.rl - term) & d.p_right][["r", "b", "p", "rl"]].rename(columns={"b": "x"})
    cand = pd.concat([left, right], ignore_index=True).sort_values(["r", "x"]).reset_index(drop=True)
    if cand.empty:
        return assessed, pd.DataFrame(columns=cols)
    cand["c"] = ((cand.r != cand.r.shift()) | (cand.x - cand.x.shift() > cluster)).cumsum()
    cl = cand.groupby("c").agg(read_id=("r", "first"), x=("x", "median"), support=("p", "nunique"),
                               rl=("rl", "first")).reset_index()
    iv = d[["r", "a", "b"]].sort_values(["r", "a"]).reset_index(drop=True)
    iv["bmax"] = iv.groupby("r").b.cummax()
    q = cl[["c", "read_id", "x"]].assign(key=lambda z: (z.x - span).astype(np.int64)).sort_values("key")
    iv["key"] = iv.a.astype(np.int64)
    m = pd.merge_asof(q, iv[["r", "key", "bmax"]].sort_values("key"), on="key", left_by="read_id", right_by="r",
                      direction="backward")
    spanned = m.set_index("c").bmax.reindex(cl.c).to_numpy() >= cl.x.to_numpy() + span
    bp = cl.loc[~spanned].copy()
    bp["edge"] = np.minimum(bp.x, bp.rl - bp.x)
    return assessed, bp[cols].reset_index(drop=True)


def adapter_hits(reads, barcodes: dict | None = None, edge: int = 150, win: int = 500, step: int = 400):
    """reads: iterable of (read_id, sequence). Returns DataFrame read_id, pattern, start, end, ed for hits >= edge
    from both read ends. Needs the optional dependency edlib."""
    import edlib

    pats = [("lsk", LSK_CORE, 5), ("lsk", _rc(LSK_CORE), 5), ("rad", RAD_PRIMER, 10), ("rad", _rc(RAD_PRIMER), 10)]
    for name, seq in (barcodes or {}).items():
        q = NB_FRONT + seq + NB_REAR
        pats += [(f"nb:{name}", q, 6), (f"nb:{name}", _rc(q), 6)]
    rows = []
    for rid, s in reads:
        n = len(s)
        for name, q, k in pats:
            last = -10**9
            o = 0
            while o < n:
                seg = s[o:o + win]
                if len(seg) >= len(q) // 2:
                    r = edlib.align(q, seg, mode="HW", task="locations", k=k)
                    if r["editDistance"] >= 0:
                        a, b = r["locations"][0]
                        a, b = o + a, o + b + 1
                        if a >= edge and b <= n - edge and a - last > 50:
                            rows.append((rid, name, a, b, r["editDistance"]))
                            last = a
                if o + win >= n:
                    break
                o += step
    return pd.DataFrame(rows, columns=["read_id", "pattern", "start", "end", "ed"])


def read_sequences(path):
    """FASTA or FASTQ, plain or gzip: yields (read_id, sequence)."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        first = f.read(1)
        f.seek(0)
        if first == "@":
            for h in f:
                s = next(f).strip().upper()
                next(f)
                next(f)
                yield h[1:].split()[0], s
            return
        name, seq = None, []
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip().upper())
        if name is not None:
            yield name, "".join(seq)


def chimera_flags(paf=None, reads=None, barcodes=None, min_support: int = 2, min_edge: int = 500) -> pd.DataFrame:
    """One row per judged read: read_id, overlap_assessable, n_breakpoints, flag_overlap, n_adapters, flag_adapter,
    flag (either). paf: read_paf output or None; reads: iterable of (read_id, sequence) or None."""
    out = []
    if paf is not None:
        assessed, bp = overlap_breakpoints(paf)
        ok = bp[(bp.support >= min_support) & (bp.edge >= min_edge)]
        n = ok.groupby("read_id").size()
        out.append(pd.DataFrame({"read_id": sorted(assessed)}).assign(
            overlap_assessable=True, n_breakpoints=lambda z: z.read_id.map(n).fillna(0).astype(int)))
    if reads is not None:
        h = adapter_hits(reads, barcodes)
        out.append(h.groupby("read_id").size().rename("n_adapters").reset_index())
    if not out:
        raise ValueError("give a PAF, reads, or both")
    t = out[0]
    for x in out[1:]:
        t = t.merge(x, on="read_id", how="outer")
    t["overlap_assessable"] = t.get("overlap_assessable", pd.Series(False, index=t.index)).fillna(False).astype(bool)
    t["n_breakpoints"] = t.get("n_breakpoints", pd.Series(0, index=t.index)).fillna(0).astype(int)
    t["n_adapters"] = t.get("n_adapters", pd.Series(0, index=t.index)).fillna(0).astype(int)
    t["flag_overlap"] = t.n_breakpoints > 0
    t["flag_adapter"] = t.n_adapters > 0
    t["flag"] = t.flag_overlap | t.flag_adapter
    return t.sort_values("read_id").reset_index(drop=True)
