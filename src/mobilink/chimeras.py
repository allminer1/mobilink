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
  A PAF given as a path is read in one pass and never held in memory (overlap_breakpoints_stream): each overlap gives
  one small record per read, the records are spilled to bucket files by read, and the breakpoints are found bucket by
  bucket. The result is the same as with the table in memory.
Adapters inside reads (optional; needs edlib): the ligation adapter core (all Kit 14 kits) and the rapid primer, plus
  native barcodes with their flanks if a barcode table is given, found in 500 bp windows anywhere >= 150 bp from both
  read ends. Specific (0.07-0.2% of clean mock reads) but current basecallers split most such reads already.
"""

from __future__ import annotations

import gzip
import os
import tempfile

import numpy as np
import pandas as pd

LSK_CORE = "CCTGTACTTCGTTCAGTTACGTATTGC"
RAD_PRIMER = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
NB_FRONT, NB_REAR = "AAGGTTAA", "CAGCACCT"
_COMP = str.maketrans("ACGTN", "TGCAN")
_BP_COLS = ["read_id", "x", "support", "edge"]
# one record per overlap and read: read, overlap start and end on it, its length, partner, and whether the partner
# continues past the left and the right end of the overlap
_REC = np.dtype([("r", "<i4"), ("a", "<i4"), ("b", "<i4"), ("rl", "<i4"), ("p", "<i4"), ("pl", "?"), ("pr", "?")])


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
    d = _side_records(_usable(paf, min_block, min_density), term)
    return set(d.r), _breakpoints(d, term, cluster, span)


def _usable(paf: pd.DataFrame, min_block: int, min_density: float) -> pd.DataFrame:
    return paf[(paf.block >= min_block) & (paf.nmatch >= min_density * paf.block) & (paf.q != paf.t)]


def _side_records(o: pd.DataFrame, term: int) -> pd.DataFrame:
    """Each usable overlap seen from both of its reads: r (the read), a, b (overlap on it), rl (its length), p (the
    partner), p_left / p_right (the partner continues past that end of the overlap by more than term)."""
    plus = (o.strand == "+").to_numpy()
    fw = pd.DataFrame({"r": o.q.to_numpy(), "a": o.qs.to_numpy(), "b": o.qe.to_numpy(), "rl": o.qlen.to_numpy(),
                       "p": o.t.to_numpy(),
                       "p_left": np.where(plus, o.ts > term, o.te < o.tlen - term),
                       "p_right": np.where(plus, o.te < o.tlen - term, o.ts > term)})
    rv = pd.DataFrame({"r": o.t.to_numpy(), "a": o.ts.to_numpy(), "b": o.te.to_numpy(), "rl": o.tlen.to_numpy(),
                       "p": o.q.to_numpy(),
                       "p_left": np.where(plus, o.qs > term, o.qe < o.qlen - term),
                       "p_right": np.where(plus, o.qe < o.qlen - term, o.qs > term)})
    return pd.concat([fw, rv], ignore_index=True)


def _breakpoints(d: pd.DataFrame, term: int, cluster: int, span: int) -> pd.DataFrame:
    """Breakpoints from side records. Every step is within one read, so the records may be processed in any
    partition that keeps the records of a read together."""
    cols = _BP_COLS
    left = d[(d.a > term) & d.p_left][["r", "a", "p", "rl"]].rename(columns={"a": "x"})
    right = d[(d.b < d.rl - term) & d.p_right][["r", "b", "p", "rl"]].rename(columns={"b": "x"})
    cand = pd.concat([left, right], ignore_index=True).sort_values(["r", "x"]).reset_index(drop=True)
    if cand.empty:
        return pd.DataFrame(columns=cols)
    cand["c"] = ((cand.r != cand.r.shift()) | (cand.x - cand.x.shift() > cluster)).cumsum()
    cl = cand.groupby("c").agg(read_id=("r", "first"), x=("x", "median"), support=("p", "nunique"),
                               rl=("rl", "first")).reset_index()
    iv = d[["r", "a", "b"]].sort_values(["r", "a"]).reset_index(drop=True)
    iv["bmax"] = iv.groupby("r").b.cummax()
    # overlaps of a read that start at the same position: keep the row with the furthest end reached so far, so that
    # the veto below does not depend on the order of the rows (before v1.0 it could miss a spanning overlap in a tie)
    iv = iv.drop_duplicates(["r", "a"], keep="last")
    q = cl[["c", "read_id", "x"]].assign(key=lambda z: (z.x - span).astype(np.int64)).sort_values("key")
    iv["key"] = iv.a.astype(np.int64)
    m = pd.merge_asof(q, iv[["r", "key", "bmax"]].sort_values("key"), on="key", left_by="read_id", right_by="r",
                      direction="backward")
    spanned = m.set_index("c").bmax.reindex(cl.c).to_numpy() >= cl.x.to_numpy() + span
    bp = cl.loc[~spanned].copy()
    bp["edge"] = np.minimum(bp.x, bp.rl - bp.x)
    return bp[cols].reset_index(drop=True)


def overlap_breakpoints_stream(
    path,
    min_block: int = 2000,
    min_density: float = 0.25,
    term: int = 200,
    cluster: int = 200,
    span: int = 500,
    chunksize: int = 2_000_000,
    buckets: int = 64,
    tmp_dir=None,
):
    """overlap_breakpoints for a PAF file (gzip ok) that is too large to hold in memory; same result.
    The file is read once in chunks. Read names become integer codes, each usable overlap gives two 22-byte side
    records, and the records go to `buckets` temporary files by read code (in tmp_dir, default: the system's
    temporary directory; removed at the end). Memory: one chunk, the read names, and the records of one bucket.
    Disk: 44 bytes per usable overlap."""
    from .loci import PAF_COLS

    names: dict[str, int] = {}
    parts = []
    with tempfile.TemporaryDirectory(dir=tmp_dir, prefix="mobilink_chimeras_") as tmp:
        files = [open(os.path.join(tmp, f"{k}.bin"), "wb") for k in range(buckets)]
        try:
            try:
                reader = pd.read_csv(path, sep="\t", header=None, usecols=range(11), names=PAF_COLS[:11],
                                     dtype={"q": str, "t": str, "strand": str}, chunksize=chunksize)
            except pd.errors.EmptyDataError:  # no overlaps at all
                reader = ()
            for p in reader:
                p = _usable(p, min_block, min_density)
                if p.empty:
                    continue
                local, uniq = pd.factorize(pd.concat([p.q, p.t], ignore_index=True))
                code = np.fromiter((names.setdefault(x, len(names)) for x in uniq), np.int64, len(uniq))[local]
                d = _side_records(p.assign(q=code[: len(p)], t=code[len(p) :]), term)
                rec = np.empty(len(d), _REC)
                for f, c in zip(_REC.names, ("r", "a", "b", "rl", "p", "p_left", "p_right"), strict=True):
                    rec[f] = d[c].to_numpy()
                k = rec["r"] % buckets
                order = np.argsort(k, kind="stable")
                rec, edges = rec[order], np.searchsorted(k[order], np.arange(buckets + 1))
                for j in range(buckets):
                    if edges[j + 1] > edges[j]:
                        files[j].write(rec[edges[j] : edges[j + 1]].tobytes())
        finally:
            for fh in files:
                fh.close()
        for j in range(buckets):
            rec = np.fromfile(os.path.join(tmp, f"{j}.bin"), _REC)
            if rec.size:
                d = pd.DataFrame({"r": rec["r"], "a": rec["a"], "b": rec["b"], "rl": rec["rl"], "p": rec["p"],
                                  "p_left": rec["pl"], "p_right": rec["pr"]})
                parts.append(_breakpoints(d, term, cluster, span))
    ids = np.array(list(names), dtype=object)  # insertion order = code order
    parts = [x for x in parts if len(x)]
    if not parts:
        return set(names), pd.DataFrame(columns=_BP_COLS)
    bp = pd.concat(parts, ignore_index=True)
    bp["read_id"] = ids[bp.read_id.to_numpy(np.int64)]
    return set(names), bp.sort_values(["read_id", "x"]).reset_index(drop=True)


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


def chimera_flags(paf=None, reads=None, barcodes=None, min_support: int = 2, min_edge: int = 500,
                  tmp_dir=None) -> pd.DataFrame:
    """One row per judged read: read_id, overlap_assessable, n_breakpoints, flag_overlap, n_adapters, flag_adapter,
    flag (either). paf: a PAF file (path; read in one pass, see overlap_breakpoints_stream), a read_paf table, or
    None; reads: iterable of (read_id, sequence) or None; tmp_dir: where the temporary files of a PAF path go."""
    out = []
    if paf is not None:
        if isinstance(paf, (str, os.PathLike)):
            assessed, bp = overlap_breakpoints_stream(paf, tmp_dir=tmp_dir)
        else:
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
