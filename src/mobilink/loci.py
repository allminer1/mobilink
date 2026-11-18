"""Locus clusters: groups of reads that cover the same locus (P2 step 2.3; amendment 04 draft, Part 2).

Real data: minimap2 -x ava-ont overlaps among the reads that carry a call (PAF). Two reads join if the overlap's
alignment block is >= min_block bp and >= min_free bp of the overlap lie outside MGE calls on both reads. The second
rule stops copies of one IS element or transposon at different places from merging unrelated loci: two insertion
sites share the element but not its flanks, so their overlap lies (almost) entirely inside the MGE call.
Clusters are the connected components; a read without accepted overlaps is its own cluster.
Simulation: reads on the same replicon whose intervals overlap by >= min_block bp (exact, from positions).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

PAF_COLS = ["q", "qlen", "qs", "qe", "strand", "t", "tlen", "ts", "te", "nmatch", "block", "mapq"]


def read_paf(path, min_block: int = 2000, keep=None, chunksize: int = 5_000_000) -> pd.DataFrame:
    """The first 12 PAF columns; self-hits and blocks shorter than min_block dropped. keep: a set of read ids; if given,
    only overlaps with both reads in it are kept (e.g. the reads that carry a class-A call), read in chunks so that a
    large PAF does not have to fit in memory."""
    parts = []
    try:
        reader = pd.read_csv(path, sep="\t", header=None, usecols=range(12), names=PAF_COLS, chunksize=chunksize)
    except pd.errors.EmptyDataError:  # no overlaps at all
        return pd.DataFrame(columns=PAF_COLS)
    for p in reader:
        p = p[(p.q != p.t) & (p.block >= min_block)]
        if keep is not None:
            p = p[p.q.isin(keep) & p.t.isin(keep)]
        parts.append(p)
    if not parts:
        return pd.DataFrame(columns=PAF_COLS)
    return pd.concat(parts, ignore_index=True)


def _inside(read, s, e, calls: pd.DataFrame) -> np.ndarray:
    """Bases of [s, e) on `read` covered by the calls (molecule_id, start, end), per row (calls assumed disjoint)."""
    q = pd.DataFrame({"m": read, "s": s, "e": e, "i": np.arange(len(read))})
    c = calls[["molecule_id", "start", "end"]].rename(columns={"molecule_id": "m"})
    j = q.merge(c, on="m")
    ov = np.clip(np.minimum(j.e, j.end) - np.maximum(j.s, j.start), 0, None)
    out = np.zeros(len(read))
    np.add.at(out, j.i.to_numpy(), ov.to_numpy(float))
    return out


def accepted_overlaps(paf: pd.DataFrame, mge_calls: pd.DataFrame | None = None, min_free: int = 1000) -> pd.DataFrame:
    """PAF rows that pass the MGE guard; adds free_q and free_t (overlap bases outside MGE calls)."""
    p = paf.copy()
    p["free_q"] = (p.qe - p.qs).astype(float)
    p["free_t"] = (p.te - p.ts).astype(float)
    if mge_calls is not None and len(mge_calls):
        p["free_q"] -= _inside(p.q.to_numpy(), p.qs.to_numpy(), p.qe.to_numpy(), mge_calls)
        p["free_t"] -= _inside(p.t.to_numpy(), p.ts.to_numpy(), p.te.to_numpy(), mge_calls)
    return p[(p.free_q >= min_free) & (p.free_t >= min_free)].reset_index(drop=True)


def locus_clusters(molecule_ids, paf: pd.DataFrame, mge_calls: pd.DataFrame | None = None, min_free: int = 1000):
    """Cluster code per molecule (0..G-1). paf as read_paf (already filtered by block length)."""
    ids = pd.Index(molecule_ids)
    acc = accepted_overlaps(paf, mge_calls, min_free)
    i, j = ids.get_indexer(acc.q), ids.get_indexer(acc.t)
    ok = (i >= 0) & (j >= 0)
    n = len(ids)
    g = coo_matrix((np.ones(ok.sum()), (i[ok], j[ok])), shape=(n, n))
    _, lab = connected_components(g, directed=False)
    return lab


def clusters_from_positions(replicon, start, length, min_block: int = 2000) -> np.ndarray:
    """Exact loci for simulated reads: reads on one replicon overlapping by >= min_block join (connected components).
    Sweep by start with union-find: the earlier read with the largest end has the largest overlap with the current
    read, and any two earlier reads that both overlap it by >= min_block overlap each other by >= min_block (both
    contain [start, start + min_block]), so joining the current read to that one read is enough."""
    rep, s = np.asarray(replicon), np.asarray(start, dtype=np.int64)
    e = s + np.asarray(length, dtype=np.int64)
    o = np.lexsort((s, rep))
    parent = np.arange(rep.size)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    best, best_end, cur_rep = -1, -np.inf, None
    for k in o:
        if rep[k] != cur_rep:
            best, best_end, cur_rep = k, e[k], rep[k]
            continue
        if min(best_end, e[k]) - s[k] >= min_block:
            ra, rb = find(k), find(best)
            if ra != rb:
                parent[ra] = rb
        if e[k] > best_end:
            best, best_end = k, e[k]
    roots = np.array([find(x) for x in range(rep.size)])
    return pd.factorize(roots)[0].astype(np.int64)
