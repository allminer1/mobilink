"""Unit loci: the A units on different reads that are copies of the same locus (P2 step 2.7, probe v4).

Reads that see the same side of the same locus share one distance to the nearest B, so under no linkage their
outcomes are one random draw: the null-variance score test (mobilink.stats.nvtest) needs these clumps.

Real data (unit_loci_from_paf): minimap2 -x ava-ont overlaps among carrier reads (read_paf, blocks >= 1,000 bp),
kept if nmatch / block >= min_density (homologous genes at different loci align with low match density; P2 step
2.3: median 0.92 within a locus, 0.15 between loci). Each A unit on the query that lies >= 80% inside the overlap is
mapped to the target linearly (scaled by the two overlap lengths; reversed for strand "-"); it joins a target unit
that overlaps the mapped interval by >= min_frac of the shorter. Joins carry a parity (1 for strand "-"): the right
side on one read is the left side on the other. A locus whose joins disagree on parity (inverted repeats, chimeras)
is flagged and its two sides are merged into one clump, which can only overstate the variance.
Simulation (unit_loci_from_positions): genome intervals of the units (reads are forward), joined if they overlap
by >= min_frac of the shorter.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def _find(parent, par, x):
    """Root of x and the parity of x relative to it (path compression)."""
    path = []
    while parent[x] != x:
        path.append(x)
        x = parent[x]
    root, acc = x, 0
    for y in reversed(path):  # nearest to root first
        acc ^= par[y]
        par[y] = acc
        parent[y] = root
    return root


def parity_components(n: int, a: np.ndarray, b: np.ndarray, p: np.ndarray):
    """Union-find with parity over n nodes and edges (a, b, parity). Returns (component, parity to root,
    conflict per node: its component has an edge whose parity disagrees)."""
    parent, par = np.arange(n), np.zeros(n, np.int8)
    bad_root = np.zeros(n, bool)
    for x, y, q in zip(a.tolist(), b.tolist(), p.tolist(), strict=True):
        rx, ry = _find(parent, par, x), _find(parent, par, y)
        px, py = par[x] if x != rx else 0, par[y] if y != ry else 0
        if rx == ry:
            if (px ^ py) != q:
                bad_root[rx] = True
            continue
        parent[ry] = rx
        par[ry] = px ^ py ^ q
        bad_root[rx] |= bad_root[ry]
    root = np.array([_find(parent, par, x) for x in range(n)])
    rel = np.where(root == np.arange(n), 0, par).astype(np.int8)
    return pd.factorize(root)[0], rel, bad_root[root]


def unit_loci_from_paf(
    units: pd.DataFrame, paf: pd.DataFrame, min_frac: float = 0.5, min_density: float = 0.25, min_inside: float = 0.8
):
    """units: molecule_id, start, end (one row per A unit occurrence); paf: mobilink.loci.read_paf output.
    Returns a DataFrame aligned with units: locus, parity (0/1 relative to the locus), conflict (bool)."""
    u = units[["molecule_id", "start", "end"]].reset_index(drop=True).assign(u=lambda x: np.arange(len(x)))
    mols = set(u.molecule_id)
    p = paf[(paf.nmatch >= min_density * paf.block) & paf.q.isin(mols) & paf.t.isin(mols)]
    p = p[["q", "qs", "qe", "strand", "t", "ts", "te"]].reset_index(drop=True)
    a_list, b_list, par_list = [], [], []
    step = 2_000_000
    for i in range(0, len(p), step):
        c = p.iloc[i : i + step].merge(u, left_on="q", right_on="molecule_id")
        inside = np.minimum(c.end, c.qe) - np.maximum(c.start, c.qs)
        c = c[inside >= min_inside * (c.end - c.start)]
        scale = (c.te - c.ts) / (c.qe - c.qs).clip(lower=1)
        x0, x1 = (c.start - c.qs) * scale, (c.end - c.qs) * scale
        rev = (c.strand == "-").to_numpy()
        c = c.assign(ms=np.where(rev, c.te - x1, c.ts + x0), me=np.where(rev, c.te - x0, c.ts + x1), rev=rev)
        c = c[["t", "ms", "me", "u", "rev"]].merge(u.rename(columns={"u": "v"}), left_on="t", right_on="molecule_id")
        ov = np.minimum(c.me, c.end) - np.maximum(c.ms, c.start)
        short = np.minimum(c.me - c.ms, c.end - c.start)
        c = c[ov >= min_frac * short]
        a_list.append(c.u.to_numpy()), b_list.append(c.v.to_numpy()), par_list.append(c.rev.to_numpy(np.int8))
    a = np.concatenate(a_list) if a_list else np.array([], int)
    b = np.concatenate(b_list) if b_list else np.array([], int)
    q = np.concatenate(par_list) if par_list else np.array([], np.int8)
    comp, rel, bad = parity_components(len(u), a, b, q)
    return pd.DataFrame({"locus": comp, "parity": rel, "conflict": bad})


def unit_loci_from_positions(replicon, gstart, gend, min_frac: float = 0.5) -> np.ndarray:
    """Simulation: connected components of unit genome intervals on one replicon overlapping by >= min_frac of
    the shorter (reads are forward, so the parity is 0)."""
    rep, gs, ge = np.asarray(replicon), np.asarray(gstart), np.asarray(gend)
    n = rep.size
    o = np.lexsort((gs, rep))
    rows, cols = [], []
    for k in range(n):
        i = o[k]
        j = k + 1
        while j < n and rep[o[j]] == rep[i] and gs[o[j]] < ge[i]:
            m = o[j]
            ov = min(ge[i], ge[m]) - max(gs[i], gs[m])
            if ov >= min_frac * min(ge[i] - gs[i], ge[m] - gs[m]):
                rows.append(i), cols.append(m)
            j += 1
    g = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    return connected_components(g, directed=False)[1]
