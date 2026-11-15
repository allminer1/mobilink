"""Null-variance score tests for the net linked test (P2 step 2.7, probes v3-v4; amendment 04 draft v2).

Why not the observed variance: per host or per locus cluster, O - E is skewed (rare large positive terms where a
locus sees a B, many small negative ones), so sum c / sqrt(sum c^2) has SD 1.1-1.9 and a negative mean under no
linkage (P2_07d-f). The variance is computed under the null instead.

Sides: for each A-unit side, the room to the molecule end on that side and the distance to the far end of the
nearest B on it (as linkstat.stats.netlinked.sides, with the side recorded). Chance: F_h, the weighted
Kaplan-Meier CDF of the host's anchor sides (anchors placed within host strata).
  net_w (window w): y_i = 1[nearest B within tau_i], tau_i = min(w, room_i); E = sum F_h(tau_i).
  mol (same molecule): per unit occurrence y = 1[a B on either side within the room];
       P(y = 0) = (1 - F_h(room_L)) (1 - F_h(room_R)).
Clumps: sides that see the same genome side of the same locus share one distance D ~ F_h (net_w clump = unit
locus x oriented side; mol clump = unit locus). Then
  net_w: Var_k = sum_ij F(min(tau_i, tau_j)) - E_k^2          (F increasing: F(min tau) = min F)
  mol:   Cov_ij = (1 - max F_L)(1 - max F_R) - q_i q_j.
Orientation: a side's genome side is its read side XOR the unit's parity (linkstat.unitloci). Loci whose overlaps
disagree on orientation (conflict) get one shared draw for all their outcomes (comonotone; it maximises the clump
variance for the given marginals, so it is conservative).
Reference (P2 step 2.7, run 6 -> D45): the exact null distribution of O under this clump model, by Monte Carlo
(n_mc draws; each clump draws its shared distances by inverse CDF, U ~ U(0, 1)), with mid-p values:
  p_up = (#{O* > O} + (#{O* = O} + 1) / 2) / (n_mc + 1), p_low likewise, p = min(1, 2 min(p_up, p_low)).
Why not the normal: at few effective events O is skewed and discrete; in run 6 (1,800 datasets under no linkage)
the normal net_w test rejected 4.7% on the enrichment side at 2.5% when E^2 / V < 5, and almost never on the
depletion side. z = (O - E) / sqrt(sum_k Var_k) and E^2 / Var stay in the output; below 10 events are flagged.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def sides_dir(lengths, rid, start, end, is_a, is_b, cut) -> pd.DataFrame:
    """Per A side: i (feature row), t, e (as netlinked.sides), side (0 right, 1 left), room (to the molecule end)."""
    ia, ib = np.flatnonzero(is_a), np.flatnonzero(is_b)
    ib = ib[np.argsort(rid[ib], kind="stable")]
    b_rid = rid[ib]
    out = []
    for i in ia:
        r, s, e, L = rid[i], start[i], end[i], lengths[rid[i]]
        lo, hi = np.searchsorted(b_rid, r, "left"), np.searchsorted(b_rid, r, "right")
        bs, be = start[ib[lo:hi]], end[ib[lo:hi]]
        c = cut[ib[lo:hi]]
        okr, okl = (bs >= e) & ~c, (be <= s) & ~c
        right, left = be[okr] - e, s - bs[okl]
        for side, (d, room) in enumerate(((right, L - e), (left, s))):
            if d.size:
                out.append((i, float(d.min()), 1, side, float(room)))
            elif room > 0:
                out.append((i, float(room), 0, side, float(room)))
    return pd.DataFrame(out, columns=["i", "t", "e", "side", "room"])


def km_curve(t, e, wt):
    """Weighted Kaplan-Meier: event times and F = 1 - S at them (events before censorings at ties)."""
    o = np.lexsort((-e, t))
    t, e, wt = t[o], e[o], wt[o]
    at_risk = wt[::-1].cumsum()[::-1]
    use = e == 1
    return t[use], 1 - np.cumprod(1 - wt[use] / at_risk[use])


def cdf_at(curve, x):
    tt, ff = curve
    x = np.asarray(x, float)
    if tt.size == 0:
        return np.zeros(x.shape)
    k = np.searchsorted(tt, x, side="right")
    return np.where(k > 0, ff[np.maximum(k - 1, 0)], 0.0)


def clump_var_w(f: np.ndarray) -> float:
    """Var of sum_i 1[D <= tau_i] for one clump, from F(tau_i): sum_ij min(F_i, F_j) - (sum F)^2."""
    s = np.sort(f)
    n = s.size
    return float(np.sum(s * (2 * (n - np.arange(n)) - 1)) - s.sum() ** 2)


def clump_var_mol(fl: np.ndarray, fr: np.ndarray) -> float:
    """Var of sum_i 1[D_L <= tau_Li or D_R <= tau_Ri] for one clump, from F(tau_L), F(tau_R)."""
    q = (1 - fl) * (1 - fr)
    p00 = (1 - np.maximum.outer(fl, fl)) * (1 - np.maximum.outer(fr, fr))
    return float(np.sum(p00 - np.outer(q, q)))


def mid_p(o: float, sims: np.ndarray) -> tuple[float, float, float]:
    """Monte Carlo mid-p values of an observed count o against null draws: (two-sided, upper, lower)."""
    n = sims.size
    if n == 0:
        return np.nan, np.nan, np.nan
    gt, eq = int(np.sum(sims > o + 0.5)), int(np.sum(np.abs(sims - o) <= 0.5))
    lt = n - gt - eq
    up, low = (gt + (eq + 1) / 2) / (n + 1), (lt + (eq + 1) / 2) / (n + 1)
    return min(1.0, 2 * min(up, low)), up, low


def _chunks(n_items: int, n_mc: int, budget: float = 1e7):
    step = max(1, int(budget // max(n_items, 1)))
    return [(b, min(step, n_mc - b)) for b in range(0, n_mc, step)]


def mc_counts_w(f: np.ndarray, clump: np.ndarray, n_mc: int, rng) -> np.ndarray:
    """Null draws of sum_i 1[D_k(i) <= tau_i], one shared D per clump: 1[U_k <= F(tau_i)] with U_k ~ U(0, 1)."""
    keep = f > 0
    f, clump = f[keep], pd.factorize(clump[keep])[0]
    out = np.zeros(n_mc)
    if f.size == 0:
        return out
    k = int(clump.max()) + 1
    for b0, b in _chunks(f.size, n_mc):
        u = rng.random((k, b), dtype=np.float32)
        out[b0 : b0 + b] = (u[clump] <= f[:, None]).sum(0)
    return out


def mc_counts_mol(fl, fr, clump, shared, n_mc: int, rng) -> np.ndarray:
    """Null draws of sum_i 1[D_Lk <= tau_Li or D_Rk <= tau_Ri] (genome-oriented F per occurrence); clumps with
    shared = True use one draw for both sides and all occurrences: 1[U_k <= 1 - (1 - F_L)(1 - F_R)]."""
    keep = (fl > 0) | (fr > 0)
    fl, fr, shared = fl[keep], fr[keep], shared[keep]
    clump = pd.factorize(clump[keep])[0]
    out = np.zeros(n_mc)
    if fl.size == 0:
        return out
    k = int(clump.max()) + 1
    p = 1 - (1 - fl) * (1 - fr)
    for b0, b in _chunks(fl.size, n_mc, 5e6):
        ul = rng.random((k, b), dtype=np.float32)[clump]
        ur = rng.random((k, b), dtype=np.float32)[clump]
        y = np.where(shared[:, None], ul <= p[:, None], (ul <= fl[:, None]) | (ur <= fr[:, None]))
        out[b0 : b0 + b] = y.sum(0)
    return out


def _z(o, e, v, min_eff, sims=None):
    z = (o - e) / np.sqrt(v) if v > 0 else np.nan
    eff = e * e / v if v > 0 else np.nan
    if sims is None:  # normal reference (kept for comparison)
        p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
        up = float(stats.norm.sf(z)) if np.isfinite(z) else np.nan
    else:
        p, up, _ = mid_p(o, sims)
    return dict(O=o, E=e, V=v, z=z, p=p, p_enrich=up, eff=eff, few=not (eff >= min_eff))


def score_tests(
    sd: pd.DataFrame,
    side_host,
    curves: dict,
    unit_of_side,
    locus,
    parity,
    conflict,
    w: int,
    min_eff: float = 10.0,
    n_mc: int = 9_999,
    seed: int = 0,
) -> dict:
    """sd: sides_dir table; side_host: host label per side; curves: host -> km_curve of its anchor sides;
    unit_of_side: A-unit index per side; locus / parity / conflict: per A unit. n_mc = 0: normal reference.
    Returns {"w": ..., "mol": ...}."""
    rng = np.random.default_rng(seed)
    room = sd.room.to_numpy()
    tau = np.minimum(w, room)
    fw, ff = np.zeros(len(sd)), np.zeros(len(sd))
    for h, ix in pd.Series(np.arange(len(sd))).groupby(np.asarray(side_host)).indices.items():
        if h in curves:
            fw[ix] = cdf_at(curves[h], tau[ix])
            ff[ix] = cdf_at(curves[h], room[ix])
    u = np.asarray(unit_of_side)
    loc, par, bad = np.asarray(locus)[u], np.asarray(parity)[u], np.asarray(conflict)[u]
    oside = np.where(bad, 0, sd.side.to_numpy() ^ par)  # conflicting loci: both sides in one clump
    yw = ((sd.e.to_numpy() == 1) & (sd.t.to_numpy() <= tau)).astype(float)
    fwk, ck = fw, pd.factorize(pd.MultiIndex.from_arrays([loc, oside]))[0]
    vw = sum(clump_var_w(fwk[ix]) for ix in pd.Series(np.arange(ck.size)).groupby(ck).indices.values())
    sims_w = mc_counts_w(fwk, ck, n_mc, rng) if n_mc else None
    res_w = _z(float(yw.sum()), float(fwk.sum()), float(vw), min_eff, sims_w)
    # same molecule: one row per unit occurrence (feature row i), sides as columns (a missing side: F = 0, no event)
    t = pd.DataFrame({"i": sd.i.to_numpy(), "side": sd.side.to_numpy(), "F": ff, "e": sd.e.to_numpy(), "lc": loc})
    t = t.assign(par=par, bad=bad)
    occ = t.pivot_table(index="i", columns="side", values=["F", "e"], aggfunc="first").fillna(0.0)
    f_right = occ[("F", 0)].to_numpy() if ("F", 0) in occ else np.zeros(len(occ))
    f_left = occ[("F", 1)].to_numpy() if ("F", 1) in occ else np.zeros(len(occ))
    ev = sum(occ[("e", s)].to_numpy() for s in (0, 1) if ("e", s) in occ)
    first = t.drop_duplicates("i").set_index("i").reindex(occ.index)
    lo, po = first["lc"].to_numpy(), first["par"].to_numpy().astype(int)
    sh = pd.Series(first["bad"].to_numpy().astype(bool)).groupby(lo).transform("any").to_numpy()  # per locus
    fl, fr = np.where(po == 1, f_right, f_left), np.where(po == 1, f_left, f_right)  # genome-oriented sides
    p_occ = 1 - (1 - fl) * (1 - fr)
    vm = 0.0
    for ix in pd.Series(np.arange(lo.size)).groupby(lo).indices.values():
        vm += clump_var_w(p_occ[ix]) if sh[ix[0]] else clump_var_mol(fl[ix], fr[ix])
    sims_m = mc_counts_mol(fl, fr, lo, sh, n_mc, rng) if n_mc else None
    res_m = _z(float((ev > 0).sum()), float(p_occ.sum()), float(vm), min_eff, sims_m)
    return {"w": res_w, "mol": res_m}
