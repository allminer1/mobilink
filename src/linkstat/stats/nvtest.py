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
z = (O - E) / sqrt(sum_k Var_k), two-sided normal p. Effective events E^2 / Var; below 10 the p is flagged.
Calibration (no linkage; probe v3, 1,080 datasets with exact clumps): 5.9%, SD of z 1.07.
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


def _z(o, e, v, min_eff):
    z = (o - e) / np.sqrt(v) if v > 0 else np.nan
    eff = e * e / v if v > 0 else np.nan
    return dict(
        O=o,
        E=e,
        V=v,
        z=z,
        p=float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan,
        eff=eff,
        few=not (eff >= min_eff),
    )


def score_tests(
    sd: pd.DataFrame, side_host, curves: dict, unit_of_side, locus, parity, conflict, w: int, min_eff: float = 10.0
) -> dict:
    """sd: sides_dir table; side_host: host label per side; curves: host -> km_curve of its anchor sides;
    unit_of_side: A-unit index per side; locus / parity / conflict: per A unit. Returns {"w": ..., "mol": ...}."""
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
    vw = pd.Series(fw).groupby([loc, oside]).apply(lambda x: clump_var_w(x.to_numpy())).sum()
    res_w = _z(float(yw.sum()), float(fw.sum()), float(vw), min_eff)
    # same molecule: one row per unit occurrence (feature row i), sides as columns (a missing side: F = 0, no event)
    t = pd.DataFrame({"i": sd.i.to_numpy(), "side": sd.side.to_numpy(), "F": ff, "e": sd.e.to_numpy(), "lc": loc})
    occ = t.pivot_table(index="i", columns="side", values=["F", "e"], aggfunc="first").fillna(0.0)
    fr = occ[("F", 0)].to_numpy() if ("F", 0) in occ else np.zeros(len(occ))
    fl = occ[("F", 1)].to_numpy() if ("F", 1) in occ else np.zeros(len(occ))
    ev = sum(occ[("e", s)].to_numpy() for s in (0, 1) if ("e", s) in occ)
    lo = t.drop_duplicates("i").set_index("i")["lc"].reindex(occ.index).to_numpy()
    p_occ = 1 - (1 - fl) * (1 - fr)
    vm = sum(clump_var_mol(fl[ix], fr[ix]) for ix in pd.Series(np.arange(lo.size)).groupby(lo).indices.values())
    res_m = _z(float((ev > 0).sum()), float(p_occ.sum()), float(vm), min_eff)
    return {"w": res_w, "mol": res_m}
