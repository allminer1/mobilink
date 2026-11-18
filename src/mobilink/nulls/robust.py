"""Locus-aware molecule-level test (P2 step 2.4; amendment 04 draft, Part 2).

Within the N2-ev-u cells (1 kb bins of the modelled class's unit room x strata), O - E is a sum of molecule scores
    u_i = (A_i - a_c / n_c) (B_i - b_c / n_c),     sum_i u_i = O - E,
with A the modelled class, B the fixed class and a_c, b_c, n_c the cell counts. Under the null (A placed at random
within cells) E[u_i] = 0. Reads that cover the same locus have correlated scores, so the variance is the cluster-robust
    V_cr = G / (G - 1) * sum_g (sum_{i in g} u_i)^2         (g = locus clusters, G of them)
and never below the independence (hypergeometric) variance V_perm. z = (O - E) / sqrt(max(V_cr, V_perm)), with a
continuity correction of 0.5. The reported p value is the larger of the normal cluster-robust p and the exact N2-ev-u
p: each is valid in its regime (many loci / one read per locus), so the larger one is valid in both.
deff = V_cr / V_perm is the design effect (about 1 without pseudo-replication).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
from scipy import stats

from .binned import cells, length_bins


@dataclass
class RobustResult:
    observed: int
    expected: float
    var_perm: float
    var_cr: float
    deff: float
    clusters: int
    z: float
    p_two_cr: float
    p_enriched_cr: float
    p_depleted_cr: float
    p_two_sided: float  # max(cluster-robust, exact)
    p_enriched: float
    p_depleted: float

    def as_dict(self):
        return asdict(self)


def cluster_robust_test(
    h_mod, h_fix, room, clusters, strata=None, exact=None, width: int = 1000, cap_quantile: float = 0.995
) -> RobustResult:
    """h_mod, h_fix: boolean per molecule; room: the bin variable (unit room of the modelled class); clusters: locus
    cluster code per molecule; exact: the N2-ev-u BinnedResult on the same cells (for the max rule), or None."""
    A, B = np.asarray(h_mod, float), np.asarray(h_fix, float)
    c = cells(length_bins(room, width, cap_quantile), strata)
    n = np.bincount(c).astype(float)
    a, b = np.bincount(c, A), np.bincount(c, B)
    u = (A - (a / n)[c]) * (B - (b / n)[c])
    obs = int((A * B).sum())
    E = float(obs - u.sum())
    vp = float(np.sum(np.where(n > 1, a * b * (n - a) * (n - b) / (n**2 * np.maximum(n - 1, 1)), 0.0)))
    g = pd.factorize(np.asarray(clusters))[0]
    G = int(g.max()) + 1 if g.size else 0
    sg = np.bincount(g, u)
    vcr = float(G / (G - 1) * np.sum(sg**2)) if G > 1 else np.nan
    v = max(vp, vcr) if np.isfinite(vcr) else vp
    d = obs - E
    if v > 0:
        z = d / np.sqrt(v)
        p_hi = float(stats.norm.sf((d - 0.5) / np.sqrt(v)))
        p_lo = float(stats.norm.cdf((d + 0.5) / np.sqrt(v)))
    else:
        z, p_hi, p_lo = np.nan, 1.0, 1.0
    p2 = min(1.0, 2 * min(p_hi, p_lo))
    if exact is not None:
        e_hi, e_lo, e2 = exact.p_enriched, exact.p_depleted, exact.p_two_sided
    else:
        e_hi = e_lo = e2 = 0.0
    return RobustResult(
        observed=obs,
        expected=E,
        var_perm=vp,
        var_cr=vcr,
        deff=vcr / vp if vp > 0 and np.isfinite(vcr) else np.nan,
        clusters=G,
        z=float(z),
        p_two_cr=p2,
        p_enriched_cr=p_hi,
        p_depleted_cr=p_lo,
        p_two_sided=max(p2, e2),
        p_enriched=max(p_hi, e_hi),
        p_depleted=max(p_lo, e_lo),
    )
