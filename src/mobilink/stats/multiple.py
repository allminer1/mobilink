"""Multiple testing: Benjamini-Hochberg adjusted p values (q values), NaN-aware."""

from __future__ import annotations

import numpy as np


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg step-up adjusted p values; NaN entries stay NaN and are not counted."""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0:
        return q
    pv = p[ok]
    o = np.argsort(pv)
    ranked = pv[o] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(adj, 1.0)
    q[ok] = out
    return q
