"""In-silico fragmentation of molecules: Poisson breakpoints with a rate tuned to a target N50.

Port of P0_04_gonogo (n50, random_bounds, rate_for_n50, fragment_random) and P1_03_estimators.fragment_partial
(frozen in amendment 02 Part B). Numerics and random-number order are kept identical to the frozen code.
"""

from __future__ import annotations

import numpy as np


def n50(x) -> int:
    s = np.sort(x)[::-1]
    return int(s[np.searchsorted(np.cumsum(s), s.sum() / 2)])


def random_bounds(lengths, rate: float, seed: int):
    """Piece boundaries in concatenated coordinates for Poisson breakpoints at `rate` per bp.
    Returns (offsets of the molecules, sorted unique boundaries including every molecule start and the total end)."""
    rng = np.random.default_rng(seed)
    off = np.concatenate(([0], np.cumsum(lengths)[:-1]))
    k = rng.poisson(lengths * rate)
    rr = np.repeat(np.arange(len(lengths)), k)
    pos = off[rr] + 1 + (rng.random(k.sum()) * (lengths[rr] - 1)).astype(np.int64)
    return off, np.unique(np.concatenate((off, pos, [off[-1] + lengths[-1]])))


def rate_for_n50(lengths, target: int, seed: int, max_reads: int = 300_000, iters: int = 14) -> float:
    """Breakpoint rate whose fragments have N50 closest to `target` (bisection on log rate), tuned on a random
    subsample of at most `max_reads` molecules (the fragment N50 depends on the length distribution, not the count)."""
    if len(lengths) > max_reads:
        lengths = np.random.default_rng(seed).choice(lengths, max_reads, replace=False)
    lo, hi = np.log(1e-7), np.log(1e-1)
    for _ in range(iters):
        mid = (lo + hi) / 2
        _, b = random_bounds(lengths, np.exp(mid), seed)
        if n50(np.diff(b)) > target:
            lo = mid
        else:
            hi = mid
    return float(np.exp((lo + hi) / 2))


def fragment_whole(lengths, rid, start, end, cls, bounds):
    """Cut molecules at the given boundaries; keep only features that stay whole.
    Returns (piece lengths, piece index, start, end, cls)."""
    off, b = bounds
    gs, ge = off[rid] + start, off[rid] + end
    piece = np.searchsorted(b, gs, "right") - 1
    whole = (np.searchsorted(b, ge - 1, "right") - 1) == piece
    return (np.diff(b), piece[whole], (gs - b[piece])[whole], (ge - b[piece])[whole], cls[whole])


def fragment_partial(lengths, rid, start, end, cls, bounds, min_frac: float = 0.8, return_keep: bool = False):
    """As fragment_whole, but a feature cut by a breakpoint is kept on the piece holding the larger part if that part
    is >= min_frac of its length (truncated to the piece), as a caller's subject-cover rule would on real reads."""
    off, b = bounds
    gs, ge = off[rid] + start, off[rid] + end
    p0 = np.searchsorted(b, gs, "right") - 1
    p1 = np.searchsorted(b, ge - 1, "right") - 1
    glen = (end - start).astype(float)
    left_part = b[p0 + 1] - gs
    right_part = ge - b[p1]
    use_start = (p0 == p1) | (left_part >= right_part)
    piece = np.where(use_start, p0, p1)
    s_new = np.where(use_start, gs, b[p1]) - b[piece]
    e_new = np.where(use_start, np.minimum(ge, b[p0 + 1]), ge) - b[piece]
    keep = (p0 == p1) | ((e_new - s_new) >= min_frac * glen)
    out = (np.diff(b), piece[keep], s_new[keep], e_new[keep], cls[keep])
    return (out, keep) if return_keep else out
