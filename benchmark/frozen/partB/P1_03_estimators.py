#!/usr/bin/env python3
"""P1 step 1.3: candidate net-linked estimators (development; not frozen).

Distance definitions for an A side:
  start  prototype: event time = distance to the START of the nearest B; censored at (read end - mean B length)
  end    event time = distance to the far END of the nearest whole B; censored at the read end.
         On a read, a B is seen exactly when the read reaches its far end, so "event seen" and "not yet censored"
         are the same condition and censoring by the read end is independent of where B lies.
Chance terms:
  poisson   exp(-lambda_B * w) as in the prototype (start definition only)
  anchored  the same Kaplan-Meier estimate for random anchors (lengths drawn from the A genes, placed uniformly on the
            reads where they fit, anchors overlapping a B by >= 50% dropped), fragmented together with A and B
Weights:
  none | ipw  ipw weights each A gene by 1 / P(a gene of its length lies whole on a read of this read set),
              P(g) proportional to sum_r max(0, L_r - g); anchors get the same weights from their own lengths.
net(w) = 1 - (1 - F_A(w)) / (1 - F_chance(w)).
"""
import numpy as np
import pandas as pd


def sides(lengths, rid, start, end, is_a, is_b, mode="end", len_b=None, margin=None, cut=None):
    """Per A side: (a_index, time, event).
    margin (end mode only): a B whose far end lies within `margin` bp of the read end is treated as cut by the read end
    (a partial call), so it censors the side instead of ending it; None = every B call is an event.
    cut (end mode only): boolean mask over all features; B calls marked cut censor instead of ending (see cut_mask)."""
    ia, ib = np.flatnonzero(is_a), np.flatnonzero(is_b)
    ib = ib[np.argsort(rid[ib], kind="stable")]
    b_rid = rid[ib]
    out_i, out_t, out_e = [], [], []
    for i in ia:
        r, s, e, L = rid[i], start[i], end[i], lengths[rid[i]]
        lo, hi = np.searchsorted(b_rid, r, "left"), np.searchsorted(b_rid, r, "right")
        bs, be = start[ib[lo:hi]], end[ib[lo:hi]]
        if mode == "start":
            right, left = bs[bs >= e] - e, s - be[be <= s]
            rooms = (L - e, s) if len_b is None else (L - e - len_b, s - len_b)
        else:
            okr, okl = bs >= e, be <= s
            if margin is not None:
                okr &= be < L - margin; okl &= bs > margin
            if cut is not None:
                c = cut[ib[lo:hi]]; okr &= ~c; okl &= ~c
            right, left = be[okr] - e, s - bs[okl]                   # far end of B beyond A, either side
            rooms = (L - e, s)
        for d, room in ((right, rooms[0]), (left, rooms[1])):
            if d.size:
                out_i.append(i); out_t.append(float(d.min())); out_e.append(1)
            elif room > 0:
                out_i.append(i); out_t.append(float(room)); out_e.append(0)
    return np.array(out_i, dtype=np.int64), np.array(out_t), np.array(out_e, dtype=int)


def km(times, events, w, weights=None, min_at_risk=30):
    """Weighted Kaplan-Meier F(w) = 1 - S(w); NaN if fewer than min_at_risk sides (unweighted) reach w."""
    if times.size == 0 or np.sum(times >= w) < min_at_risk:
        return np.nan
    wt = np.ones_like(times) if weights is None else weights
    o = np.lexsort((-events, times))                           # events before censorings at tied times
    t, ev, wt = times[o], events[o], wt[o]
    at_risk = wt[::-1].cumsum()[::-1]
    use = (ev == 1) & (t <= w)
    return float(1 - np.prod(1 - wt[use] / at_risk[use]))


def p_whole(lengths, g):
    """Relative probability that a gene of length g lies whole on a read of this read set."""
    L = np.sort(lengths.astype(float)); c = np.concatenate(([0.0], np.cumsum(L[::-1])))[::-1]  # c[k] = sum L[k:]
    k = np.searchsorted(L, g, "right")
    return (c[k] - g * (L.size - k)) / L.sum()


def cut_mask(lengths, rid, start, end, scov, near=50, scov_max=95.0):
    """B calls cut by a read end: ending exactly at the read (or fragment) end, or within `near` bp of it with subject
    cover < scov_max (DIAMOND partial alignments at read ends)."""
    d_end = np.minimum(start, lengths[rid] - end)
    return (d_end <= 0) | ((d_end <= near) & (scov < scov_max))


def fragment_partial(lengths, rid, start, end, cls, bounds, min_frac=0.8, return_keep=False):
    """Random fragmentation as P0_04.fragment_random, but a feature cut by a breakpoint is kept on the piece holding
    the larger part if that part is >= min_frac of its length (truncated to the piece), mimicking the caller's
    subject-cover rule (DIAMOND, 80%) on native reads. Returns (piece lengths, piece, start, end, cls)."""
    off, b = bounds
    gs, ge = off[rid] + start, off[rid] + end
    p0 = np.searchsorted(b, gs, "right") - 1
    p1 = np.searchsorted(b, ge - 1, "right") - 1
    glen = (end - start).astype(float)
    left_part = b[p0 + 1] - gs                                    # bases on the start piece
    right_part = ge - b[p1]                                       # bases on the end piece
    use_start = (p0 == p1) | (left_part >= right_part)
    piece = np.where(use_start, p0, p1)
    s_new = np.where(use_start, gs, b[p1]) - b[piece]
    e_new = np.where(use_start, np.minimum(ge, b[p0 + 1]), ge) - b[piece]
    keep = (p0 == p1) | ((e_new - s_new) >= min_frac * glen)
    out = (np.diff(b), piece[keep], s_new[keep], e_new[keep], cls[keep])
    return (out, keep) if return_keep else out


def anchors(lengths, a_len, n, rng):
    """n random anchors: length drawn from a_len, read chosen with weight (L - len)+ (rejection sampling)."""
    rid, st, ln = [], [], []
    p_read = lengths / lengths.sum()
    while sum(map(len, rid)) < n:
        m = 2 * n
        r = rng.choice(len(lengths), m, p=p_read); g = rng.choice(a_len, m)
        ok = rng.random(m) < np.clip(lengths[r] - g, 0, None) / lengths[r]
        r, g = r[ok], g[ok]
        rid.append(r); ln.append(g); st.append((rng.random(r.size) * (lengths[r] - g)).astype(np.int64))
    rid, st, ln = (np.concatenate(x)[:n] for x in (rid, st, ln))
    return rid, st, st + ln
