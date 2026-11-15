"""Net linked fraction of class A with class B at window w (amendment 02 Part B estimator).

For each side of each A feature, the event time is the distance to the far end of the nearest B feature on that side;
a side without such a B is censored at the molecule end. A B call cut by a molecule end (ending exactly at it, or
within 50 bp of it with subject cover < 95%) censors the side instead of ending it.
F(w) is a weighted Kaplan-Meier estimate; each A side is weighted by 1 / P(a gene of its length lies whole on a molecule
of this set). The chance term is the same estimate for random anchors (lengths drawn from the A features, placed on
the molecules where they fit, anchors overlapping a B by >= 50% dropped).
    net(w) = 1 - (1 - F_A(w)) / (1 - F_anchor(w)),   NaN if fewer than `min_sides` A sides reach w.
Intervals: percentile bootstrap over molecules (A sides and anchor sides resampled by molecule).
Robustness to read length: the same estimate after random fragmentation to target N50s (several breakpoint draws),
with features cut by a breakpoint kept if >= 80% of them stays on one piece.

Port of P1_03_estimators and P1_06_g1b (frozen in amendment 02 Part B); numerics and random-number order are identical.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..features import drop_same_orf
from ..fragment import fragment_partial, n50, random_bounds, rate_for_n50

SEED = 20261114
N_ANCHORS = 20_000


def sides(lengths, rid, start, end, is_a, is_b, cut=None):
    """Per A side: (A index, time, event). Event time = distance to the far end of the nearest B beyond the A feature
    on that side; censored at the molecule end. `cut`: mask over all features; B calls marked cut censor."""
    ia, ib = np.flatnonzero(is_a), np.flatnonzero(is_b)
    ib = ib[np.argsort(rid[ib], kind="stable")]
    b_rid = rid[ib]
    out_i, out_t, out_e = [], [], []
    for i in ia:
        r, s, e, L = rid[i], start[i], end[i], lengths[rid[i]]
        lo, hi = np.searchsorted(b_rid, r, "left"), np.searchsorted(b_rid, r, "right")
        bs, be = start[ib[lo:hi]], end[ib[lo:hi]]
        okr, okl = bs >= e, be <= s
        if cut is not None:
            c = cut[ib[lo:hi]]
            okr &= ~c
            okl &= ~c
        right, left = be[okr] - e, s - bs[okl]
        for d, room in ((right, L - e), (left, s)):
            if d.size:
                out_i.append(i)
                out_t.append(float(d.min()))
                out_e.append(1)
            elif room > 0:
                out_i.append(i)
                out_t.append(float(room))
                out_e.append(0)
    return np.array(out_i, dtype=np.int64), np.array(out_t), np.array(out_e, dtype=int)


def km(times, events, w, weights=None, min_at_risk=30) -> float:
    """Weighted Kaplan-Meier F(w) = 1 - S(w); NaN if fewer than min_at_risk sides (unweighted) reach w.
    Events come before censorings at tied times."""
    if times.size == 0 or np.sum(times >= w) < min_at_risk:
        return np.nan
    wt = np.ones_like(times) if weights is None else weights
    o = np.lexsort((-events, times))
    t, ev, wt = times[o], events[o], wt[o]
    at_risk = wt[::-1].cumsum()[::-1]
    use = (ev == 1) & (t <= w)
    return float(1 - np.prod(1 - wt[use] / at_risk[use]))


def events_within(times, events, w) -> int:
    """Number of events at distances <= w."""
    return int(np.sum((events == 1) & (times <= w)))


def zero_event_upper(n_reaching_w: int, alpha: float = 0.05) -> float:
    """Exact (Clopper-Pearson) upper confidence limit of F(w) when no event occurs within w among n sides that reach
    w: 1 - (alpha / 2) ** (1 / n). A resampling interval cannot show this uncertainty: every resample has F(w) = 0."""
    if n_reaching_w <= 0:
        return 1.0
    return float(1 - (alpha / 2) ** (1 / n_reaching_w))


def p_whole(lengths, g):
    """Relative probability that a gene of length g lies whole on a molecule of this set:
    sum_r (L_r - g)+ / sum_r L_r."""
    L = np.sort(lengths.astype(float))
    c = np.concatenate(([0.0], np.cumsum(L[::-1])))[::-1]
    k = np.searchsorted(L, g, "right")
    return (c[k] - g * (L.size - k)) / L.sum()


def cut_mask(lengths, rid, start, end, scov, near: int = 50, scov_max: float = 95.0):
    """Calls cut by a molecule end: ending exactly at it, or within `near` bp of it with subject cover < scov_max."""
    d_end = np.minimum(start, lengths[rid] - end)
    return (d_end <= 0) | ((d_end <= near) & (scov < scov_max))


def anchors(lengths, a_len, n: int, rng):
    """n random anchors: length drawn from a_len, molecule chosen with weight (L - len)+ (rejection sampling)."""
    rid, st, ln = [], [], []
    p_read = lengths / lengths.sum()
    while sum(map(len, rid)) < n:
        m = 2 * n
        r = rng.choice(len(lengths), m, p=p_read)
        g = rng.choice(a_len, m)
        ok = rng.random(m) < np.clip(lengths[r] - g, 0, None) / lengths[r]
        r, g = r[ok], g[ok]
        rid.append(r)
        ln.append(g)
        st.append((rng.random(r.size) * (lengths[r] - g)).astype(np.int64))
    rid, st, ln = (np.concatenate(x)[:n] for x in (rid, st, ln))
    return rid, st, st + ln


def side_table(lengths, rid, start, end, cls, cut, who):
    i, t, e = sides(lengths, rid, start, end, cls == who, cls == "B", cut=cut)
    return dict(t=t, e=e, w=1 / p_whole(lengths, (end - start)[i]), read=rid[i])


def net(A, R, w, ia=slice(None), ir=slice(None), min_sides: int = 30) -> float:
    if np.sum(A["t"][ia] >= w) < min_sides:
        return np.nan
    fa = km(A["t"][ia], A["e"][ia], w, A["w"][ia], min_at_risk=1)
    fr = km(R["t"][ir], R["e"][ir], w, R["w"][ir], min_at_risk=1)
    return 1 - (1 - fa) / (1 - fr)


def bootstrap(A, R, w, n_boot: int, rng, min_sides: int = 30):
    """95% percentile interval, resampling molecules for the A sides and the anchor sides independently."""
    ga, gr = pd.factorize(A["read"])[0], pd.factorize(R["read"])[0]
    by_a = np.split(np.argsort(ga, kind="stable"), np.cumsum(np.bincount(ga))[:-1])
    by_r = np.split(np.argsort(gr, kind="stable"), np.cumsum(np.bincount(gr))[:-1])
    out = np.empty(n_boot)
    for b in range(n_boot):
        ia = np.concatenate([by_a[k] for k in rng.integers(0, len(by_a), len(by_a))])
        ir = np.concatenate([by_r[k] for k in rng.integers(0, len(by_r), len(by_r))])
        out[b] = net(A, R, w, ia, ir, min_sides)
    return np.nanquantile(out, [0.025, 0.975])


def analyse(lengths, rid, start, end, cls, scov, w, min_sides: int = 30):
    """Side tables and the summary (N50, molecules with A, naive share, sides, sides reaching w, net) of one version."""
    cut = cut_mask(lengths, rid, start, end, scov)
    A = side_table(lengths, rid, start, end, cls, cut, "A")
    R = side_table(lengths, rid, start, end, cls, cut, "RND")
    has_a = np.zeros(len(lengths), bool)
    has_a[rid[cls == "A"]] = True
    has_b = np.zeros(len(lengths), bool)
    has_b[rid[cls == "B"]] = True
    naive = (has_a & has_b).sum() / has_a.sum() if has_a.sum() else np.nan
    return (
        A,
        R,
        dict(
            n50=n50(lengths),
            reads_A=int(has_a.sum()),
            naive=naive,
            sides=int(A["t"].size),
            sides_w=int((A["t"] >= w).sum()),
            net=net(A, R, w, min_sides=min_sides),
        ),
    )


def prepare(lengths: pd.DataFrame, fa: pd.DataFrame, fb: pd.DataFrame, n_anchors: int = N_ANCHORS, seed: int = SEED):
    """Combine A features, B features and random anchors into arrays (molecule-indexed).
    lengths: molecule_id, length (zero lengths removed); fa, fb: molecule_id, start, end, scov.
    Returns (data, scov, info) with data = (lengths, rid, start, end, cls) and cls in {"A", "RND", "B"}."""
    ln = lengths.length.to_numpy(np.int64)
    idx = pd.Index(lengths.molecule_id)
    fa, fb, n_dual = drop_same_orf(fa, fb)
    rng = np.random.default_rng(seed)
    r_rid, r_st, r_en = anchors(ln, (fa.end - fa.start).to_numpy(np.int64), n_anchors, rng)
    rnd = pd.DataFrame({"molecule_id": idx[r_rid], "start": r_st, "end": r_en, "scov": 100.0})
    rnd, _, _ = drop_same_orf(rnd, fb)
    q = pd.concat([fa.assign(c="A"), rnd.assign(c="RND"), fb.assign(c="B")], ignore_index=True)
    rid = idx.get_indexer(q.molecule_id)
    if (rid < 0).any():
        raise ValueError(f"{int((rid < 0).sum())} features name a molecule missing from the lengths table")
    data = (ln, rid, q.start.to_numpy(np.int64), q.end.to_numpy(np.int64), q.c.to_numpy(str))
    return data, q.scov.to_numpy(float), dict(n_dual=n_dual, A_calls=len(fa), B_calls=len(fb))


def length_robustness(
    lengths: pd.DataFrame,
    fa: pd.DataFrame,
    fb: pd.DataFrame,
    w: int = 2000,
    targets=(3000, 2000),
    draws: int = 20,
    n_boot: int = 500,
    seed: int = SEED,
    n_anchors: int = N_ANCHORS,
    min_sides: int = 30,
) -> dict:
    """Net linked fraction on the intact molecules and after random fragmentation to each target N50.
    Intact: estimate and 95% bootstrap interval (seed + 1). Per target: mean over `draws` breakpoint draws
    (seeds seed + 0 ... seed + draws - 1); half-width sqrt(h^2 + (1.96 sd)^2), h = bootstrap half-width at the first
    draw (seed + 2) and sd = between-draw sd. Keys as in the G1b table (intact_*, r<target>_*)."""
    data, scov, info = prepare(lengths, fa, fb, n_anchors, seed)
    ln = data[0]
    row = dict(w=w, **info)
    A, R, s = analyse(*data, scov, w, min_sides)
    lo, hi = bootstrap(A, R, w, n_boot, np.random.default_rng(seed + 1), min_sides)
    row.update({f"intact_{k}": v for k, v in s.items()}, intact_lo=lo, intact_hi=hi)
    for t in targets:
        rate = rate_for_n50(ln, t, seed)
        vals = []
        for k in range(draws):
            dd, keep = fragment_partial(*data, random_bounds(ln, rate, seed + k), return_keep=True)
            A, R, s = analyse(*dd, scov[keep], w, min_sides)
            vals.append(s)
            if k == 0:
                lo, hi = bootstrap(A, R, w, n_boot, np.random.default_rng(seed + 2), min_sides)
                h0 = (hi - lo) / 2
        v = pd.DataFrame(vals)
        sd = v.net.std()
        half = float(np.sqrt(h0**2 + (1.96 * sd) ** 2))
        row.update(
            {
                f"r{t}_n50": int(v.n50.mean()),
                f"r{t}_naive": v.naive.mean(),
                f"r{t}_sides_w": int(v.sides_w.min()),
                f"r{t}_net": v.net.mean(),
                f"r{t}_net_sd_draws": sd,
                f"r{t}_half": half,
                f"r{t}_lo": v.net.mean() - half,
                f"r{t}_hi": v.net.mean() + half,
            }
        )
    return row
