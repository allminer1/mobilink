"""The net linked test (amendment 04 draft v2; D39): H0 net linked = 0, with host-stratified anchors.

Effect sizes:
  net_w   = 1 - (1 - F_A(w)) / (1 - F_R(w))       far-end distance within w (frozen Part B estimator)
  net_mol = 1 - (1 - M_A) / (1 - M_R)             M = weighted share of sides whose molecule carries any B call
with A the A-call sides (gene-length weights 1 / P(whole)) and R the anchor sides. Molecule-level co-occurrence is
thus compared with random positions on the same molecules, so read lengths and excluded volume (anchors that overlap
a B call by >= 50% are dropped, as A calls are) enter empirically, without an exposure model.
Anchors: with strata (host labels per molecule), each stratum gets anchors in proportion to its A calls, placed on
its own molecules with weight (L - len)+; without strata, over all molecules. Host co-variation (A-rich hosts are
also B-rich) then cancels (P2 step 2.4, tables/P2/P2_04d_netlinked_null.tsv). Anchors touching any B call are
rejected and redrawn in their stratum (the frozen Part B rule drops only >= 50% overlaps; a partly overlapping anchor
skips that B in the side distances, which biases net linked upward where B is dense).
A units (units_d, default 1 kb; None = genes as in Part B): A calls within d of each other, with no B call
between them, form one unit, and the distances and anchors are those of units (anchor lengths drawn from the unit
spans). Under the null the units, not their member genes, are placed at random; with A in islands the gene-level net
is negative under no linkage (-0.004 in tables/P2/P2_07_feasibility.tsv: island members are flanked by A, and an
island needs a longer free stretch than one gene). For A without clusters units are genes.
Inference: units are strata if given (A sides and anchors of a stratum go together), else locus clusters, else
molecules. SE from the delete-one-unit jackknife (jackknife=True, strata or clusters, <= 2,000 units; recommended
for few, unequal clusters), else the SD of n_boot bootstrap replicates over units; both SEs are reported.
z = estimate / SE with a t reference on G_eff - 1 degrees of freedom, G_eff = (sum n_g)^2 / sum n_g^2 the
effective number of units by A sides (hosts are very unequal in size); 95% interval estimate +- t * SE
(percentile interval reported as *_plo, *_phi).
Flags: "few units" (< 30 units carry A sides), "sides < 30" (sides reaching w), "dense classes" (A and B
calls cover > 20% of the bases of the molecules that carry them: in simulation the read-level estimate drifts upward
once genes fill about a third of the genome, +0.012 at 30x and +0.13 at 55x the G2b density; it is unbiased up to
14x).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from ..features import drop_same_orf
from . import netlinked as nl


def _touches_b(m, s, e, b_by_m):
    """Boolean per anchor: overlaps any B call on its molecule (b_by_m: molecule index -> (starts, ends))."""
    out = np.zeros(m.size, bool)
    for k in np.flatnonzero(np.isin(m, list(b_by_m))):
        bs, be = b_by_m[m[k]]
        out[k] = bool(np.any((bs < e[k]) & (be > s[k])))
    return out


def stratified_anchors(lengths: np.ndarray, strata, a_strata, a_len, n: int, rng, b_by_m=None):
    """Anchors per stratum in proportion to the A calls in it (at least 1 per stratum with A); strata None = one
    stratum. With b_by_m, anchors that touch a B call are rejected and redrawn in the same stratum until its target
    count is met (real calls do not overlap B; a partly overlapping anchor would skip that B and look past it)."""
    st = np.zeros(lengths.size, int) if strata is None else np.asarray(strata)
    ast = np.zeros(len(a_len), int) if strata is None else np.asarray(a_strata)
    rr, rs, re = [], [], []
    counts = pd.Series(ast).value_counts()
    for h, n_h in counts.items():
        m = np.flatnonzero(st == h)
        if m.size == 0:
            continue
        target = max(1, int(round(n * n_h / counts.sum())))
        got_r, got_s, got_e, k = [], [], [], 0
        for _ in range(50):
            r_rid, r_st, r_en = nl.anchors(lengths[m], a_len, max(target - k, 1) * 2, rng)
            mm = m[r_rid]
            keep = ~_touches_b(mm, r_st, r_en, b_by_m) if b_by_m else np.ones(mm.size, bool)
            got_r.append(mm[keep]), got_s.append(r_st[keep]), got_e.append(r_en[keep])
            k += int(keep.sum())
            if k >= target:
                break
        rr.append(np.concatenate(got_r)[:target])
        rs.append(np.concatenate(got_s)[:target])
        re.append(np.concatenate(got_e)[:target])
    return np.concatenate(rr), np.concatenate(rs), np.concatenate(re)


def prepare(
    lengths: pd.DataFrame, fa: pd.DataFrame, fb: pd.DataFrame, strata=None, n_anchors=nl.N_ANCHORS, seed=nl.SEED
):
    """Features and anchors for the net linked test: anchors stratified by host (if strata are given), and anchors
    touching a B call rejected and redrawn (amendment 04 draft v2). The frozen Part B anchors are netlinked.prepare."""
    ln = lengths.length.to_numpy(np.int64)
    idx = pd.Index(lengths.molecule_id)
    fa, fb, n_dual = drop_same_orf(fa, fb)
    bm = idx.get_indexer(fb.molecule_id)
    b_by_m = {k: (g.start.to_numpy(), g.end.to_numpy()) for k, g in fb.assign(m=bm).groupby("m")}
    rng = np.random.default_rng(seed)
    a_m = idx.get_indexer(fa.molecule_id)
    st = None if strata is None else np.asarray(strata)
    r_rid, r_st, r_en = stratified_anchors(
        ln, st, None if st is None else st[a_m], (fa.end - fa.start).to_numpy(np.int64), n_anchors, rng, b_by_m
    )
    rnd = pd.DataFrame({"molecule_id": idx[r_rid], "start": r_st, "end": r_en, "scov": 100.0})
    q = pd.concat([fa.assign(c="A"), rnd.assign(c="RND"), fb.assign(c="B")], ignore_index=True)
    rid = idx.get_indexer(q.molecule_id)
    if (rid < 0).any():
        raise ValueError(f"{int((rid < 0).sum())} features name a molecule missing from the lengths table")
    data = (ln, rid, q.start.to_numpy(np.int64), q.end.to_numpy(np.int64), q.c.to_numpy(str))
    return data, q.scov.to_numpy(float), dict(n_dual=n_dual, A_calls=len(fa), B_calls=len(fb))


def _mol_net(A, R, has_b, ia=slice(None), ir=slice(None)):
    ma = np.average(has_b[A["read"][ia]], weights=A["w"][ia]) if np.size(A["w"][ia]) else np.nan
    mr = np.average(has_b[R["read"][ir]], weights=R["w"][ir]) if np.size(R["w"][ir]) else np.nan
    return 1 - (1 - ma) / (1 - mr) if mr < 1 else np.nan


def a_units(fa: pd.DataFrame, fb: pd.DataFrame, d: int) -> pd.DataFrame:
    """A calls merged into units: consecutive A calls on a molecule join if their gap is <= d and no B call lies
    between them (a unit never spans a B call). Span first start - last end; scov the minimum of its calls."""
    a = fa[["molecule_id", "start", "end", "scov"]].assign(isb=False)
    b = fb[["molecule_id", "start", "end"]].assign(scov=np.nan, isb=True)
    q = pd.concat([a, b], ignore_index=True).sort_values(["molecule_id", "start", "end"], kind="stable")
    q = q.reset_index(drop=True)
    same = (q.molecule_id == q.molecule_id.shift()).to_numpy()
    prev_b = q.isb.shift(fill_value=True).to_numpy()
    run_end = q.assign(e=np.where(q.isb, -np.inf, q.end)).groupby("molecule_id", sort=False).e.cummax().to_numpy()
    prev_end = np.r_[-np.inf, run_end[:-1]]
    cont = same & ~prev_b & (q.start.to_numpy() - prev_end <= d) & ~q.isb.to_numpy()
    qa = q[~q.isb.to_numpy()]
    uid = np.cumsum(~cont[~q.isb.to_numpy()])
    g = qa.assign(uid=uid).groupby("uid", sort=True)
    return pd.DataFrame(
        {"molecule_id": g.molecule_id.first(), "start": g.start.min(), "end": g.end.max(), "scov": g.scov.min()}
    ).reset_index(drop=True)


def expected_events(A, R, w):
    """Observed and null-expected A events within w (a B call within w of the side): O counts the A sides with an
    event at t <= w; E is the log-rank expectation, the sum over anchor event times t <= w of the anchors' weighted
    hazard times the number of A sides still at risk (t_A >= t). Unweighted counts, for the minimum-event rule."""
    ta, ea = A["t"], A["e"]
    o = int(np.sum((ea == 1) & (ta <= w)))
    tr, er, wr = R["t"], R["e"], R["w"]
    k = np.lexsort((-er, tr))
    tr, er, wr = tr[k], er[k], wr[k]
    at_risk = wr[::-1].cumsum()[::-1]
    use = (er == 1) & (tr <= w)
    ta_sorted = np.sort(ta)
    y_a = ta_sorted.size - np.searchsorted(ta_sorted, tr[use], side="left")
    return o, float(np.sum(wr[use] / at_risk[use] * y_a))


def sign_flip(pv: np.ndarray, n_flip: int, rng, alpha: float = 0.05, chunk: int = 1000):
    """Randomization test and interval from jackknife pseudo-values (one per host or cluster).
    H0: value = delta. Statistic |sum(pv - delta)|; null distribution from random signs per unit. It is exact if
    the units' deviations are independent and symmetric, whatever the unit sizes: a dominant unit's sign carries no
    evidence, so the test rests on the other units. A t reference on G_eff - 1 over-rejects when a few hosts carry
    the jackknife variance (P2 step 2.7, run 4: SE coefficient of variation 0.6-0.8 across datasets).
    Returns (p at delta = 0, lower, upper): the interval is {delta: p(delta) > alpha}, found by bisection."""
    pv = pv[np.isfinite(pv)]
    g = pv.size
    if g < 3:
        return np.nan, np.nan, np.nan
    a, r = np.empty(n_flip), np.empty(n_flip)
    for i in range(0, n_flip, chunk):
        e = rng.integers(0, 2, size=(min(chunk, n_flip - i), g)).astype(np.float64) * 2 - 1
        a[i : i + e.shape[0]], r[i : i + e.shape[0]] = e @ pv, e.sum(1)
    tot = pv.sum()

    def p_at(d):
        s = abs(tot - g * d)
        return (1 + np.sum(np.abs(a - d * r) >= s * (1 - 1e-12))) / (n_flip + 1)

    centre = tot / g
    spread = float(np.sqrt(np.sum((pv - centre) ** 2) / (g * (g - 1)))) or 1e-9
    bounds = []
    for sgn in (-1, 1):
        inner, outer = centre, centre + sgn * 4 * spread
        while p_at(outer) > alpha and abs(outer - centre) < 1e3 * spread:
            outer = centre + 2 * (outer - centre)
        for _ in range(40):
            mid = (inner + outer) / 2
            inner, outer = (mid, outer) if p_at(mid) > alpha else (inner, mid)
        bounds.append((inner + outer) / 2)
    return float(p_at(0.0)), float(bounds[0]), float(bounds[1])


def net_linked_test(
    lengths: pd.DataFrame,
    fa: pd.DataFrame,
    fb: pd.DataFrame,
    strata=None,
    clusters=None,
    w: int = 2000,
    n_boot: int = 500,
    seed: int = nl.SEED,
    min_sides: int = 30,
    units_d: int | None = 1000,
    jackknife: bool = True,
    n_flip: int = 4999,
    min_events: float = 10.0,
    dense_cover: float = 0.1,
) -> dict:
    """lengths: molecule_id, length; fa, fb: molecule_id, start, end[, scov]; strata / clusters: one label per
    molecule (lengths order). Returns estimates, bootstrap SE, 95% intervals and two-sided p values."""
    fa = fa.assign(scov=fa["scov"] if "scov" in fa else 100.0)[["molecule_id", "start", "end", "scov"]]
    fb = fb.assign(scov=fb["scov"] if "scov" in fb else 100.0)[["molecule_id", "start", "end", "scov"]]
    if units_d is not None:  # measure from A units: islands of A are what the null places at random
        fa, fb, _ = drop_same_orf(fa, fb)
        fa = a_units(fa, fb, units_d)
    data, scov, info = prepare(lengths[["molecule_id", "length"]], fa, fb, strata, seed=seed)
    ln, rid, st, en, cls = data
    A, R, s = nl.analyse(*data, scov, w, min_sides)
    has_b = np.zeros(ln.size, float)
    has_b[rid[cls == "B"]] = 1.0
    est_w, est_m = s["net"], _mol_net(A, R, has_b)
    if strata is not None:
        unit, unit_kind = pd.factorize(np.asarray(strata))[0], "strata"
    elif clusters is not None:
        unit, unit_kind = pd.factorize(np.asarray(clusters))[0], "clusters"
    else:
        unit, unit_kind = np.arange(ln.size), "molecules"
    ua, ur = unit[A["read"]], unit[R["read"]]
    units = np.unique(ua)
    by_a = pd.Series(np.arange(ua.size)).groupby(ua).apply(np.asarray).to_dict()
    by_r = pd.Series(np.arange(ur.size)).groupby(ur).apply(np.asarray).to_dict()
    empty = np.array([], dtype=np.int64)
    rng = np.random.default_rng(seed + 2)
    bw, bm = np.full(n_boot, np.nan), np.full(n_boot, np.nan)
    joint = strata is not None
    r_units = units if joint else np.unique(ur)
    for b in range(n_boot):
        pick = rng.choice(units, units.size)
        ia = np.concatenate([by_a[u] for u in pick])
        pr = pick if joint else rng.choice(r_units, r_units.size)
        ir = np.concatenate([by_r.get(u, empty) for u in pr])
        if ir.size == 0:
            continue
        bw[b] = nl.net(A, R, w, ia, ir, min_sides)
        bm[b] = _mol_net(A, R, has_b, ia, ir)
    # delete-one-unit jackknife (units carrying A sides; with strata the unit's anchors are removed with it)
    jw, jm = np.full(units.size, np.nan), np.full(units.size, np.nan)
    if jackknife and unit_kind != "molecules" and 2 < units.size <= 2000:
        a_unit_of = ua
        r_unit_of = ur
        for j, u in enumerate(units):
            ia = np.flatnonzero(a_unit_of != u)
            ir = np.flatnonzero(r_unit_of != u) if joint else np.arange(ur.size)
            jw[j] = nl.net(A, R, w, ia, ir, min_sides)
            jm[j] = _mol_net(A, R, has_b, ia, ir)
    n_u = np.bincount(pd.factorize(ua)[0]).astype(float)
    g_eff = float(n_u.sum() ** 2 / np.sum(n_u**2)) if n_u.size else np.nan
    out = dict(
        w=w,
        unit=unit_kind,
        units_with_A=int(units.size),
        units_effective=g_eff,
        **info,
        sides=int(A["t"].size),
        sides_w=int((A["t"] >= w).sum()),
        a_units_d=units_d,
    )
    df = max(g_eff - 1, 1.0) if np.isfinite(g_eff) else 1.0
    q = float(stats.t.ppf(0.975, df))
    o_w, e_w = expected_events(A, R, w)
    m_r = np.average(has_b[R["read"]], weights=R["w"]) if R["w"].size else np.nan
    out.update(
        net_w_events=o_w,
        net_w_expected=e_w,
        net_mol_events=int(has_b[A["read"]].sum()),
        net_mol_expected=float(A["t"].size * m_r),
    )
    frng = np.random.default_rng(seed + 3)
    for name, est, bs, js in (("net_w", est_w, bw, jw), ("net_mol", est_m, bm, jm)):
        jo = js[np.isfinite(js)]
        pv = jo.size * est - (jo.size - 1) * jo if jo.size > 2 and np.isfinite(est) else np.array([])
        out[f"{name}_p_flip"], out[f"{name}_flo"], out[f"{name}_fhi"] = sign_flip(pv, n_flip, frng)
        ok = bs[np.isfinite(bs)]
        se_boot = float(np.std(ok, ddof=1)) if ok.size > 2 else np.nan
        jo = js[np.isfinite(js)]
        se_jk = float(np.sqrt((jo.size - 1) / jo.size * np.sum((jo - jo.mean()) ** 2))) if jo.size > 2 else np.nan
        se = se_jk if np.isfinite(se_jk) else se_boot
        out[f"{name}_se_boot"] = se_boot
        out[f"{name}_se_jk"] = se_jk
        z = est / se if se and np.isfinite(se) and se > 0 else np.nan
        out.update(
            {
                name: est,
                f"{name}_se": se,
                f"{name}_lo": est - q * se if np.isfinite(se) else np.nan,
                f"{name}_hi": est + q * se if np.isfinite(se) else np.nan,
                f"{name}_plo": float(np.quantile(ok, 0.025)) if ok.size else np.nan,
                f"{name}_phi": float(np.quantile(ok, 0.975)) if ok.size else np.nan,
                f"{name}_p": float(2 * stats.t.sf(abs(z), df)) if np.isfinite(z) else np.nan,
                f"{name}_p_normal": float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan,
            }
        )
    carriers = np.unique(rid[(cls == "A") | (cls == "B")])
    called = float(np.sum((en - st)[(cls == "A") | (cls == "B")]))
    out["call_cover_carriers"] = called / float(ln[carriers].sum()) if carriers.size else np.nan
    out["call_cover_all"] = called / float(ln.sum()) if ln.size else np.nan
    flags = []
    if out["call_cover_all"] > dense_cover:
        flags.append(f"dense classes (> {dense_cover:g} of all bases called; read-level net biased upward)")
    if units.size < 30:
        flags.append("few units")
    if not e_w >= min_events:
        flags.append(f"few events (expected A events within w {e_w:.1f} < {min_events:g})")
    if out["sides_w"] < min_sides:
        flags.append("sides < 30")
    if strata is None:
        flags.append("no host strata (anchors over all molecules)")
    out["flags"] = "; ".join(flags)
    return out
