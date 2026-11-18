"""The net linked test (amendment 04 draft v2; D39, D42, D44): H0 no linkage, with host-stratified anchors.

Effect sizes (descriptive, with host-jackknife intervals):
  net_w   = 1 - (1 - F_A(w)) / (1 - F_R(w))      far-end distance within w (frozen Part B estimator), A units, with
                                                 the margin rule of net_mol: units (and anchors) ending >= d from
                                                 both molecule ends, weight 1 / p_whole(span + 2 d). Without it, a
                                                 partly covered island still forms a unit whose outer side carries
                                                 the weight of a whole unit of the observed span, so unlinked
                                                 islands outweighed linked singletons: net_w 0.026 against truth
                                                 0.032 with A islands and linkage in soil (run 6, probes v6-v7, D45).
  net_mol = 1 - (1 - M_A) / (1 - M_R)            M = weighted share of A units (or anchors) on a molecule that carries a
                                                 B call; margin rule: units ending >= d from both molecule ends, anchors
                                                 with the same margin (a unit cut by a molecule end looks complete but
                                                 sits on a read chosen by its span, which biased net_mol by -0.01 with A
                                                 in islands; P2 step 2.7, probes v1-v2)
Tests (registered p values): null-variance score tests, linkstat.stats.nvtest, over A-unit sides against the host's
anchor CDF, with clumps from unit loci (linkstat.unitloci): net_w (window w, all A-unit sides; clump = unit
locus x oriented side) and same molecule (clump = unit locus); p values from the exact null
distribution of the clump model by Monte Carlo (n_mc draws, mid-p). Run 6 (normal reference): 3.7% and 4.3%.
Anchors: with strata (host labels per molecule), each stratum gets anchors in proportion to its A units
(anchors_per_unit, default 40, at least netlinked.N_ANCHORS in all), placed on its own molecules with weight
(L - len)+; without strata, over all molecules. Anchors touching any B call are redrawn in their stratum.
A units (units_d, default 1 kb; None = genes as in Part B): A calls within d of each other, with no B call between
them. With A in islands the gene-level net is negative under no linkage (-0.004, P2_07_feasibility.tsv).
Flags: few effective events (E^2 / V < min_eff, per test), dense classes (A and B calls cover > dense_cover of all
molecule bases: in simulated genomes without other genes, crowded A and B genes pack against each other and create
association under no linkage), no host strata, no locus clumps (each unit occurrence its own clump: anti-conservative
when reads cover a locus more than once), locus conflicts (> 10% of units in loci whose overlaps disagree on
orientation; their sides are merged, which is conservative), few units (< 30 hosts carry A units), sides < 30.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from ..features import drop_same_orf
from ..unitloci import unit_loci_from_paf, unit_loci_from_positions
from . import netlinked as nl
from . import nvtest as nv


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
    lengths: pd.DataFrame,
    fa: pd.DataFrame,
    fb: pd.DataFrame,
    strata=None,
    n_anchors=nl.N_ANCHORS,
    seed=nl.SEED,
    same_orf: float | None = 0.5,
):
    """Features and anchors for the net linked test: anchors stratified by host (if strata are given), and anchors
    touching a B call rejected and redrawn (amendment 04 draft v2). The frozen Part B anchors are netlinked.prepare.
    same_orf: overlap (share of the shorter) at which A and B calls count as one ORF and are removed; None = off
    (annotation sensitivity only; the registered rule is 0.5)."""
    ln = lengths.length.to_numpy(np.int64)
    idx = pd.Index(lengths.molecule_id)
    fa, fb, n_dual = drop_same_orf(fa, fb, same_orf) if same_orf is not None else (fa, fb, 0)
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


def loci_from_paf(paf: pd.DataFrame):
    """Unit-locus function for net_linked_test from ava overlaps (linkstat.loci.read_paf with min_block 1000)."""
    return lambda units: unit_loci_from_paf(units, paf)


def loci_by_molecule(units: pd.DataFrame) -> pd.DataFrame:
    """Unit-locus function for contigs (loci="molecule"): all A units of one molecule are one clump, both sides
    merged. Units on one contig share the event "the contig carries B", so they are not independent draws under no
    linkage. The estimates do not change; the tests get the wider null variance."""
    return pd.DataFrame({"locus": pd.factorize(units.molecule_id)[0], "parity": 0, "conflict": True})


def loci_from_positions(lengths: pd.DataFrame):
    """Unit-locus function for simulated reads (lengths with replicon and pos columns; reads are forward)."""
    idx = pd.Index(lengths.molecule_id)
    rep, pos = lengths.replicon.to_numpy(), lengths.pos.to_numpy()

    def f(units):
        m = idx.get_indexer(units.molecule_id)
        lab = unit_loci_from_positions(rep[m], pos[m] + units.start.to_numpy(), pos[m] + units.end.to_numpy())
        return pd.DataFrame({"locus": lab, "parity": 0, "conflict": False})

    return f


def _margin_mol(ln, idx, fa, fb, st, d, n_per, seed):
    """Margin-rule inputs for net_mol: kept A units (outcome, weight, host) and their anchors."""
    um = idx.get_indexer(fa.molecule_id)
    keep = (fa.start.to_numpy() >= d) & (ln[um] - fa.end.to_numpy() >= d)
    am = um[keep]
    span = (fa.end - fa.start).to_numpy(np.int64)[keep]
    has_b = np.zeros(ln.size, bool)
    bm = idx.get_indexer(fb.molecule_id)
    has_b[bm] = True
    if am.size == 0:
        return None
    lnm = np.maximum(ln - 2 * d, 0)
    b_by_m = {k: (g.start.to_numpy() - d, g.end.to_numpy() - d) for k, g in fb.assign(m=bm).groupby("m")}
    r_m, r_s, r_e = stratified_anchors(
        lnm, st, st[am], span, max(nl.N_ANCHORS, n_per * am.size), np.random.default_rng(seed + 5), b_by_m
    )
    return dict(
        ya=has_b[am].astype(float),
        wa=1 / nl.p_whole(lnm, span),
        ha=st[am],
        yr=has_b[r_m].astype(float),
        wr=1 / nl.p_whole(lnm, r_e - r_s),
        hr=st[r_m],
    )


def _mol_from_sums(sa, wa, sr, wr):
    ma, mr = sa / wa, sr / wr
    return 1 - (1 - ma) / (1 - mr) if wa > 0 and wr > 0 and mr < 1 else np.nan


def _ratio_from_sums(sa, wa, sr, wr):
    """log(M_A / M_R): the same-molecule share of A units over that of the anchors (amendment 04, 2.6: the "O/E" of
    H3 and H4)."""
    return np.log((sa / wa) / (sr / wr)) if wa > 0 and wr > 0 and sa > 0 and sr > 0 else np.nan


def _ci(est, se, df):
    q = float(stats.t.ppf(0.975, df))
    return (est - q * se, est + q * se) if np.isfinite(se) else (np.nan, np.nan)


def net_linked_test(
    lengths: pd.DataFrame,
    fa: pd.DataFrame,
    fb: pd.DataFrame,
    strata=None,
    loci=None,
    w: int = 2000,
    seed: int = nl.SEED,
    min_sides: int = 30,
    units_d: int | None = 1000,
    anchors_per_unit: int = 40,
    n_boot: int = 200,
    min_eff: float = 10.0,
    dense_cover: float = 0.05,
    n_mc: int = 9_999,
    same_orf: float | None = 0.5,
) -> dict:
    """lengths: molecule_id, length[, replicon, pos]; fa, fb: molecule_id, start, end[, scov]; strata: one host label
    per molecule (lengths order) or None; loci: None, "positions" (simulated reads), "molecule" (contigs: one clump
    per molecule), or a function units -> DataFrame (locus, parity, conflict) such as loci_from_paf(paf). same_orf: the same-ORF overlap rule (0.5 registered; None =
    off, for annotation sensitivity). Returns estimates, intervals and the registered p values."""
    fa = fa.assign(scov=fa["scov"] if "scov" in fa else 100.0)[["molecule_id", "start", "end", "scov"]]
    fb = fb.assign(scov=fb["scov"] if "scov" in fb else 100.0)[["molecule_id", "start", "end", "scov"]]
    fa, fb, n_dual = drop_same_orf(fa, fb, same_orf) if same_orf is not None else (fa, fb, 0)
    if units_d is not None:  # measure from A units: islands of A are what the null places at random
        fa = a_units(fa, fb, units_d)
    fa = fa.reset_index(drop=True)
    ln_df = lengths[["molecule_id", "length"]]
    idx = pd.Index(lengths.molecule_id)
    ln = lengths.length.to_numpy(np.int64)
    st = np.zeros(ln.size, int) if strata is None else pd.factorize(np.asarray(strata))[0]
    n_anch = max(nl.N_ANCHORS, anchors_per_unit * len(fa))
    data, scov, info = prepare(
        ln_df, fa, fb, None if strata is None else st, n_anchors=n_anch, seed=seed, same_orf=same_orf
    )
    info["n_dual"] = n_dual
    _, rid, s_, e_, cls = data
    cut = nl.cut_mask(ln, rid, s_, e_, scov)
    span = e_ - s_
    sd = nv.sides_dir(ln, rid, s_, e_, cls == "A", cls == "B", cut)
    sr = nv.sides_dir(ln, rid, s_, e_, cls == "RND", cls == "B", cut)

    def side_dict(t, keep, add):
        i = t.i.to_numpy()[keep]
        return dict(t=t.t.to_numpy()[keep], e=t.e.to_numpy()[keep], w=1 / nl.p_whole(ln, span[i] + add), read=rid[i])

    margin = units_d or 0  # margin rule (units and anchors ending >= d from both molecule ends)
    inner = np.minimum(s_, ln[rid] - e_)
    keep_a, keep_r = inner[sd.i.to_numpy()] >= margin, inner[sr.i.to_numpy()] >= margin
    A, R = side_dict(sd, keep_a, 2 * margin), side_dict(sr, keep_r, 2 * margin)  # net_w effect size
    R_all = side_dict(sr, np.ones(len(sr), bool), 0)  # chance curves for the tests: all anchor sides
    # registered tests: null variance with unit-locus clumps
    if loci is None:
        lt = pd.DataFrame({"locus": np.arange(len(fa)), "parity": 0, "conflict": False})
    elif isinstance(loci, str) and loci == "positions":
        lt = loci_from_positions(lengths)(fa)
    elif isinstance(loci, str) and loci == "molecule":
        lt = loci_by_molecule(fa)
    elif isinstance(loci, str):
        raise ValueError(f"loci: unknown value {loci!r} (use 'positions', 'molecule' or a function)")
    else:
        lt = loci(fa)
    hr = st[R_all["read"]]
    curves = {
        h: nv.km_curve(R_all["t"][ix], R_all["e"][ix], R_all["w"][ix])
        for h, ix in pd.Series(hr).groupby(hr).indices.items()
    }
    tests = nv.score_tests(
        sd,
        st[rid[sd.i.to_numpy()]],
        curves,
        sd.i.to_numpy(),
        lt.locus.to_numpy(),
        lt.parity.to_numpy(),
        lt.conflict.to_numpy(),
        w,
        min_eff,
        n_mc=n_mc,
        seed=seed + 7,
    )
    # effect sizes: net_w (frozen estimator) and net_mol, both with the margin rule; host-jackknife intervals
    est_w = nl.net(A, R, w, min_sides=min_sides)
    mm = _margin_mol(ln, idx, fa, fb, st, units_d or 0, anchors_per_unit // 2, seed)
    if mm is not None:
        SA = pd.Series(mm["wa"] * mm["ya"]).groupby(mm["ha"]).sum()
        WA = pd.Series(mm["wa"]).groupby(mm["ha"]).sum()
        SR = pd.Series(mm["wr"] * mm["yr"]).groupby(mm["hr"]).sum()
        WR = pd.Series(mm["wr"]).groupby(mm["hr"]).sum()
        est_m = _mol_from_sums(SA.sum(), WA.sum(), SR.sum(), WR.sum())
        est_lr = _ratio_from_sums(SA.sum(), WA.sum(), SR.sum(), WR.sum())
        m_a, m_r = float(SA.sum() / WA.sum()), float(SR.sum() / WR.sum())
    else:
        est_m = est_lr = m_a = m_r = np.nan
    ua, ur = st[A["read"]], st[R["read"]]
    hosts = np.unique(ua)
    n_u = np.bincount(pd.factorize(ua)[0]).astype(float)
    g_eff = float(n_u.sum() ** 2 / np.sum(n_u**2)) if n_u.size else np.nan
    df = max(g_eff - 1, 1.0) if np.isfinite(g_eff) else 1.0
    se_w = se_m = se_lr = np.nan
    if strata is not None and 2 < hosts.size:
        jw = np.array([nl.net(A, R, w, np.flatnonzero(ua != h), np.flatnonzero(ur != h), min_sides) for h in hosts])
        jo = jw[np.isfinite(jw)]
        se_w = float(np.sqrt((jo.size - 1) / jo.size * np.sum((jo - jo.mean()) ** 2))) if jo.size > 2 else np.nan
        if mm is not None:
            jm = np.array(
                [
                    _mol_from_sums(
                        SA.sum() - SA.get(h, 0),
                        WA.sum() - WA.get(h, 0),
                        SR.sum() - SR.get(h, 0),
                        WR.sum() - WR.get(h, 0),
                    )
                    for h in hosts
                ]
            )
            jo = jm[np.isfinite(jm)]
            se_m = float(np.sqrt((jo.size - 1) / jo.size * np.sum((jo - jo.mean()) ** 2))) if jo.size > 2 else np.nan
            jl = np.array(
                [
                    _ratio_from_sums(
                        SA.sum() - SA.get(h, 0),
                        WA.sum() - WA.get(h, 0),
                        SR.sum() - SR.get(h, 0),
                        WR.sum() - WR.get(h, 0),
                    )
                    for h in hosts
                ]
            )
            jo = jl[np.isfinite(jl)]
            se_lr = float(np.sqrt((jo.size - 1) / jo.size * np.sum((jo - jo.mean()) ** 2))) if jo.size > 2 else np.nan
    else:  # no strata: bootstrap over molecules
        rng = np.random.default_rng(seed + 2)
        by_a = pd.Series(np.arange(A["read"].size)).groupby(A["read"]).indices
        mols = np.array(list(by_a))
        bs = []
        for _ in range(n_boot):
            ia = np.concatenate([by_a[m] for m in rng.choice(mols, mols.size)])
            bs.append(nl.net(A, R, w, ia, slice(None), min_sides))
        bs = np.array(bs)
        se_w = float(np.nanstd(bs, ddof=1)) if np.isfinite(bs).sum() > 2 else np.nan
    out = dict(
        w=w,
        hosts_with_A=int(hosts.size),
        hosts_effective=g_eff,
        **info,
        a_units_d=units_d,
        anchors=n_anch,
        units=len(fa),
        unit_loci=int(lt.locus.nunique()),
        loci_conflict_share=float(lt.conflict.mean()) if len(lt) else 0.0,
        margin=margin,
        sides_all=int(len(sd)),
        sides=int(A["t"].size),
        sides_w=int((A["t"] >= w).sum()),
    )
    for name, est, se, key in (("net_w", est_w, se_w, "w"), ("net_mol", est_m, se_m, "mol")):
        lo, hi = _ci(est, se, df)
        r = tests[key]
        out.update(
            {
                name: est,
                f"{name}_se": se,
                f"{name}_lo": lo,
                f"{name}_hi": hi,
                f"{name}_p": r["p"],
                f"{name}_p_enrich": r["p_enrich"],
                f"{name}_z": r["z"],
                f"{name}_O": r["O"],
                f"{name}_E": r["E"],
                f"{name}_V": r["V"],
                f"{name}_eff": r["eff"],
            }
        )
    # same-molecule ratio M_A / M_R with its 90% host-jackknife interval on the log scale (for the registered TOST)
    q90 = float(stats.t.ppf(0.95, df))
    ok = np.isfinite(est_lr) and np.isfinite(se_lr)
    out.update(
        mol_share_A=m_a,
        mol_share_R=m_r,
        mol_ratio=float(np.exp(est_lr)) if np.isfinite(est_lr) else np.nan,
        mol_log_ratio_se=se_lr,
        mol_ratio_lo90=float(np.exp(est_lr - q90 * se_lr)) if ok else np.nan,
        mol_ratio_hi90=float(np.exp(est_lr + q90 * se_lr)) if ok else np.nan,
    )
    called = float(np.sum((e_ - s_)[(cls == "A") | (cls == "B")]))
    out["call_cover_all"] = called / float(ln.sum()) if ln.size else np.nan
    flags = []
    for name, key in (("net_w", "w"), ("net_mol", "mol")):
        if tests[key]["few"]:
            flags.append(f"{name}: few effective events (E^2/V {tests[key]['eff']:.1f} < {min_eff:g})")
    if out["call_cover_all"] > dense_cover:
        flags.append(f"dense classes (> {dense_cover:g} of all bases called)")
    if strata is None:
        flags.append("no host strata (anchors over all molecules)")
    if loci is None:
        flags.append("no locus clumps (each unit occurrence is its own clump)")
    if isinstance(loci, str) and loci == "molecule":
        flags.append("molecule clumps (all A units of a molecule are one clump)")
    elif out["loci_conflict_share"] > 0.1:
        flags.append("locus conflicts (> 10% of units; sides merged)")
    if strata is not None and hosts.size < 30:
        flags.append("few units (< 30 hosts carry A units)")
    if out["sides_w"] < min_sides:
        flags.append("sides < 30")
    out["flags"] = "; ".join(flags)
    return out
