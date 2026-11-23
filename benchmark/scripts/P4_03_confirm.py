"""P4 (amendment 07, Parts B to D): the between-environment analysis. Run once on the final list of accepted runs.

usage: python -I P4_03_confirm.py RUNS.tsv OUT_DIR [--procs 8] [--skip-cut] [--min-projects 6]
  RUNS.tsv: run_accession, env (class), project, for the accepted runs; the per-run tables are in public/<run>/.
Per run (Part B): the per-run test and the cut columns of the frozen amendment 05 script, unchanged (imported).
Per class (Part C):
  full       all runs pooled (pool_runs), net_test_pairs with unit loci and 9,999 draws, as amendment 05 Part C;
  jackknife  the pooled ARG:MGE estimate without each project in turn (no loci, 99 draws: only the estimate is used);
             SE = jackknife SE over projects, t interval with projects - 1 degrees of freedom;
  projects   runs of a project pooled (no loci, 99 draws), for the random-effects summary (secondary);
  cut        the pooled estimate after the reads of every run are cut to N50 3 kb (one draw; pieces keep the stratum).
Rules (Part D): H3a (wastewater - soil, human gut - soil; Welch t, one-sided, Holm), H3b (omnibus chi-square over the
class estimates; pairwise contrasts), H1c (within-project slopes on log2 N50), H4b (ARG:MRG same molecule per class).
Outputs in OUT_DIR: per_run.tsv, partB.tsv, pooled.tsv, jackknife.tsv, projects.tsv, cut.tsv, classes.tsv,
contrasts.tsv, h1c.tsv, verdict.tsv, verdict.txt.
"""
import argparse
import importlib.util
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats

L = "/mnt/rmrfc/workflow/iskander/linkstat"
KDB = f"{L}/db/kraken2/k2_standard_16gb_20241228"
PAIRS = [("ARG", "MGE"), ("ARG", "MRG"), ("MRG", "MGE")]
SEED = 20261114
MIN_PROJECTS = 6
H3A = [("wastewater", "soil"), ("human gut", "soil")]
H1C_CLASSES = ("human gut", "livestock", "wastewater", "activated sludge")


def frozen():
    spec = importlib.util.spec_from_file_location("p3_confirm", f"{L}/frozen/amendment05/P3_03_confirm.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def per_run_w(args):
    return frozen().per_run(args)


def part_b_w(args):
    return frozen().part_b(args)


def load_light(run):
    """The loader of the frozen script without the overlaps: lengths, registered features, Kraken2 species strata."""
    from mobilink.hosts import read_taxonomy, strata_for_molecules

    P = f"{L}/public/{run}/{run}"
    ln = pd.read_csv(f"{P}.lengths.tsv.gz", sep="\t", dtype={"read_id": str}).rename(columns={"read_id": "molecule_id"})
    ln = ln[ln.length > 0].reset_index(drop=True)
    f = pd.read_csv(f"{P}.core.features.tsv.gz", sep="\t", dtype={"read_id": str}).rename(columns={"read_id": "molecule_id"})
    f = f[f["class"].isin(["ARG", "MRG", "MGE"]) & f.molecule_id.isin(ln.molecule_id)]
    f = f[["molecule_id", "start", "end", "class", "scov"]].reset_index(drop=True)
    k = pd.read_csv(f"{P}.kraken.tsv.gz", sep="\t", dtype={"read_id": str})
    st = strata_for_molecules(ln.molecule_id.to_numpy(), k, read_taxonomy(f"{KDB}/ktaxonomy.tsv"), rank="S")
    return ln, f, st


def fast_net_w(lengths, fa, fb, strata, w=2000, units_d=1000, anchors_per_unit=40, same_orf=0.5, min_sides=30):
    """The net_w estimate of mobilink.stats.nettest.net_linked_test and nothing else: the same calls in the same order
    (same-ORF rule, A units, stratified anchors with the default seed, sides, margin rule, weighted Kaplan-Meier), without
    the tests and without the delete-one-host jackknife, which the replicates do not use. The value is identical to
    net_linked_test(...)["net_w"]; rules() checks it against the registered pooling for every class."""
    from mobilink.features import drop_same_orf
    from mobilink.stats import netlinked as nl
    from mobilink.stats import nvtest as nv
    from mobilink.stats.nettest import a_units, prepare

    fa = fa.assign(scov=fa["scov"] if "scov" in fa else 100.0)[["molecule_id", "start", "end", "scov"]]
    fb = fb.assign(scov=fb["scov"] if "scov" in fb else 100.0)[["molecule_id", "start", "end", "scov"]]
    fa, fb, _ = drop_same_orf(fa, fb, same_orf)
    fa = a_units(fa, fb, units_d).reset_index(drop=True)
    ln = lengths.length.to_numpy(np.int64)
    st = pd.factorize(np.asarray(strata))[0]
    data, scov, _ = prepare(lengths[["molecule_id", "length"]], fa, fb, st, n_anchors=max(nl.N_ANCHORS, anchors_per_unit * len(fa)),
                            seed=nl.SEED, same_orf=same_orf)
    _, rid, s_, e_, cls = data
    cut = nl.cut_mask(ln, rid, s_, e_, scov)
    span = e_ - s_
    sd = nv.sides_dir(ln, rid, s_, e_, cls == "A", cls == "B", cut)
    sr = nv.sides_dir(ln, rid, s_, e_, cls == "RND", cls == "B", cut)

    def side_dict(t, keep, add):
        i = t.i.to_numpy()[keep]
        return dict(t=t.t.to_numpy()[keep], e=t.e.to_numpy()[keep], w=1 / nl.p_whole(ln, span[i] + add), read=rid[i])

    inner = np.minimum(s_, ln[rid] - e_)
    A = side_dict(sd, inner[sd.i.to_numpy()] >= units_d, 2 * units_d)
    R = side_dict(sr, inner[sr.i.to_numpy()] >= units_d, 2 * units_d)
    return dict(net_w=nl.net(A, R, w, min_sides=min_sides), units=len(fa), sides_w=int((A["t"] >= w).sum()),
                events_w=int(((A["e"] == 1) & (A["t"] <= w)).sum()))


def raw_arg_mge(data):
    """Pooled ARG:MGE net_w without the reporting rule (a replicate must not drop out when it has no event)."""
    from mobilink.pairs import pool_runs

    ln, f, st, _ = pool_runs(data)
    return fast_net_w(ln, f[f["class"] == "ARG"], f[f["class"] == "MGE"], st)


def class_full(args):
    env, runs, out = args
    path = f"{out}/pooled/{env.replace(' ', '_')}.tsv"
    if os.path.exists(path):
        return path
    from mobilink.pairs import net_test_pairs, pool_runs

    fz, data = frozen(), []
    for run in runs:
        ln, f, st, paf = fz.load_run(run)
        data.append(dict(run=run, lengths=ln, features=f, strata=st, paf=paf))
    ln, f, st, paf = pool_runs(data)
    t = net_test_pairs(ln, f, PAIRS, strata=st, paf=paf, oe=False)
    t.insert(0, "env", env)
    t["runs"] = len(runs)
    t.to_csv(path + ".tmp", sep="\t", index=False)
    os.replace(path + ".tmp", path)
    return path


def class_jack(args):
    env, sel, out, skip_cut = args
    pj, pc, pp = (f"{out}/{d}/{env.replace(' ', '_')}.tsv" for d in ("jackknife", "cut", "projects"))
    if all(os.path.exists(x) for x in (pj, pp)) and (skip_cut or os.path.exists(pc)):
        return pj
    from mobilink.fragment import fragment_partial, n50, random_bounds, rate_for_n50
    from mobilink.pairs import net_test_pairs, pool_runs

    data = {}
    for r in sel.itertuples():
        ln, f, st = load_light(r.run_accession)
        data[r.run_accession] = dict(run=r.run_accession, lengths=ln, features=f, strata=st)
    proj = sel.groupby("project").run_accession.apply(list)
    rows = [dict(env=env, left_out="none", runs=len(data), **raw_arg_mge(list(data.values())))]
    for p, rr in proj.items():
        rest = [d for k, d in data.items() if k not in rr]
        rows.append(dict(env=env, left_out=p, runs=len(rest), **(raw_arg_mge(rest) if rest else {})))
    pd.DataFrame(rows).to_csv(pj, sep="\t", index=False)
    # project-level estimates (reporting rule applied)
    rows = []
    for p, rr in proj.items():
        ln, f, st, _ = pool_runs([data[k] for k in rr])
        t = net_test_pairs(ln, f, [("ARG", "MGE")], strata=st, loci=None, n_mc=99, oe=False).iloc[0]
        rows.append(dict(env=env, project=p, runs=len(rr), net_w=t.get("net_w"), net_w_se=t.get("net_w_se"), net_mol=t.get("net_mol"),
                         units=t.get("units"), flags=t.get("flags")))
    pd.DataFrame(rows).to_csv(pp, sep="\t", index=False)
    if not skip_cut:
        cut = []
        for k, d in data.items():
            ln, f = d["lengths"], d["features"]
            L_ = ln.length.to_numpy(np.int64)
            if n50(L_) <= 3000:
                cut.append(d)
                continue
            rid = pd.Index(ln.molecule_id).get_indexer(f.molecule_id)
            off, b = random_bounds(L_, rate_for_n50(L_, 3000, SEED), SEED)
            (pl, piece, s, e, cl), keep = fragment_partial(L_, rid, f.start.to_numpy(np.int64), f.end.to_numpy(np.int64),
                                                           f["class"].to_numpy(str), (off, b), return_keep=True)
            parent = np.searchsorted(off, b[:-1], "right") - 1
            cut.append(dict(run=k, lengths=pd.DataFrame({"molecule_id": np.arange(pl.size).astype(str), "length": pl}),
                            features=pd.DataFrame({"molecule_id": piece.astype(str), "start": s, "end": e, "class": cl,
                                                   "scov": f.scov.to_numpy(float)[keep]}),
                            strata=np.asarray(d["strata"], dtype=object)[parent]))
        pd.DataFrame([dict(env=env, **raw_arg_mge(cut))]).to_csv(pc, sep="\t", index=False)
    return pj


def jack_se(full, reps):
    reps = np.asarray(reps, float)
    reps = reps[np.isfinite(reps)]
    k = len(reps)
    if k < 2 or not np.isfinite(full):
        return np.nan, k
    return float(np.sqrt((k - 1) / k * np.sum((reps - reps.mean()) ** 2))), k


def dl_hk(y, se):
    """DerSimonian-Laird mean with the Hartung-Knapp interval."""
    y, se = np.asarray(y, float), np.asarray(se, float)
    ok = np.isfinite(y) & np.isfinite(se) & (se > 0)
    y, v = y[ok], se[ok] ** 2
    k = len(y)
    if k < 2:
        return dict(k=k)
    w = 1 / v
    fixed = np.sum(w * y) / w.sum()
    q = float(np.sum(w * (y - fixed) ** 2))
    tau2 = max(0.0, (q - (k - 1)) / (w.sum() - np.sum(w ** 2) / w.sum()))
    ws = 1 / (v + tau2)
    mu = float(np.sum(ws * y) / ws.sum())
    se_hk = float(np.sqrt(np.sum(ws * (y - mu) ** 2) / ((k - 1) * ws.sum())))
    t = stats.t.ppf(0.975, k - 1)
    return dict(k=k, re_mean=mu, re_lo95=mu - t * se_hk, re_hi95=mu + t * se_hk, tau=float(np.sqrt(tau2)))


def slope_fe(d, y, x="x", group="g", cluster="project"):
    """OLS slope of y on x with one intercept per group; cluster-robust (CR1) standard error; t with clusters - 1 df."""
    d = d[np.isfinite(d[y]) & np.isfinite(d[x])]
    d = d[d.groupby(group)[y].transform("size") >= 2]
    n, G = len(d), d[cluster].nunique()
    if n < 4 or G < 3:
        return dict(runs=n, projects=G)
    xc, yc = d[x] - d.groupby(group)[x].transform("mean"), d[y] - d.groupby(group)[y].transform("mean")
    sxx = float((xc ** 2).sum())
    b = float((xc * yc).sum() / sxx)
    e = yc - b * xc
    k = d[group].nunique() + 1
    meat = float(((xc * e).groupby(d[cluster]).sum() ** 2).sum())
    se = float(np.sqrt(meat / sxx ** 2 * G / (G - 1) * (n - 1) / max(n - k, 1)))
    t90 = stats.t.ppf(0.95, G - 1)
    return dict(runs=n, projects=G, slope=b, se=se, lo90=b - t90 * se, hi90=b + t90 * se, p_positive=float(stats.t.sf(b / se, G - 1)))


def rules(sel, out, skip_cut):
    fz = frozen()
    lines, verdicts = [], []
    pr = pd.concat([pd.read_csv(f"{out}/per_run/{r}.tsv", sep="\t") for r in sel.run_accession], ignore_index=True)
    pr = pr.drop(columns=[c for c in ("env", "project", "run_accession") if c in pr]).merge(sel[["run_accession", "env", "project"]], left_on="run", right_on="run_accession")
    pb = pd.concat([pd.read_csv(f"{out}/partB/{r}.tsv", sep="\t") for r in sel.run_accession if os.path.exists(f"{out}/partB/{r}.tsv")], ignore_index=True)
    pb = pb.drop(columns=[c for c in ("env", "project", "run_accession") if c in pb]).merge(sel[["run_accession", "env", "project"]], left_on="run", right_on="run_accession")
    envs = sorted(sel.env.unique())
    rd = lambda d: pd.concat([pd.read_csv(f"{out}/{d}/{e.replace(' ', '_')}.tsv", sep="\t") for e in envs if os.path.exists(f"{out}/{d}/{e.replace(' ', '_')}.tsv")], ignore_index=True)  # noqa: E731
    po, jk, pj = rd("pooled"), rd("jackknife"), rd("projects")
    for t, name in ((pr, "per_run"), (pb, "partB"), (po, "pooled"), (jk, "jackknife"), (pj, "projects")):
        t.to_csv(f"{out}/{name}.tsv", sep="\t", index=False)
    cut = None if skip_cut else rd("cut")
    if cut is not None:
        cut.to_csv(f"{out}/cut.tsv", sep="\t", index=False)
    lines.append(f"runs analysed: {sel.run_accession.nunique()} in {sel.project.nunique()} projects ("
                 + ", ".join(f"{e} {int((sel.env == e).sum())} runs / {sel[sel.env == e].project.nunique()} projects" for e in envs) + ")")

    # Part C: class estimates
    rows = []
    for e in envs:
        g = jk[jk.env == e]
        full = float(g[g.left_out == "none"].net_w.iloc[0])
        se, k = jack_se(full, g[g.left_out != "none"].net_w)
        n_proj = sel[sel.env == e].project.nunique()
        u = float(g[g.left_out == "none"].units.iloc[0]) - g[g.left_out != "none"].units.to_numpy(float)   # A units per project
        eff = float(u.sum() ** 2 / np.sum(u ** 2)) if np.sum(u ** 2) > 0 else np.nan
        a = po[(po.env == e) & (po.class_a == "ARG") & (po.class_b == "MGE")].iloc[0]
        if np.isfinite(a.net_w) and abs(a.net_w - full) > 1e-9:
            lines.append(f"!! class {e}: the estimate-only path gives {full:.10f}, the registered pooling {a.net_w:.10f}")
        t = stats.t.ppf(0.975, k - 1) if k >= 2 else np.nan
        row = dict(env=e, projects=n_proj, projects_effective=eff, runs=int((sel.env == e).sum()), analysed=bool(n_proj >= MIN_PROJECTS),
                   net_w=full, se_jack=se, lo95=full - t * se, hi95=full + t * se, df=k - 1,
                   net_w_registered_pooling=a.net_w, net_w_p=a.net_w_p, units=int(a.units),
                   naive_native=float(np.average(pb[pb.env == e].intact_naive.fillna(0), weights=pb[pb.env == e].intact_reads_A.clip(lower=0))) if (pb.env == e).any() else np.nan,
                   n50_median=float(pb[pb.env == e].intact_n50.median()) if (pb.env == e).any() else np.nan,
                   **{k2: v for k2, v in dl_hk(pj[pj.env == e].net_w, pj[pj.env == e].net_w_se).items()},
                   projects_with_estimate=int(np.isfinite(pj[pj.env == e].net_w).sum()))
        if cut is not None and (cut.env == e).any():
            row["net_w_cut_3kb"] = float(cut[cut.env == e].net_w.iloc[0])
        rows.append(row)
    C = pd.DataFrame(rows).set_index("env")
    C.to_csv(f"{out}/classes.tsv", sep="\t")
    for e, r in C.iterrows():
        lines.append(f"class {e}: net_w {r.net_w:.4f} ({r.lo95:.4f} to {r.hi95:.4f}; jackknife over {int(r.df) + 1} projects, {r.projects_effective:.1f} effective by ARG units)"
                     + ("" if r.analysed else " [not analysed: fewer than 6 projects]"))

    # contrasts
    def contrast(a, b):
        ra, rb = C.loc[a], C.loc[b]
        d, s = ra.net_w - rb.net_w, np.hypot(ra.se_jack, rb.se_jack)
        df = s ** 4 / (ra.se_jack ** 4 / ra.df + rb.se_jack ** 4 / rb.df) if s > 0 else np.nan
        t = d / s if s > 0 else np.nan
        return dict(a=a, b=b, diff=d, se=s, df=df, t=t, p_one_sided=float(stats.t.sf(t, df)), p_two_sided=float(2 * stats.t.sf(abs(t), df)))

    ok = [e for e in envs if C.analysed[e] and np.isfinite(C.se_jack[e])]
    K = pd.DataFrame([contrast(a, b) for i, a in enumerate(ok) for b in ok[i + 1:]])
    if len(K):
        K["p_two_sided_holm"] = holm(K.p_two_sided.to_numpy())
    K.to_csv(f"{out}/contrasts.tsv", sep="\t", index=False)

    # H3a
    if "soil" not in ok or not any(a in ok for a, _ in H3A):
        v, detail = "not determined", "soil or both other classes are not analysed"
    else:
        hc = [contrast(a, b) for a, b in H3A if a in ok]
        ph = holm(np.array([h["p_one_sided"] for h in hc]))
        sig = int(np.sum(ph <= 0.05))
        v = "supported" if (sig == len(H3A) and len(hc) == len(H3A)) else ("partly supported" if sig >= 1 else "not supported")
        detail = "; ".join(f"{h['a']} - soil {h['diff']:.4f} (SE {h['se']:.4f}, t {h['t']:.2f}, df {h['df']:.1f}, Holm p {p:.4g})" for h, p in zip(hc, ph))
        # jackknife sign check for gate G5
        for a, b in H3A:
            if a in ok:
                ja, jb = jk[(jk.env == a) & (jk.left_out != "none")].net_w, jk[(jk.env == b) & (jk.left_out != "none")].net_w
                worst = min(ja.min() - C.net_w[b], C.net_w[a] - jb.max())
                detail += f"; smallest {a} - soil difference without one project {worst:.4f}"
    lines.append(f"H3a: {v} | {detail}")
    verdicts.append(dict(hypothesis="H3a", verdict=v, detail=detail))

    # H3b
    if len(ok) >= 3:
        w = 1 / C.se_jack[ok] ** 2
        m = float(np.sum(w * C.net_w[ok]) / w.sum())
        q = float(np.sum(w * (C.net_w[ok] - m) ** 2))
        p = float(stats.chi2.sf(q, len(ok) - 1))
        v = "supported" if p <= 0.05 else "not supported"
        detail = f"Q {q:.2f}, df {len(ok) - 1}, p {p:.4g}; classes by net_w: " + " > ".join(C.net_w[ok].sort_values(ascending=False).index)
        if "net_w_cut_3kb" in C:
            detail += "; cut to 3 kb: " + " > ".join(C.net_w_cut_3kb[ok].sort_values(ascending=False).index)
        detail += "; by naive share: " + " > ".join(C.naive_native[ok].sort_values(ascending=False).index)
    else:
        v, detail = "not determined", "fewer than 3 classes analysed"
    lines.append(f"H3b: {v} | {detail}")
    verdicts.append(dict(hypothesis="H3b", verdict=v, detail=detail))

    # H1c
    a = pr[(pr.class_a == "ARG") & (pr.class_b == "MGE")][["run", "env", "project", "net_w"]].merge(pb[["run", "intact_naive", "intact_n50", "intact_reads_A"]], on="run")
    a = a[a.env.isin(H1C_CLASSES) & (a.intact_reads_A >= 50)].assign(x=lambda d: np.log2(d.intact_n50), g=lambda d: d.env + "|" + d.project)
    sn, sw = slope_fe(a, "intact_naive"), slope_fe(a, "net_w")
    pd.DataFrame([dict(outcome="naive share", **sn), dict(outcome="net_w", **sw)]).to_csv(f"{out}/h1c.tsv", sep="\t", index=False)
    if "slope" not in sn or "slope" not in sw:
        v, detail = "not determined", "too few projects with two or more runs"
    else:
        i = sn["p_positive"] <= 0.05
        inside, outside = (sw["lo90"] >= -0.04 and sw["hi90"] <= 0.04), (sw["lo90"] > 0.04 or sw["hi90"] < -0.04)
        v = "not supported" if not i else ("supported" if inside else ("contradicted" if outside else "not determined for the net linked fraction"))
        detail = (f"naive share slope {sn['slope']:.4f} per doubling of N50 (SE {sn['se']:.4f}, one-sided p {sn['p_positive']:.4g}; {sn['runs']} runs, "
                  f"{sn['projects']} projects); net_w slope {sw['slope']:.4f} (90% {sw['lo90']:.4f} to {sw['hi90']:.4f}; {sw['runs']} runs, {sw['projects']} projects)")
    lines.append(f"H1c: {v} | {detail}")
    verdicts.append(dict(hypothesis="H1c", verdict=v, detail=detail))

    # H4b
    m = po[(po.class_a == "ARG") & (po.class_b == "MRG") & po.env.isin(ok)].copy()
    if len(m):
        m["q"] = fz.bh(m.net_mol_p.to_numpy(float))
        for r in m.itertuples():
            o = fz.classify(r.q, r.net_mol, r.mol_ratio_lo90, r.mol_ratio_hi90)
            lines.append(f"H4b {r.env}: {o} | same molecule {r.net_mol_O:.0f} observed / {r.net_mol_E:.1f} expected, q {r.q:.4g}, share ratio {r.mol_ratio:.2f} (90% {r.mol_ratio_lo90:.2f}-{r.mol_ratio_hi90:.2f})")
            verdicts.append(dict(hypothesis=f"H4b {r.env}", verdict=o, detail=""))
    pd.DataFrame(verdicts).to_csv(f"{out}/verdict.tsv", sep="\t", index=False)
    open(f"{out}/verdict.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]))
    out = np.empty(len(p))
    out[o] = adj
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs")
    ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--skip-cut", action="store_true")
    ap.add_argument("--min-projects", type=int, default=6, help="6 is the registered value; lower only for the feasibility run on development data")
    a = ap.parse_args()
    global MIN_PROJECTS
    MIN_PROJECTS = a.min_projects
    sel = pd.read_csv(a.runs, sep="\t", dtype=str)[["run_accession", "env", "project"]]
    for d in ("per_run", "partB", "pooled", "jackknife", "projects", "cut"):
        os.makedirs(f"{a.out}/{d}", exist_ok=True)
    envs = sorted(sel.env.unique())
    with Pool(a.procs) as p:
        p.map(per_run_w, [(r, a.out) for r in sel.run_accession], chunksize=1)
        p.map(part_b_w, [(r, a.out) for r in sel.run_accession], chunksize=1)
    with Pool(min(a.procs, 2 * len(envs))) as p:
        p.map(class_jack, [(e, sel[sel.env == e], a.out, a.skip_cut) for e in envs], chunksize=1)
        p.map(class_full, [(e, sel[sel.env == e].run_accession.tolist(), a.out) for e in envs], chunksize=1)
    rules(sel, a.out, a.skip_cut)


if __name__ == "__main__":
    main()
