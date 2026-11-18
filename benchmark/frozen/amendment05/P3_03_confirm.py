"""P3 step 3.3: the confirmatory analysis of amendment 05 (Parts B, C and D), from the tables of P3_02.

Input: a run list (run_accession, env, project, unit; the analysed runs after any recorded replacement) and
public/<run>/ tables. Steps:
  1. per run: N50; the registered net linked test (linkstat.pairs.net_test_pairs: Kraken2 species strata, unit loci
     from the carrier overlaps, registered classes from <run>.core.features.tsv.gz, pairs ARG:MGE, ARG:MRG, MRG:MGE,
     w = 2 kb, 9,999 Monte Carlo draws, default seed; BH across the three pairs);
  2. per run: the frozen Part B driver (frozen/partB/P1_06_g1b.py, mobility core) for the A3 columns (naive share,
     native against N50 3 kb);
  3. per class: the pooled registered test (linkstat.pairs.pool_runs; strata run x host; BH across the three pairs);
  4. the decision rules of Part D -> OUT_DIR/verdict.txt and verdict.tsv.
Rules (Part D):
  H1a  counting runs (A3 `counts`); supported if >= 2/3 of them have naive share native / 3 kb >= 1.5; < 6 counting
       runs: not determined.
  H1b  Spearman correlation of naive share (ARG:MGE, Part B `intact_naive`) with subsample N50 within each class, summed
       with weights n - 1; one-sided permutation p (N50 permuted within class, 9,999 permutations, seed 20261116);
       supported at p <= 0.05.
  H2   (under linkage) supported if >= 2/3 of the counting runs are stable (|native - 3 kb| <= max(0.02, 0.25 |native|))
       with overlapping intervals; < 6 counting runs: not determined. (Under no linkage: Spearman of the ARG:MRG
       net_mol with N50, reported.)
  H4   soil runs, ARG:MRG, same-molecule test. Per run (BH across the soil runs) and pooled (BH across the pairs):
       above chance (q <= 0.05, net_mol > 0), below chance (q <= 0.05, net_mol < 0), at chance (90% jackknife interval
       of M_A / M_R within 0.8-1.25), else not determined. Supported if the pooled result is at or below chance and at
       most a quarter of the soil runs (3 of 12) are above chance; contradicted if the pooled result is above chance;
       else not determined.
  H3   descriptive: pooled ARG:MGE net_w and net_mol per class.
usage: python -I P3_03_confirm.py RUN_LIST OUT_DIR [--procs 12] [--h4-class soil] [--steps 1,2,3,4]
"""
import argparse
import os
import subprocess
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

L = "/mnt/rmrfc/workflow/iskander/linkstat"
KDB = f"{L}/db/kraken2/k2_standard_16gb_20241228"
PAIRS = [("ARG", "MGE"), ("ARG", "MRG"), ("MRG", "MGE")]
A_CLASSES = ("ARG", "MRG")
PERM_SEED, N_PERM = 20261116, 9_999


def n50(x):
    x = np.sort(np.asarray(x))[::-1]
    c = np.cumsum(x)
    return int(x[np.searchsorted(c, c[-1] / 2)])


def load_run(run):
    """lengths, features (registered classes), Kraken2 species strata, carrier overlaps among A-class molecules."""
    from linkstat.hosts import read_taxonomy, strata_for_molecules
    from linkstat.loci import read_paf

    P = f"{L}/public/{run}/{run}"
    ln = pd.read_csv(f"{P}.lengths.tsv.gz", sep="\t", dtype={"read_id": str}).rename(columns={"read_id": "molecule_id"})
    ln = ln[ln.length > 0].reset_index(drop=True)
    f = pd.read_csv(f"{P}.core.features.tsv.gz", sep="\t", dtype={"read_id": str}).rename(
        columns={"read_id": "molecule_id"})
    f = f[f["class"].isin(["ARG", "MRG", "MGE"]) & f.molecule_id.isin(ln.molecule_id)]
    f = f[["molecule_id", "start", "end", "class", "scov"]].reset_index(drop=True)
    k = pd.read_csv(f"{P}.kraken.tsv.gz", sep="\t", dtype={"read_id": str})
    st = strata_for_molecules(ln.molecule_id.to_numpy(), k, read_taxonomy(f"{KDB}/ktaxonomy.tsv"), rank="S")
    keep = set(f.loc[f["class"].isin(A_CLASSES), "molecule_id"])
    paf = read_paf(f"{P}.ava.paf.gz", min_block=1000, keep=keep)
    return ln, f, st, paf


def per_run(args):
    run, out = args
    path = f"{out}/per_run/{run}.tsv"
    if os.path.exists(path):
        return path
    from linkstat.pairs import net_test_pairs

    ln, f, st, paf = load_run(run)
    t = net_test_pairs(ln, f, PAIRS, strata=st, paf=paf)
    t.insert(0, "run", run)
    t["subsample_n50"], t["subsample_bases"] = n50(ln.length), int(ln.length.sum())
    for c in ("ARG", "MRG", "MGE"):
        t[f"calls_{c}"] = int((f["class"] == c).sum())
    ad = f"{L}/public/{run}/{run}.adapters.tsv.gz"
    if os.path.exists(ad):  # kit from the read starts: ligation adapter, rapid primer, native barcode
        a = pd.read_csv(ad, sep="\t", dtype={"read_id": str})
        s = a[a.region == "start"]
        for nm, m in (("lsk", s.pattern == "lsk"), ("rad", s.pattern == "rad"), ("nb", s.pattern.str.startswith("nb:"))):
            t[f"start_{nm}"] = s.read_id[m].nunique() / len(ln)
    t.to_csv(path + ".tmp", sep="\t", index=False)
    os.replace(path + ".tmp", path)
    return path


def part_b(args):
    run, out = args
    path = f"{out}/partB/{run}.tsv"
    if not os.path.exists(path):
        r = subprocess.run([sys.executable, "-I", f"{L}/frozen/partB/P1_06_g1b.py", "--runs", run, "--out", path],
                           capture_output=True, text=True)
        if r.returncode != 0:
            open(f"{out}/partB/{run}.err", "w").write(r.stdout + r.stderr)
            return None
    return path


def pooled(args):
    env, runs, out = args
    path = f"{out}/pooled/{env.replace(' ', '_')}.tsv"
    if os.path.exists(path):
        return path
    from linkstat.pairs import net_test_pairs, pool_runs

    data = []
    for run in runs:
        ln, f, st, paf = load_run(run)
        data.append(dict(run=run, lengths=ln, features=f, strata=st, paf=paf))
    ln, f, st, paf = pool_runs(data)
    t = net_test_pairs(ln, f, PAIRS, strata=st, paf=paf, oe=False)
    t.insert(0, "env", env)
    t["runs"] = len(runs)
    t.to_csv(path + ".tmp", sep="\t", index=False)
    os.replace(path + ".tmp", path)
    return path


def bh(p):
    from linkstat.stats.multiple import bh as _bh
    return _bh(np.asarray(p, float))


def classify(q, net_mol, lo90, hi90):
    if np.isfinite(q) and q <= 0.05 and np.isfinite(net_mol) and net_mol > 0:
        return "above chance"
    if np.isfinite(q) and q <= 0.05 and np.isfinite(net_mol) and net_mol < 0:
        return "below chance"
    if np.isfinite(lo90) and np.isfinite(hi90) and lo90 >= 0.8 and hi90 <= 1.25:
        return "at chance"
    return "not determined"


def perm_test(x, n50s, env, rng):
    """Sum over classes of (n - 1) x Spearman(x, N50); one-sided permutation p (N50 permuted within class)."""
    def stat(v):
        s = 0.0
        for e in np.unique(env):
            m = (env == e) & np.isfinite(x)
            if m.sum() >= 3 and np.ptp(x[m]) > 0 and np.ptp(v[m]) > 0:
                s += (m.sum() - 1) * spearmanr(x[m], v[m]).statistic
        return s
    obs = stat(n50s)
    perm = np.empty(N_PERM)
    for i in range(N_PERM):
        v = n50s.copy()
        for e in np.unique(env):
            m = np.flatnonzero(env == e)
            v[m] = v[rng.permutation(m)]
        perm[i] = stat(v)
    return obs, float((1 + np.sum(perm >= obs - 1e-12)) / (N_PERM + 1))


def rules(sel, out, h4_class):
    lines, rows = [], []
    pr = pd.concat([pd.read_csv(f"{out}/per_run/{r}.tsv", sep="\t") for r in sel.run_accession], ignore_index=True)
    pr = pr.merge(sel[["run_accession", "env", "project"]], left_on="run", right_on="run_accession")
    pb = pd.concat([pd.read_csv(f"{out}/partB/{r}.tsv", sep="\t") for r in sel.run_accession
                    if os.path.exists(f"{out}/partB/{r}.tsv")], ignore_index=True)
    pb = pb.merge(sel[["run_accession", "env"]], left_on="run", right_on="run_accession")
    po = pd.concat([pd.read_csv(f"{out}/pooled/{e.replace(' ', '_')}.tsv", sep="\t") for e in sorted(sel.env.unique())],
                   ignore_index=True)
    pr.to_csv(f"{out}/per_run.tsv", sep="\t", index=False)
    pb.to_csv(f"{out}/partB.tsv", sep="\t", index=False)
    po.to_csv(f"{out}/pooled.tsv", sep="\t", index=False)
    lines.append(f"runs analysed: {sel.run_accession.nunique()} ({', '.join(f'{e} {n}' for e, n in sel.env.value_counts().sort_index().items())})")

    # H1a and H2 (A3 rule on the counting runs)
    c = pb[pb.counts.astype(bool)]
    n_c = len(c)
    for name, ok, label in (("H1a", c.naive_ok.astype(bool), "naive share native / 3 kb >= 1.5"),
                            ("H2 (under linkage)", c.stable.astype(bool) & c.overlap.astype(bool),
                             "net linked stable, intervals overlap")):
        share = float(ok.mean()) if n_c else np.nan
        v = "not determined" if n_c < 6 else ("supported" if share >= 2 / 3 else "not supported")
        lines.append(f"{name}: {v} | counting runs {n_c} (classes: {c.env.nunique()}), {label} in {int(ok.sum())} "
                     f"({share:.2f})")
        rows.append(dict(hypothesis=name, verdict=v, counting_runs=n_c, meeting=int(ok.sum()), share=share))
    if n_c:
        lines.append(f"  naive ratio median {c.naive_ratio.median():.2f} (range {c.naive_ratio.min():.2f}-"
                     f"{c.naive_ratio.max():.2f}); |native - 3 kb| median {c['diff'].median():.4f}")

    # H1b
    x = pb.set_index("run")
    rng = np.random.default_rng(PERM_SEED)
    obs, p = perm_test(x.intact_naive.to_numpy(float), x.intact_n50.to_numpy(float), x.env.to_numpy(str), rng)
    rho = {e: spearmanr(g.intact_naive, g.intact_n50).statistic for e, g in pb.groupby("env") if len(g) >= 3}
    v = "supported" if p <= 0.05 else "not supported"
    lines.append(f"H1b: {v} | weighted sum of within-class Spearman {obs:.2f}, one-sided permutation p {p:.4f}; "
                 + ", ".join(f"{e} rho {r:.2f}" for e, r in rho.items()))
    rows.append(dict(hypothesis="H1b", verdict=v, statistic=obs, p=p))

    # H2 under no linkage (descriptive)
    am = pr[(pr.class_a == "ARG") & (pr.class_b == "MRG") & pr.net_mol.notna()]
    if len(am) >= 3:
        r = spearmanr(am.net_mol, am.subsample_n50)
        lines.append(f"H2 (no linkage, descriptive): Spearman of ARG:MRG net_mol with N50 over {len(am)} runs "
                     f"{r.statistic:.2f} (p {r.pvalue:.3f})")

    # H4
    s = pr[(pr.env == h4_class) & (pr.class_a == "ARG") & (pr.class_b == "MRG")].copy()
    if len(s):
        s["q_runs"] = bh(s.net_mol_p)
        s["h4_class"] = [classify(q, m, lo, hi) for q, m, lo, hi in
                         zip(s.q_runs, s.net_mol, s.mol_ratio_lo90, s.mol_ratio_hi90)]
        s[["run", "project", "net_mol", "net_mol_lo", "net_mol_hi", "net_mol_p", "q_runs", "mol_ratio",
           "mol_ratio_lo90", "mol_ratio_hi90", "h4_class", "flags"]].to_csv(f"{out}/H4_runs.tsv", sep="\t", index=False)
        g = po[(po.env == h4_class) & (po.class_a == "ARG") & (po.class_b == "MRG")].iloc[0]
        pooled_class = classify(g.net_mol_q, g.net_mol, g.mol_ratio_lo90, g.mol_ratio_hi90)
        above = int((s.h4_class == "above chance").sum())
        if pooled_class == "above chance":
            v = "contradicted"
        elif pooled_class in ("at chance", "below chance") and above <= len(s) / 4:
            v = "supported"
        else:
            v = "not determined"
        lines.append(f"H4: {v} | pooled {h4_class} ARG:MRG {pooled_class} (net_mol {g.net_mol:.4f}, q {g.net_mol_q:.4f}, "
                     f"M_A/M_R {g.mol_ratio:.3f}, 90% {g.mol_ratio_lo90:.3f}-{g.mol_ratio_hi90:.3f}); runs: "
                     + ", ".join(f"{k} {n}" for k, n in s.h4_class.value_counts().items()))
        rows.append(dict(hypothesis="H4", verdict=v, pooled=pooled_class, runs_above=above, runs=len(s)))

    # H3 (descriptive)
    for _, g in po[(po.class_a == "ARG") & (po.class_b == "MGE")].iterrows():
        lines.append(f"H3 (descriptive) {g.env}: ARG:MGE net_w {g.net_w:.4f} ({g.net_w_lo:.4f} to {g.net_w_hi:.4f}), "
                     f"net_mol {g.net_mol:.4f} ({g.net_mol_lo:.4f} to {g.net_mol_hi:.4f}), q {g.net_w_q:.4g} / "
                     f"{g.net_mol_q:.4g}")
    pd.DataFrame(rows).to_csv(f"{out}/verdict.tsv", sep="\t", index=False)
    open(f"{out}/verdict.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_list")
    ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=12)
    ap.add_argument("--h4-class", default="soil")
    ap.add_argument("--steps", default="1,2,3,4")
    a = ap.parse_args()
    steps = set(a.steps.split(","))
    sel = pd.read_csv(a.run_list, sep="\t", dtype=str)
    for d in ("per_run", "partB", "pooled"):
        os.makedirs(f"{a.out}/{d}", exist_ok=True)
    runs = sel.run_accession.tolist()
    if "1" in steps:
        with Pool(a.procs) as p:
            p.map(per_run, [(r, a.out) for r in runs], chunksize=1)
    if "2" in steps:
        with Pool(a.procs) as p:
            res = p.map(part_b, [(r, a.out) for r in runs], chunksize=1)
        bad = [r for r, x in zip(runs, res) if x is None]
        if bad:
            sys.exit(f"Part B driver failed for: {', '.join(bad)}")
    if "3" in steps:
        jobs = [(e, g.run_accession.tolist(), a.out) for e, g in sel.groupby("env")]
        with Pool(min(a.procs, len(jobs))) as p:
            p.map(pooled, jobs, chunksize=1)
    if "4" in steps:
        rules(sel, a.out, a.h4_class)


if __name__ == "__main__":
    main()
