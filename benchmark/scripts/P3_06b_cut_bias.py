"""P3 step 3.5, exploratory (outside the registration): why is the net linked fraction slightly lower after the reads
are cut to 3 kb (0.5-6% in 8 of 8 counting runs with a net linked fraction >= 0.06)?

Hypothesis: the A sides are weighted by 1 / P(a gene of this length lies whole on a piece), but a gene cut by a
breakpoint is kept when >= 80% of it stays on one piece. Long genes are then kept more often than their weight assumes,
so they are over-weighted in the cut version. If long ARG are linked to MGE less often than short ARG, the cut
estimate falls. Checks, per run, with the functions and seeds of the frozen Part B script (frozen/partB/P1_06_g1b.py):
  registered   the registered 3 kb version (80% rule), 20 draws: must reproduce r3000_net of the registered table
  whole        the same breakpoints, but only genes that stay whole are kept (retention then matches the weights)
  unweighted   the 80% rule with all weights 1 (intact and cut)
  by length    the intact estimate for ARG shorter and longer than the median ARG call of the run
Nothing registered is recomputed or overwritten. Output: tables/P3/P3_06b_cut_bias.tsv.
usage: python -I P3_06b_cut_bias.py [--procs 8]
"""
import argparse
import importlib.util
from multiprocessing import Pool

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
RUNS = ["SRR23278716", "SRR23278658", "SRR23278671", "SRR29980985", "SRR33225186", "DRR351862", "SRR16329124",
        "SRR29786052"]
W, DRAWS = 2000, 20


def g1b():
    spec = importlib.util.spec_from_file_location("g1b", f"{L}/frozen/partB/P1_06_g1b.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def net_of(m, d, scov, unweighted=False, a_mask=None):
    """net(w) of one version; a_mask: boolean over the features, A calls outside it are dropped."""
    ln, rid, st, en, cl = d
    if a_mask is not None:
        k = (cl != "A") | a_mask
        rid, st, en, cl, scov = rid[k], st[k], en[k], cl[k], scov[k]
    cut = m.est.cut_mask(ln, rid, st, en, scov)
    A, R = m.side_table(ln, rid, st, en, cl, cut, "A"), m.side_table(ln, rid, st, en, cl, cut, "RND")
    if unweighted:
        A["w"], R["w"] = np.ones_like(A["w"]), np.ones_like(R["w"])
    return m.net(A, R, W), int((A["t"] >= W).sum())


def one(run):
    m = g1b()
    ld = pd.read_csv(f"{L}/public/{run}/{run}.lengths.tsv.gz", sep="\t", dtype={"read_id": str})
    ld = ld[ld.length > 0].reset_index(drop=True)
    lengths, idx = ld.length.to_numpy(np.int64), pd.Index(ld.read_id)
    f = pd.read_csv(f"{L}/public/{run}/{run}.primary.features.tsv.gz", sep="\t", dtype={"read_id": str})
    fa = f[f["class"] == "ARG"][["read_id", "start", "end", "scov"]]
    fb = m.mge_calls(run, "MC")
    fa, fb, _ = m.p04.drop_same_orf(fa, fb)
    rng = np.random.default_rng(m.SEED)
    r_rid, r_st, r_en = m.est.anchors(lengths, (fa.end - fa.start).to_numpy(np.int64), m.N_ANCH, rng)
    rnd = pd.DataFrame({"read_id": idx[r_rid], "start": r_st, "end": r_en, "scov": 100.0})
    rnd, _, _ = m.p04.drop_same_orf(rnd, fb)
    q = pd.concat([fa.assign(c="A"), rnd.assign(c="RND"), fb.assign(c="B")], ignore_index=True)
    d = (lengths, idx.get_indexer(q.read_id), q.start.to_numpy(np.int64), q.end.to_numpy(np.int64), q.c.to_numpy(str))
    scov = q.scov.to_numpy(float)
    glen = d[3] - d[2]
    med = float(np.median(glen[d[4] == "A"]))
    short, long_ = glen < med, glen >= med
    row = dict(run=run, A_calls=len(fa), median_A_len=med)
    row["intact"], _ = net_of(m, d, scov)
    row["intact_unweighted"], _ = net_of(m, d, scov, unweighted=True)
    row["intact_short_A"], row["sides_w_short"] = net_of(m, d, scov, a_mask=short)
    row["intact_long_A"], row["sides_w_long"] = net_of(m, d, scov, a_mask=long_)
    rate = m.p04.rate_for_n50(lengths, 3000, m.SEED)
    v = {k: [] for k in ("registered", "whole", "unweighted", "short_A", "long_A")}
    kept_long = {"registered": [], "whole": []}
    for k in range(DRAWS):
        b = m.p04.random_bounds(lengths, rate, m.SEED + k)
        dd, keep = m.est.fragment_partial(*d, b, return_keep=True)
        v["registered"].append(net_of(m, dd, scov[keep])[0])
        v["unweighted"].append(net_of(m, dd, scov[keep], unweighted=True)[0])
        v["short_A"].append(net_of(m, dd, scov[keep], a_mask=short[keep])[0])
        v["long_A"].append(net_of(m, dd, scov[keep], a_mask=long_[keep])[0])
        kept_long["registered"].append((long_ & (d[4] == "A"))[keep].sum() / max(1, (d[4] == "A")[keep].sum()))
        dw, keep_w = m.est.fragment_partial(*d, b, min_frac=1.0, return_keep=True)
        v["whole"].append(net_of(m, dw, scov[keep_w])[0])
        kept_long["whole"].append((long_ & (d[4] == "A"))[keep_w].sum() / max(1, (d[4] == "A")[keep_w].sum()))
    for k, x in v.items():
        row[f"r3000_{k}"], row[f"r3000_{k}_sd_draws"] = float(np.nanmean(x)), float(np.nanstd(x, ddof=1))
    row["share_long_A_kept_registered"] = float(np.mean(kept_long["registered"]))
    row["share_long_A_kept_whole"] = float(np.mean(kept_long["whole"]))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=8)
    a = ap.parse_args()
    with Pool(a.procs) as p:
        t = pd.DataFrame(p.map(one, RUNS, chunksize=1))
    reg = pd.concat([pd.read_csv(f"{L}/tables/P3/{d}/partB.tsv", sep="\t") for d in ("confirm", "confirm_ext")])
    reg = reg[(reg.b_class == "MC") & (reg.w == W)].drop_duplicates("run").set_index("run")
    t["registered_intact"], t["registered_r3000"] = reg.intact_net.reindex(t.run).to_numpy(), reg.r3000_net.reindex(t.run).to_numpy()
    t["reproduced"] = (np.abs(t.intact - t.registered_intact) < 1e-9) & (np.abs(t.r3000_registered - t.registered_r3000) < 1e-9)
    for k in ("registered", "whole", "unweighted"):
        base = t.intact_unweighted if k == "unweighted" else t.intact
        t[f"rel_change_{k}"] = (t[f"r3000_{k}"] - base) / base
    t.to_csv(f"{L}/tables/P3/P3_06b_cut_bias.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(t[["run", "reproduced", "intact", "r3000_registered", "r3000_whole", "rel_change_registered", "rel_change_whole",
             "rel_change_unweighted", "intact_short_A", "intact_long_A", "r3000_short_A", "r3000_long_A", "median_A_len",
             "share_long_A_kept_registered", "share_long_A_kept_whole"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
