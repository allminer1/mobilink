#!/usr/bin/env python3
"""P1 step 1.6: the G1b estimator and gate columns (amendment 02 Part A3; this script is Part B's implementation).

Per run (one run = one unit):
  A  ARG calls, primary rule (CARD protein homolog, identity >= 80%, subject cover >= 80%), all ARGs
  B  MGE calls of the mobility core (MC, db/mobility_core.tsv.gz): mobileOG hits with identity >= 60% and subject cover
     >= 80%, one call per ORF as in P0 (P0_02c collapse over all mobileOG hits), kept if the call's subject is in MC.
     --b-class D3 uses the registered P0 class instead (secondary).
  Same-ORF rule between A and B as registered (>= 50% of the shorter).
Net linked fraction at w (default 2 kb):
  sides     for each A side, event time = distance to the far end of the nearest B; censored at the read end;
            a B call cut by a read end (ending exactly at a read or fragment end, or within 50 bp of it with subject
            cover < 95%) censors the side instead of ending it
  F(w)      Kaplan-Meier, A sides weighted by 1 / P(a gene of the observed length lies whole on a read of this read set)
  chance    the same weighted estimate for 20,000 random anchors (lengths drawn from the A calls, placed uniformly on the
            reads where they fit, anchors overlapping a B call by >= 50% dropped), fragmented together with A and B
  net(w)    1 - (1 - F_A) / (1 - F_anchor); NaN if fewer than 30 A sides reach w
Versions: intact, and random fragmentation to N50 3 kb and 2 kb (P0 breakpoint rates; 20 draws, seeds 20261114 + 0..19);
  a call cut by a breakpoint is kept on the piece holding >= 80% of it (the caller's subject-cover rule), truncated.
Intervals: 95% percentile bootstrap over reads (A sides and anchor sides resampled by read), 500 replicates. The 3 kb
  (and 2 kb) value is the mean over the 20 draws; its half-width is sqrt(h^2 + (1.96 sd)^2), h = the bootstrap half-width
  at seed 20261114 and sd = the between-draw sd.
Gate columns (A3): counts (N50 >= 9 kb, >= 30 sides reach w in both versions, 3 kb half-width <= 0.10), naive share
  ratio intact / 3 kb (mean over draws) >= 1.5, |net intact - net 3 kb| <= max(0.02, 0.25 |net intact|), intervals overlap.
usage: P1_06_g1b.py --runs RUN [RUN ...] --out FILE [--b-class MC|D3] [--draws 20] [--boot 500]
"""
import argparse
import importlib.util

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


p04 = load("p04", f"{L}/frozen/partB/P0_04_gonogo.py")
p02c = load("p02c", f"{L}/frozen/partB/P0_02c_features.py")
est = load("est", f"{L}/frozen/partB/P1_03_estimators.py")
SEED, N_ANCH, TARGETS = 20261114, 20000, (3000, 2000)


def mge_calls(run, b_class):
    """MGE calls (one per ORF over all mobileOG hits at 60/80), restricted to the mobility core or to D3."""
    d = p02c.collapse(p02c.read_hits(f"{L}/public/{run}/{run}.mobileog.all.tsv.gz", 60.0))
    if b_class == "MC":
        mc = set(pd.read_csv(f"{L}/frozen/partB/mobility_core.tsv.gz", sep="\t", usecols=["sseqid"]).sseqid)
        d = d[d.sseqid.isin(mc)]
    else:
        lk = pd.read_csv(f"{L}/db/mobileog_all.lookup.tsv.gz", sep="\t", dtype=str, usecols=["sseqid", "category"])
        d = d.merge(lk, on="sseqid", how="left")
        d = d[d.category.isin(p02c.MGE_PRIMARY)]
    return d.rename(columns={"qseqid": "read_id"})[["read_id", "start", "end", "scov"]]


def side_table(ln, rid, st, en, cl, cut, who):
    i, t, e = est.sides(ln, rid, st, en, cl == who, cl == "B", "end", cut=cut)
    return dict(t=t, e=e, w=1 / est.p_whole(ln, (en - st)[i]), read=rid[i])


def net(A, R, W, ia=slice(None), ir=slice(None)):
    if np.sum(A["t"][ia] >= W) < 30:
        return np.nan
    fa = est.km(A["t"][ia], A["e"][ia], W, A["w"][ia], min_at_risk=1)
    fr = est.km(R["t"][ir], R["e"][ir], W, R["w"][ir], min_at_risk=1)
    return 1 - (1 - fa) / (1 - fr)


def boot(A, R, W, B, rng):
    ga, gr = pd.factorize(A["read"])[0], pd.factorize(R["read"])[0]
    by_a = np.split(np.argsort(ga, kind="stable"), np.cumsum(np.bincount(ga))[:-1])
    by_r = np.split(np.argsort(gr, kind="stable"), np.cumsum(np.bincount(gr))[:-1])
    out = np.empty(B)
    for b in range(B):
        ia = np.concatenate([by_a[k] for k in rng.integers(0, len(by_a), len(by_a))])
        ir = np.concatenate([by_r[k] for k in rng.integers(0, len(by_r), len(by_r))])
        out[b] = net(A, R, W, ia, ir)
    return np.nanquantile(out, [0.025, 0.975])


def analyse(ln, rid, st, en, cl, scov, W):
    cut = est.cut_mask(ln, rid, st, en, scov)
    A, R = side_table(ln, rid, st, en, cl, cut, "A"), side_table(ln, rid, st, en, cl, cut, "RND")
    hasA = np.zeros(len(ln), bool); hasA[rid[cl == "A"]] = True
    hasB = np.zeros(len(ln), bool); hasB[rid[cl == "B"]] = True
    naive = (hasA & hasB).sum() / hasA.sum() if hasA.sum() else np.nan
    return A, R, dict(n50=p04.n50(ln), reads_A=int(hasA.sum()), naive=naive, sides=int(A["t"].size),
                      sides_w=int((A["t"] >= W).sum()), net=net(A, R, W))


def run_one(run, b_class, W, draws, B):
    ld = pd.read_csv(f"{L}/public/{run}/{run}.lengths.tsv.gz", sep="\t", dtype={"read_id": str})
    ld = ld[ld.length > 0].reset_index(drop=True)
    lengths = ld.length.to_numpy(np.int64); idx = pd.Index(ld.read_id)
    f = pd.read_csv(f"{L}/public/{run}/{run}.primary.features.tsv.gz", sep="\t", dtype={"read_id": str})
    fa = f[f["class"] == "ARG"][["read_id", "start", "end", "scov"]]
    fb = mge_calls(run, b_class)
    fa, fb, n_dual = p04.drop_same_orf(fa, fb)
    rng = np.random.default_rng(SEED)
    r_rid, r_st, r_en = est.anchors(lengths, (fa.end - fa.start).to_numpy(np.int64), N_ANCH, rng)
    rnd = pd.DataFrame({"read_id": idx[r_rid], "start": r_st, "end": r_en, "scov": 100.0})
    rnd, _, _ = p04.drop_same_orf(rnd, fb)
    q = pd.concat([fa.assign(c="A"), rnd.assign(c="RND"), fb.assign(c="B")], ignore_index=True)
    d = (lengths, idx.get_indexer(q.read_id), q.start.to_numpy(np.int64), q.end.to_numpy(np.int64), q.c.to_numpy(str))
    scov = q.scov.to_numpy(float)
    row = dict(run=run, b_class=b_class, w=W, n_dual=n_dual, A_calls=len(fa), B_calls=len(fb))
    A, R, s = analyse(*d, scov, W)
    lo, hi = boot(A, R, W, B, np.random.default_rng(SEED + 1))
    row.update({f"intact_{k}": v for k, v in s.items()}, intact_lo=lo, intact_hi=hi)
    for t in TARGETS:
        rate = p04.rate_for_n50(lengths, t, SEED)
        vals = []
        for k in range(draws):
            (dd, keep) = est.fragment_partial(*d, p04.random_bounds(lengths, rate, SEED + k), return_keep=True)
            A, R, s = analyse(*dd, scov[keep], W)
            vals.append(s)
            if k == 0:
                lo, hi = boot(A, R, W, B, np.random.default_rng(SEED + 2))
                h0 = (hi - lo) / 2
        v = pd.DataFrame(vals)
        sd = v.net.std()
        half = float(np.sqrt(h0 ** 2 + (1.96 * sd) ** 2))
        row.update({f"r{t}_n50": int(v.n50.mean()), f"r{t}_naive": v.naive.mean(), f"r{t}_sides_w": int(v.sides_w.min()),
                    f"r{t}_net": v.net.mean(), f"r{t}_net_sd_draws": sd, f"r{t}_half": half,
                    f"r{t}_lo": v.net.mean() - half, f"r{t}_hi": v.net.mean() + half})
    n, r3 = row["intact_net"], row["r3000_net"]
    row["counts"] = bool(row["intact_n50"] >= 9000 and row["intact_sides_w"] >= 30 and row["r3000_sides_w"] >= 30
                         and row["r3000_half"] <= 0.10)
    row["naive_ratio"] = row["intact_naive"] / row["r3000_naive"] if row["r3000_naive"] > 0 else np.nan
    row["naive_ok"] = bool(row["naive_ratio"] >= 1.5)
    row["diff"] = abs(n - r3); row["tol"] = max(0.02, 0.25 * abs(n))
    row["stable"] = bool(row["diff"] <= row["tol"])
    row["overlap"] = bool(row["intact_lo"] <= row["r3000_hi"] and row["r3000_lo"] <= row["intact_hi"])
    row["meets"] = bool(row["counts"] and row["naive_ok"] and row["stable"] and row["overlap"])
    row["rel_change_3kb"] = (r3 - n) / n if abs(n) >= 0.10 else np.nan
    row["rel_change_2kb"] = (row["r2000_net"] - n) / n if abs(n) >= 0.10 else np.nan
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--b-class", default="MC", choices=["MC", "D3"])
    ap.add_argument("--w", type=int, default=2000)
    ap.add_argument("--draws", type=int, default=20)
    ap.add_argument("--boot", type=int, default=500)
    a = ap.parse_args()
    rows = []
    for r in a.runs:
        rows.append(run_one(r, a.b_class, a.w, a.draws, a.boot))
        pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
        x = rows[-1]
        print(f"[1.6] {r} {a.b_class}: N50 {x['intact_n50']}, net {x['intact_net']:.4f} [{x['intact_lo']:.4f}, "
              f"{x['intact_hi']:.4f}] -> 3 kb {x['r3000_net']:.4f} [{x['r3000_lo']:.4f}, {x['r3000_hi']:.4f}], "
              f"2 kb {x['r2000_net']:.4f}; naive ratio {x['naive_ratio']:.2f}; meets {x['meets']}", flush=True)


if __name__ == "__main__":
    main()
