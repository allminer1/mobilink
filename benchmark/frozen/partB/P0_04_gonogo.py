#!/usr/bin/env python3
"""P0 steps 0.3 and 0.4 on our own soils: length dependence and in-silico fragmentation.

Uses the frozen prototype (prototype/linkstat_gonogo.py, sha256 checked) for every estimator.
For each soil, threshold rule, ARG set (all, or without the intrinsic vanR_in_vanO_cl / rpoB2 / HelR) and pair:
  1. take the features of the two classes;
  2. drop same-ORF calls: an A feature and a B feature that overlap by >= 50% of the shorter are both removed
     (counted as n_dual);
  3. step 0.3: read-length bins on the intact reads; observed naive share and O/E against the exposure model
     with lambda fitted on all reads of the soil;
  4. step 0.4: prototype analyse() on intact reads, after cutting to fixed 1, 2, 3 and 5 kb pieces (prototype
     fragment()), and after random fragmentation to a target N50 of 1.5, 2 and 3 kb (fragment_random() below),
     for windows 1, 2 and 5 kb.

Why random fragmentation (added 2026-11-14): fixed cuts remove every piece longer than the cut, so a side can never
reach w once the cut is shorter than gene + w + gene (the preview on soils 279 and 280 gave NaN for w = 2 kb at all
cuts). A real library with an N50 of 3 kb keeps a tail of longer reads. Random breakpoints (a Poisson process along
each read, rate tuned to the target N50) keep such a tail. Genes that span a breakpoint are lost, as with fixed cuts.

usage: P0_04_gonogo.py --soils 278 279 ... [--rules primary strict relaxed mgeperm] [--out-prefix tables/P0]
"""
import argparse
import hashlib
import importlib.util

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
PROTO = f"{L}/prototype/linkstat_gonogo.py"
SITE = {"278p": "pilot", "278": "Ridder", "279": "Ridder", "280": "Ridder", "281": "Ridder",
        "304": "Pervomaisk", "307": "Pervomaisk", "308": "Pervomaisk", "310": "Pervomaisk"}
PAIRS = [("ARG", "MRG"), ("ARG", "MGE"), ("MRG", "MGE")]
CUTS = [1000, 2000, 3000, 5000]
RANDOM_N50 = [1500, 2000, 3000]
SEED = 20261114
WINDOWS = [1000, 2000, 5000]


def n50(x):
    s = np.sort(x)[::-1]
    return int(s[np.searchsorted(np.cumsum(s), s.sum() / 2)])


def random_bounds(lengths, rate, seed):
    """Piece boundaries in concatenated coordinates for Poisson breakpoints at `rate` per bp."""
    rng = np.random.default_rng(seed)
    off = np.concatenate(([0], np.cumsum(lengths)[:-1]))
    k = rng.poisson(lengths * rate)
    rr = np.repeat(np.arange(len(lengths)), k)
    pos = off[rr] + 1 + (rng.random(k.sum()) * (lengths[rr] - 1)).astype(np.int64)
    return off, np.unique(np.concatenate((off, pos, [off[-1] + lengths[-1]])))


def rate_for_n50(lengths, target, seed, max_reads=300_000, iters=14):
    """Breakpoint rate whose fragments have N50 closest to `target` (bisection on log rate).
    The fragment N50 depends on the read-length distribution, not on the read count, so the rate is tuned on a random
    subsample of at most `max_reads` reads (speed-up of 2026-11-14; 14 halvings of a 6.9-nat range give a step of
    0.0004 nats). The realised N50 on all reads is reported by the caller."""
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


def fragment_random(lengths, rid, start, end, cls, rate, seed, bounds=None):
    """Cut reads at Poisson breakpoints; keep genes that stay whole. Same output as prototype fragment().
    `bounds` = random_bounds(lengths, rate, seed), precomputed once per sample and target N50."""
    off, b = bounds if bounds is not None else random_bounds(lengths, rate, seed)
    gs, ge = off[rid] + start, off[rid] + end
    piece = np.searchsorted(b, gs, "right") - 1
    whole = (np.searchsorted(b, ge - 1, "right") - 1) == piece
    return (np.diff(b), piece[whole], (gs - b[piece])[whole], (ge - b[piece])[whole], cls[whole])
BINS = [0, 1000, 3000, 5000, 8000, 12000, np.inf]
BIN_LABELS = ["<1", "1-3", "3-5", "5-8", "8-12", ">12"]


def analyse_windows(g, d, a, b, label, key):
    """Same rows as prototype analyse(d, [(a, b)], w, label) for each w in WINDOWS, but read_level(), which does not
    depend on w, is computed once per version (speed-up of 2026-11-14; checked equal to analyse())."""
    lengths, rid, start, end, cls = d
    ha, hb = g.has_class(len(lengths), rid, cls, a), g.has_class(len(lengths), rid, cls, b)
    rl = g.read_level(lengths, ha, hb, g.mean_len(start, end, cls, a), g.mean_len(start, end, cls, b))
    return [dict(key, condition=label, pair=f"{a}:{b}", **rl, **g.linkage(lengths, rid, start, end, cls, a, b, w))
            for w in WINDOWS]


def load_proto():
    want = open(f"{PROTO}.sha256").read().split()[0]
    got = hashlib.sha256(open(PROTO, "rb").read()).hexdigest()
    if got != want:
        raise SystemExit(f"prototype changed: sha256 {got} != frozen {want}")
    spec = importlib.util.spec_from_file_location("gonogo", PROTO)
    g = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g)
    return g


def drop_same_orf(fa, fb):
    """Remove A and B features that overlap by >= 50% of the shorter (one ORF hit by two databases)."""
    m = fa.reset_index().merge(fb.reset_index(), on="read_id", suffixes=("_a", "_b"))
    if m.empty:
        return fa, fb, 0
    ov = np.minimum(m.end_a, m.end_b) - np.maximum(m.start_a, m.start_b)
    short = np.minimum(m.end_a - m.start_a, m.end_b - m.start_b)
    hit = m.loc[ov >= 0.5 * short]
    return fa.drop(index=hit.index_a.unique()), fb.drop(index=hit.index_b.unique()), int(hit.index_a.nunique())


def length_bins(g, lengths, rid, start, end, cls, a, b):
    """Step 0.3: per read-length bin, observed naive share vs the exposure model (lambda from all reads)."""
    n = len(lengths)
    ha, hb = g.has_class(n, rid, cls, a), g.has_class(n, rid, cls, b)
    la, lb = g.mean_len(start, end, cls, a), g.mean_len(start, end, cls, b)
    xa = np.clip(lengths - la, 0, None).astype(float); xb = np.clip(lengths - lb, 0, None).astype(float)
    pa = -np.expm1(-g.fit_lambda(xa, int(ha.sum())) * xa); pb = -np.expm1(-g.fit_lambda(xb, int(hb.sum())) * xb)
    df = pd.DataFrame(dict(bin=pd.cut(lengths, BINS, labels=BIN_LABELS, right=False), L=lengths, a=ha, b=hb,
                           both=ha & hb, pa=pa, pab=pa * pb))
    t = df.groupby("bin", observed=True).agg(n_reads=("L", "size"), mean_len=("L", "mean"), reads_A=("a", "sum"),
                                             reads_B=("b", "sum"), both=("both", "sum"), exp_A=("pa", "sum"),
                                             exp_both=("pab", "sum")).reset_index()
    t["naive_share"] = t.both / t.reads_A.replace(0, np.nan)
    t["model_share"] = t.exp_both / t.exp_A
    t["oe"] = t.both / t.exp_both
    return t


def samples(a):
    """(name, site, lengths loader, features loader) for own soils (--soils) and pooled public datasets (--public)."""
    out = []
    for s in a.soils or []:
        out.append((s, SITE.get(s, "own"),
                    lambda s=s: pd.read_csv(f"{a.lengths_dir}/{s}.lengths.tsv.gz", sep="\t", dtype={"read_id": str}),
                    lambda rule, s=s: pd.read_csv(f"{a.features_dir}/{s}.{rule}.features.tsv.gz", sep="\t",
                                                  dtype={"read_id": str})))
    for spec in a.public or []:          # NAME=RUN1,RUN2,...: the runs are concatenated (read IDs carry the run)
        name, runs = spec.split("=")
        runs = runs.split(",")
        rd = lambda r: f"{L}/public/{r}/{r}"
        out.append((name, "public",
                    lambda runs=runs: pd.concat([pd.read_csv(f"{rd(r)}.lengths.tsv.gz", sep="\t", dtype={"read_id": str})
                                                 for r in runs], ignore_index=True),
                    lambda rule, runs=runs: pd.concat([pd.read_csv(f"{rd(r)}.{rule}.features.tsv.gz", sep="\t",
                                                                   dtype={"read_id": str}) for r in runs],
                                                      ignore_index=True)))
    return out


def main():
    global CUTS, RANDOM_N50
    ap = argparse.ArgumentParser()
    ap.add_argument("--soils", nargs="+")
    ap.add_argument("--public", nargs="+", help="NAME=RUN1,RUN2 (runs pooled into one sample)")
    ap.add_argument("--rules", nargs="+", default=["primary", "strict", "relaxed", "mgeperm"])
    ap.add_argument("--cuts", nargs="+", type=int, default=CUTS)
    ap.add_argument("--random-n50", nargs="+", type=int, default=RANDOM_N50)
    ap.add_argument("--out-prefix", default=f"{L}/tables/P0")
    ap.add_argument("--lengths-dir", default=f"{L}/lengths")
    ap.add_argument("--features-dir", default=f"{L}/features")
    a = ap.parse_args()
    CUTS, RANDOM_N50 = a.cuts, a.random_n50
    g = load_proto()
    frag_rows, bin_rows = [], []
    for s, site, load_lengths, load_features in samples(a):
        ld = load_lengths()
        if not ld.read_id.is_unique:
            raise SystemExit(f"{s}: read IDs are not unique")
        lengths = ld.length.to_numpy(np.int64)
        index = pd.Index(ld.read_id)
        rates = {t: rate_for_n50(lengths, t, SEED) for t in RANDOM_N50}
        bounds = {t: random_bounds(lengths, rates[t], SEED) for t in RANDOM_N50}
        print(f"[P0_04] {s} native N50 {n50(lengths)}; breakpoint rate per kb (realised N50): "
              + ", ".join(f"target {t}: {r * 1e3:.3f} ({n50(np.diff(bounds[t][1]))})" for t, r in rates.items()),
              flush=True)
        for rule in a.rules:
            f = load_features(rule)
            for arg_set in ("all", "noint"):
                for ca, cb in PAIRS:
                    if arg_set == "noint" and "ARG" not in (ca, cb):
                        continue
                    fa, fb = f[f["class"] == ca], f[f["class"] == cb]
                    if arg_set == "noint":
                        fa = fa[~fa.intrinsic.astype(bool)] if ca == "ARG" else fa
                        fb = fb[~fb.intrinsic.astype(bool)] if cb == "ARG" else fb
                    fa, fb, n_dual = drop_same_orf(fa, fb)
                    p = pd.concat([fa, fb])
                    rid = index.get_indexer(p.read_id)
                    assert (rid >= 0).all()
                    data = (lengths, rid, p.start.to_numpy(np.int64), p.end.to_numpy(np.int64),
                            p["class"].to_numpy(str))
                    key = dict(soil=s, site=site, rule=rule, arg_set=arg_set, n_dual=n_dual)
                    if len(fa) == 0 or len(fb) == 0:
                        frag_rows.append(dict(key, pair=f"{ca}:{cb}", condition="intact", note="class empty"))
                        continue
                    t = length_bins(g, *data, ca, cb)
                    for r in t.to_dict("records"):
                        bin_rows.append(dict(key, pair=f"{ca}:{cb}", **r))
                    versions = ([("intact", data)] + [(f"cut to {c} bp", g.fragment(*data, c)) for c in CUTS]
                                + [(f"random N50 {t}", fragment_random(*data, rates[t], SEED, bounds[t])) for t in RANDOM_N50])
                    for label, d in versions:
                        frag_rows += analyse_windows(g, d, ca, cb, label, key)
            print(f"[P0_04] {s} {rule} done", flush=True)
    pd.DataFrame(frag_rows).to_csv(f"{a.out_prefix}_04_fragmentation.tsv", sep="\t", index=False)
    pd.DataFrame(bin_rows).to_csv(f"{a.out_prefix}_03_length_bins.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
