"""P3 step 3.5, exploratory (outside the registration): the released statistic net_w on reads cut to N50 3 kb.

The registered cut test (H2) used the gene-level Part B estimator. `mobilink test` reports net_w, the unit-level
estimator of amendment 04 (units at least 1 kb from both molecule ends, weight 1 / P(span + 2 kb fits)), in which the
retention of a unit matches its weight. Here the reads of the 8 counting runs with a net linked fraction >= 0.06 are
cut as in the registered test (same breakpoint rate and seeds, 80% rule for cut calls; 5 draws) and net_w is
estimated on the pieces, each piece in the host stratum of its read. No overlaps are used: the estimates do not
depend on the clumps. The intact value must equal the registered per-run net_w.
Output: tables/P3/P3_06d_netw_cut.tsv.
usage: python -I P3_06d_netw_cut.py [--procs 8] [--draws 5]
"""
import argparse
import importlib.util
from multiprocessing import Pool

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
RUNS = ["SRR23278716", "SRR23278658", "SRR23278671", "SRR29980985", "SRR33225186", "DRR351862", "SRR16329124",
        "SRR29786052"]
DRAWS = 5


def frozen():
    spec = importlib.util.spec_from_file_location("p3_confirm", f"{L}/frozen/amendment05/P3_03_confirm.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def one(run):
    from mobilink.fragment import fragment_partial, n50, random_bounds, rate_for_n50
    from mobilink.stats.nettest import net_linked_test

    ln, f, st, _ = frozen().load_run(run)
    fa, fb = f[f["class"] == "ARG"], f[f["class"] == "MGE"]
    keys = ("net_w", "net_w_lo", "net_w_hi", "net_mol", "units", "sides_w")
    r = net_linked_test(ln, fa, fb, strata=st, loci=None, n_mc=99)
    row = dict(run=run, intact_n50=n50(ln.length.to_numpy()), **{f"intact_{k}": r[k] for k in keys})
    q = pd.concat([fa.assign(c="A"), fb.assign(c="B")], ignore_index=True)
    lengths = ln.length.to_numpy(np.int64)
    rid = pd.Index(ln.molecule_id).get_indexer(q.molecule_id)
    d = (lengths, rid, q.start.to_numpy(np.int64), q.end.to_numpy(np.int64), q.c.to_numpy(str))
    scov = q.scov.to_numpy(float)
    rate = rate_for_n50(lengths, 3000, 20261114)
    vals = []
    for k in range(DRAWS):
        off, b = random_bounds(lengths, rate, 20261114 + k)
        (pl, piece, s, e, cl), keep = fragment_partial(*d, (off, b), return_keep=True)
        parent = np.searchsorted(off, b[:-1], "right") - 1
        pln = pd.DataFrame({"molecule_id": np.arange(pl.size), "length": pl})
        pf = pd.DataFrame({"molecule_id": piece, "start": s, "end": e, "scov": scov[keep], "c": cl})
        rr = net_linked_test(pln, pf[pf.c == "A"], pf[pf.c == "B"], strata=np.asarray(st, dtype=object)[parent], loci=None,
                             n_mc=99, seed=20261114 + k)
        vals.append(dict(n50=n50(pl), **{x: rr[x] for x in keys}))
    v = pd.DataFrame(vals)
    row.update(r3000_n50=int(v.n50.mean()), r3000_net_w=v.net_w.mean(), r3000_net_w_sd_draws=v.net_w.std(),
               r3000_net_mol=v.net_mol.mean(), r3000_units=int(v.units.mean()), r3000_sides_w=int(v.sides_w.min()))
    return row


def main():
    global DRAWS
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--draws", type=int, default=5)
    a = ap.parse_args()
    DRAWS = a.draws
    with Pool(a.procs) as p:
        t = pd.DataFrame(p.map(one, RUNS, chunksize=1))
    reg = pd.concat([pd.read_csv(f"{L}/tables/P3/{d}/per_run.tsv", sep="\t") for d in ("confirm", "confirm_ext")])
    reg = reg[(reg.class_a == "ARG") & (reg.class_b == "MGE")].set_index("run")
    t["registered_net_w"] = reg.net_w.reindex(t.run).to_numpy()
    t["intact_reproduced"] = np.abs(t.intact_net_w - t.registered_net_w) < 1e-9
    t["diff"] = t.r3000_net_w - t.intact_net_w
    t["rel_change"] = t["diff"] / t.intact_net_w
    t.to_csv(f"{L}/tables/P3/P3_06d_netw_cut.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(t.round(4).to_string(index=False))
    print(t[["diff", "rel_change"]].describe().round(4).to_string())


if __name__ == "__main__":
    main()
