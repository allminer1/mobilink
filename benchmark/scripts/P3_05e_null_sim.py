"""P3 step 3.4a: does a null show "depletion" when two classes are placed independently? (no real annotation used)

Contigs and MAG strata of the pilot's sample 278 (tables/P3/pilot/278.lengths.tsv: 42,662 contigs, 649 Mb). 270 A and
125 B genes of 1 kb are placed without any dependence between the classes:
  S0 uniform    every gene at a uniform position of the assembly;
  S1 clusters   A in pairs 3 kb apart, B in triplets 2.5 kb apart (gene clusters wider than the 1 kb unit rule);
                the cluster positions are uniform;
  S2 uneven     every gene on a contig drawn with probability length x weight, the weights Gamma(0.3) per contig and
                independent between the two classes (contigs differ in gene density; the classes are still independent).
Tests per simulated assembly: the pilot's label permutation (contigs with both classes; one-sided for depletion, as
published; 1,999 permutations) on loci and on units, and linkstat.pairs.net_test_pairs with the MAG strata
(same-molecule test and unit-level O/E, two-sided, both conditioning directions), without clumps (each unit its own
clump, the package default without overlaps) and with contig clumps ("cc": all A units of a contig are one clump).
usage: python -I P3_05e_null_sim.py OUT_TSV [--sims 100] [--procs 6]
"""
import argparse
import importlib.util
from multiprocessing import Pool

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
NA, NB, G = 270, 125, 1000


def place(rng, ln, cum, n, weights=None):
    if weights is None:
        c = np.searchsorted(cum, rng.random(n) * cum[-1])
    else:
        w = np.cumsum(ln * weights)
        c = np.searchsorted(w, rng.random(n) * w[-1])
    return c, (rng.random(n) * np.maximum(ln[c] - G, 1)).astype(int)


def simulate(scn, seed, lengths):
    rng = np.random.default_rng(seed)
    ln = lengths.length.to_numpy(float)
    cum = np.cumsum(ln)
    rows = []
    for cl, n, k, step in (("ARG", NA, 2, 3000), ("MRG", NB, 3, 2500)):
        if scn == "S1":
            c, s = place(rng, ln, cum, -(-n // k))
            c, s = np.repeat(c, k), np.repeat(s, k) + np.tile(np.arange(k) * step, len(c))
            ok = s + G <= ln[c]
            c, s = c[ok][:n], s[ok][:n]
        else:
            c, s = place(rng, ln, cum, n, rng.gamma(0.3, size=len(ln)) if scn == "S2" else None)
        rows.append(pd.DataFrame({"molecule_id": lengths.molecule_id.to_numpy()[c], "start": s, "end": s + G, "class": cl}))
    f = pd.concat(rows, ignore_index=True)
    f["host"] = f.molecule_id.map(lengths.set_index("molecule_id").host)
    return f


def one(args):
    scn, seed = args
    from linkstat.pairs import net_test_pairs

    spec = importlib.util.spec_from_file_location("p05b", f"{L}/scripts/P3_05b_pilot_test.py")
    pb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pb)
    pb.NPERM = 1999
    lengths = pd.read_csv(f"{L}/tables/P3/pilot/278.lengths.tsv", sep="\t")
    f = simulate(scn, seed, lengths)
    rng = np.random.default_rng(seed + 1)
    out = dict(scenario=scn, seed=seed)
    for name, slots in (("loci", f), ("units", pb.to_units(f))):
        obs = int((slots.groupby("molecule_id")["class"].nunique() > 1).sum())
        v = pb.null_c(slots, rng)
        out.update({f"perm_{name}_obs": obs, f"perm_{name}_exp": v.mean(),
                    f"perm_{name}_p_dep": (np.sum(v <= obs) + 1) / (pb.NPERM + 1),
                    f"perm_{name}_p_enr": (np.sum(v >= obs) + 1) / (pb.NPERM + 1)})
    for tag, loci in (("", None), ("cc_", pb.contig_clumps)):
        t = net_test_pairs(lengths[["molecule_id", "length"]], f[["molecule_id", "start", "end", "class"]],
                           [("ARG", "MRG"), ("MRG", "ARG")], strata=lengths.host.to_numpy(), seed=seed, loci=loci,
                           oe=loci is None)
        for r in t.itertuples():
            d = "ab" if r.class_a == "ARG" else "ba"
            out.update({f"{tag}mol_{d}_O": getattr(r, "net_mol_O", np.nan), f"{tag}mol_{d}_E": getattr(r, "net_mol_E", np.nan),
                        f"{tag}mol_{d}_p": getattr(r, "net_mol_p", np.nan), f"{tag}w_{d}_p": getattr(r, "net_w_p", np.nan)})
            if loci is None:
                out.update({f"oe_{d}_O": r.oe_n2u_observed, f"oe_{d}_E": r.oe_n2u_expected, f"oe_{d}_p": r.oe_n2u_p})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--sims", type=int, default=100)
    ap.add_argument("--procs", type=int, default=6)
    a = ap.parse_args()
    jobs = [(s, 20261117 + 1000 * i + k) for i, s in enumerate(("S0", "S1", "S2")) for k in range(a.sims)]
    with Pool(a.procs) as p:
        t = pd.DataFrame(p.map(one, jobs, chunksize=1))
    t.to_csv(a.out, sep="\t", index=False)
    rows = []
    for s, g in t.groupby("scenario"):
        r = dict(scenario=s, sims=len(g))
        for n in ("loci", "units"):
            r[f"perm_{n}_obs_over_exp"] = g[f"perm_{n}_obs"].sum() / g[f"perm_{n}_exp"].sum()
            r[f"perm_{n}_depletion_p05"] = (g[f"perm_{n}_p_dep"] <= 0.05).mean()
            r[f"perm_{n}_enrichment_p05"] = (g[f"perm_{n}_p_enr"] <= 0.05).mean()
        for d in ("ab", "ba"):
            for k in ("mol", "cc_mol", "oe"):
                r[f"{k}_{d}_obs_over_exp"] = g[f"{k}_{d}_O"].sum() / g[f"{k}_{d}_E"].sum()
                r[f"{k}_{d}_p05"] = (g[f"{k}_{d}_p"] <= 0.05).mean()
                r[f"{k}_{d}_p05_low"] = ((g[f"{k}_{d}_p"] <= 0.05) & (g[f"{k}_{d}_O"] < g[f"{k}_{d}_E"])).mean()
            for k in ("w", "cc_w"):
                r[f"{k}_{d}_p05"] = (g[f"{k}_{d}_p"] <= 0.05).mean()
        rows.append(r)
    s = pd.DataFrame(rows)
    s.to_csv(a.out.replace(".tsv", ".summary.tsv"), sep="\t", index=False)
    print(s.round(3).T.to_string())


if __name__ == "__main__":
    main()
