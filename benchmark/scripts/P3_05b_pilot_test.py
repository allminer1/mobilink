"""P3 step 3.4a, part 2: the published pilot's contig-level ARG:MRG test next to the linkstat test, on the same loci.

Input: tables/P3/pilot/<sample>.{lengths,features}.tsv (P3_05a; the pilot's catalogue Rv1 on the metaFlye contigs).
Statistic of the pilot: the number of contigs that carry both classes. Its nulls (ridder_study p19_stratified_nulls.py):
  A  each locus placed on a contig with probability proportional to contig length (full assembly);
  C  the loci stay where they are and the class labels are permuted (the paper's primary, "degree-preserving" null).
Here each null is run on three slot definitions, to see what moves the expectation:
  loci        every catalogue locus is a slot (the pilot);
  units       same-class loci within 1 kb on a contig, with no locus of the other class between them, are one slot
              (the unit rule of linkstat);
  units|host  null C only: labels permuted within the host stratum (the MAG; all unbinned contigs of a sample are one
              stratum).
Per sample, and summed over the four samples as in the paper (one site; no site comparison).
Then linkstat.pairs.net_test_pairs (molecule = contig, stratum = MAG; all A units of a contig are one clump, which the
no-linkage simulation P3_05e shows is needed on contigs),
per sample and pooled (pool_runs), with and without host strata, at w = 2 kb (registered) and 10 kb.
usage: python -I P3_05b_pilot_test.py IN_DIR OUT_PREFIX   -> OUT_PREFIX.nulls.tsv, OUT_PREFIX.linkstat.tsv
"""
import sys

import numpy as np
import pandas as pd

SAMPLES = ["278", "279", "280", "281"]
NPERM = 20000
SEED = 20261117


def to_units(f, d=1000):
    f = f.sort_values(["molecule_id", "start"], kind="stable").reset_index(drop=True)
    new = np.r_[True, (f.molecule_id.to_numpy()[1:] != f.molecule_id.to_numpy()[:-1])
                | (f["class"].to_numpy()[1:] != f["class"].to_numpy()[:-1])
                | (f.start.to_numpy()[1:] - f.end.to_numpy()[:-1] > d)]
    f["unit"] = np.cumsum(new)
    return f.groupby("unit", as_index=False).agg(molecule_id=("molecule_id", "first"), start=("start", "min"),
                                                 end=("end", "max"), host=("host", "first"),
                                                 n_loci=("start", "size"), **{"class": ("class", "first")})


def contig_clumps(units):
    """Unit-locus function for contigs: all A units of one contig are one clump (both sides merged). Units on one
    contig share the event "the contig carries B", so they are not independent draws under no linkage."""
    return pd.DataFrame({"locus": pd.factorize(units.molecule_id)[0], "parity": 0, "conflict": True})


def n_both(contig, is_a, n):
    a = np.bincount(contig, weights=is_a, minlength=n)
    t = np.bincount(contig, minlength=n)
    return int(((a > 0) & (a < t)).sum())


def null_c(slots, rng, by_host=False, a="ARG"):
    contig = pd.factorize(slots.molecule_id)[0]
    n = contig.max() + 1
    is_a = (slots["class"] == a).to_numpy(float)
    groups = [np.flatnonzero(slots.host.to_numpy() == h) for h in slots.host.unique()] if by_host else [np.arange(len(slots))]
    out = np.empty(NPERM, int)
    for k in range(NPERM):
        lab = is_a.copy()
        for g in groups:
            lab[g] = is_a[rng.permutation(g)]
        out[k] = n_both(contig, lab, n)
    return out


def null_a(slots, lengths, rng, a="ARG", b="MRG"):
    cum = np.cumsum(lengths.length.to_numpy(float))
    na, nb = int((slots["class"] == a).sum()), int((slots["class"] == b).sum())
    out = np.empty(NPERM, int)
    for k in range(NPERM):
        a = np.unique(np.searchsorted(cum, rng.random(na) * cum[-1]))
        b = np.unique(np.searchsorted(cum, rng.random(nb) * cum[-1]))
        out[k] = len(np.intersect1d(a, b, assume_unique=True))
    return out


def summary(sample, slots_name, null_name, obs, null, n_a, n_b):
    return dict(sample=sample, slots=slots_name, null=null_name, slots_A=n_a, slots_B=n_b, observed=obs,
                expected=null.mean(), lo95=np.percentile(null, 2.5), hi95=np.percentile(null, 97.5),
                ratio=obs / null.mean() if null.mean() > 0 else np.nan,
                p_depletion=(np.sum(null <= obs) + 1) / (NPERM + 1), p_enrichment=(np.sum(null >= obs) + 1) / (NPERM + 1))


def main():
    from linkstat.pairs import net_test_pairs, pool_runs

    d, out = sys.argv[1], sys.argv[2]
    rng = np.random.default_rng(SEED)
    rows, tot, runs = [], {}, []
    for s in SAMPLES:
        ln = pd.read_csv(f"{d}/{s}.lengths.tsv", sep="\t")
        f = pd.read_csv(f"{d}/{s}.features.tsv", sep="\t")
        runs.append(dict(run=s, lengths=ln[["molecule_id", "length"]], strata=ln.host.to_numpy(), paf=None,
                         features=f[["molecule_id", "start", "end", "class"]]))
        for name, slots in (("loci", f), ("units", to_units(f))):
            both = slots.groupby("molecule_id")["class"].nunique()
            obs = int((both > 1).sum())
            n_a, n_b = int((slots["class"] == "ARG").sum()), int((slots["class"] == "MRG").sum())
            nulls = {"A length-weighted placement": null_a(slots, ln, rng), "C label permutation": null_c(slots, rng)}
            if name == "units":
                nulls["C label permutation within host"] = null_c(slots, rng, by_host=True)
            for nn, v in nulls.items():
                rows.append(summary(s, name, nn, obs, v, n_a, n_b))
                t = tot.setdefault((name, nn), dict(obs=0, null=np.zeros(NPERM), a=0, b=0))
                t["obs"] += obs
                t["null"] += v
                t["a"] += n_a
                t["b"] += n_b
        print(f"[3.4a] nulls {s} done", flush=True)
    for (name, nn), t in tot.items():
        rows.append(summary("Ridder (sum of 4)", name, nn, t["obs"], t["null"], t["a"], t["b"]))
    pd.DataFrame(rows).to_csv(f"{out}.nulls.tsv", sep="\t", index=False)

    res = []
    pooled = pool_runs(runs)
    data = [(r["run"], r["lengths"], r["features"], r["strata"]) for r in runs] + [
        ("Ridder (pooled)", pooled[0], pooled[1], pooled[2])]
    for name, ln, f, st in data:
        for w in (2000, 10000):
            for strata_name, stv in (("MAG", st), ("none", None)):
                t = net_test_pairs(ln, f, [("ARG", "MRG"), ("MRG", "ARG")], strata=stv, w=w, oe=(w == 2000),
                                   loci=contig_clumps)
                t.insert(0, "strata", strata_name)
                if "w" not in t:
                    t.insert(0, "w", w)
                t.insert(0, "sample", name)
                res.append(t)
        print(f"[3.4a] linkstat {name} done", flush=True)
    pd.concat(res, ignore_index=True).to_csv(f"{out}.linkstat.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
