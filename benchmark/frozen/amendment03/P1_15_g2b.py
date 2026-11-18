#!/usr/bin/env python3
"""P1 step 1.15, gate G2b (amendment 03 v2, Part 2): type I error of the read-level test N0-ev under no linkage.

Generator: linkstat.simulate; no linkage; B in islands at the median mobility-core density of the G1b runs
(islands = 5, singletons = 12, island geometry as in the simulator); whole and 80% partial calls; depth 0.1 x the
community size (n_genomes x 4 Mb); strata = the source genome of each read.
Settings (index = 3 x density + profile):
  profile 0 soil 279 lengths, 1 wastewater SRR32022521, 2 gut SRR23278673;
  density 0: A 24 per 4 Mb, 200 genomes; density 1: A 6 per 4 Mb, 800 genomes.
Datasets: seeds 20261116 + setting x 100,000 + dataset. Only datasets with E >= 5 count (E of the gated test).
Per dataset: gated N0-ev with genome strata (conditional p two-sided and one-sided enrichment, mid-p); reported:
N0-ev without strata, the registered length-only N0 with genome strata, and N2-ev with genome strata.
usage: P1_15_g2b.py SETTING FIRST LAST OUT_TSV     (datasets FIRST..LAST-1)
"""
import sys

import numpy as np
import pandas as pd

from linkstat.nulls.analytic import cooccurrence_test, cooccurrence_test_ev
from linkstat.nulls.binned import binned_test
from linkstat.simulate import simulate_tables

L0 = "/mnt/rmrfc/workflow/iskander/linkstat"
PROFILES = [f"{L0}/lengths/279.lengths.tsv.gz", f"{L0}/public/SRR32022521/SRR32022521.lengths.tsv.gz",
            f"{L0}/public/SRR23278673/SRR23278673.lengths.tsv.gz"]
DENSITIES = [(24.0, 200), (6.0, 800)]


def setting(k):
    return PROFILES[k % 3], DENSITIES[k // 3]


def main():
    k, first, last, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    prof, (n_a, n_gen) = setting(k)
    rl = pd.read_csv(prof, sep="\t").length.to_numpy()
    rl = rl[rl > 0]
    g = np.random.default_rng(0)  # gene-length arrays, fixed for all settings
    al = np.maximum(150, g.lognormal(np.log(900), 0.35, 5000)).astype(int)
    bl = np.maximum(150, g.lognormal(np.log(1000), 0.5, 5000)).astype(int)
    rows = []
    for d in range(first, last):
        sim = simulate_tables(rl, al, bl, seed=20261116 + k * 100_000 + d, n_genomes=n_gen,
                              total_bases=0.1 * n_gen * 4e6, n_a=n_a, islands=5, singletons=12, clustered=True,
                              linked=0.0, min_frac=0.8, n_truth_anchors=1000)
        f, ln = sim.features, sim.lengths
        fa, fb = f[f["class"] == "ARG"], f[f["class"] == "MGE"]
        strata = ln.genome.to_numpy()
        r, cond_on, x = cooccurrence_test_ev(ln, fa, fb, strata=strata, tost_bound=None)
        u, _, _ = cooccurrence_test_ev(ln, fa, fb, strata=None, tost_bound=None)
        ids, L = ln.molecule_id, ln.length.to_numpy()
        ha, hb = ids.isin(fa.molecule_id).to_numpy(), ids.isin(fb.molecule_id).to_numpy()
        reg = cooccurrence_test(L, ha, hb, (fa.end - fa.start).to_numpy(float), (fb.end - fb.start).to_numpy(float),
                                strata=strata, tost_bound=None)
        hx, hy = (hb, ha) if cond_on == "A" else (ha, hb)  # (class whose exposure is x, conditioned class)
        n2 = binned_test(L, hx, hy, strata=strata, bin_by=x)
        rpl = f.groupby("locus").molecule_id.nunique()
        rows.append(dict(setting=k, dataset=d, conditioned=cond_on, molecules=len(ln), O=r.observed, E=r.expected,
                         oe=r.oe, p_two=r.p_two_sided, p_enr=r.p_enriched, p_two_mid=r.p_two_sided_mid,
                         p_enr_mid=r.p_enriched_mid, nostrata_oe=u.oe, nostrata_p_two=u.p_two_sided,
                         registered_oe=reg.oe, registered_p_two=reg.p_two_sided, n2ev_E=n2.expected, n2ev_oe=n2.oe,
                         n2ev_p_two=n2.p_two_sided, share_loci_2plus=float((rpl >= 2).mean())))
        if (d - first + 1) % 100 == 0:
            pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
            print(f"setting {k}: {d - first + 1} datasets", flush=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
