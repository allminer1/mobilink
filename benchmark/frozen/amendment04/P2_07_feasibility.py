#!/usr/bin/env python3
"""P2 step 2.7: feasibility of gate G3a for the net linked test (development seeds 70,000,000+; not G3 seeds).

Settings (index = 6 x profile + 2 x depth + a_mode), as the G3a draft with the host counts of draft v2:
  profile 0 soil 279, 1 wastewater SRR32022521, 2 gut SRR23278673;
  depth 0: 0.1 x with 200 hosts; 1: 5 x with 50 hosts; 2: 20 x with 50 hosts (4 Mb genomes, abundance sigma 1);
  a_mode 0: A independent (24 per genome); 1: A in islands (12 single + 3 islands per genome);
  B: 5 islands + 12 singletons; host heterogeneity density_sigma = 1; host strata = genome.
LINKED (argument): 0 (size) or 0.1 (coverage of the host-stratified truth and power).
Per dataset: linkstat.stats.nettest.net_linked_test (w = 2 kb, 200 bootstrap replicates over hosts), the truth at
2 kb (unit level, net_units, when the test uses A units; gene level net_strat otherwise), seconds.
P2_SEED0: run 1 (gene level, normal p) 70,000,000; run 2 71,000,000 (invalid: units used but t reference and unit truth missing); run 3 (A units, t on G_eff - 1, unit truth, bootstrap SE) 72,000,000; run 4 (as run 3 with the jackknife SE) 73,000,000;
  run 5 (adds the host sign-flip p and interval, expected events; P2_HOSTS01=2000 hosts at 0.1 x) 74,000,000;
  run 6 (integrated net_linked_test of D44: null-variance tests with unit-locus clumps from positions, margin-rule
  net_mol; P2_DENSITY_CAP=2.5; LINKED 0, 0.03 or 0.1) 84,000,000;
  run 7 (D45: margin rule for net_w too, Monte Carlo mid-p reference, orientation-aware same-molecule clumps;
  P2_SEED_SCHEME=2) 86,000,000.
Seed scheme 1 (runs 1-6): SEED0 + 100,000 k + 1,000,000 linked + dataset. It collides: setting k at linked 0.1 and
  setting k + 1 under no linkage share seeds (different settings, so different datasets, but shared random streams).
Seed scheme 2 (P2_SEED_SCHEME=2; run 7 and G3a): SEED0 + 1,000 k + {0: 0, 0.03: 500, 0.1: 700}[linked] + dataset,
  with dataset < 500 under no linkage and < 200 when linked: unique, SEED0 ... SEED0 + 17,899.
usage: P2_07_feasibility.py SETTING LINKED FIRST LAST OUT_TSV
"""
import sys
import time

import numpy as np
import pandas as pd

from linkstat.simulate import simulate_tables
from linkstat.stats.nettest import net_linked_test

L0 = "/mnt/rmrfc/workflow/iskander/linkstat"
PROFILES = [f"{L0}/lengths/279.lengths.tsv.gz", f"{L0}/public/SRR32022521/SRR32022521.lengths.tsv.gz",
            f"{L0}/public/SRR23278673/SRR23278673.lengths.tsv.gz"]
DEPTHS = [(0.1, int(__import__("os").environ.get("P2_HOSTS01", 200))), (5.0, 50), (20.0, 50)]  # run 5: 2,000 hosts at 0.1 x
CAP = __import__("os").environ.get("P2_DENSITY_CAP")
SEED0 = int(__import__("os").environ.get("P2_SEED0", 70_000_000))  # run 1: 70,000,000 (gene level); run 2: 71,000,000
SCHEME = int(__import__("os").environ.get("P2_SEED_SCHEME", 1))
OFFSET2 = {0.0: 0, 0.03: 500, 0.1: 700}


def seed_of(k: int, linked: float, dd: int) -> int:
    if SCHEME == 1:
        return SEED0 + k * 100_000 + int(round(linked * 1000)) * 1_000 + dd
    if dd >= (500 if linked == 0 else 200):
        raise ValueError("seed scheme 2: dataset index out of range")
    return SEED0 + k * 1_000 + OFFSET2[linked] + dd


def main():
    k, linked, first, last, out = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    prof, (depth, n_gen), a_mode = PROFILES[k // 6], DEPTHS[(k % 6) // 2], k % 2
    rl = pd.read_csv(prof, sep="\t").length.to_numpy()
    rl = rl[rl > 0]
    g = np.random.default_rng(0)
    al = np.maximum(150, g.lognormal(np.log(900), 0.35, 5000)).astype(int)
    bl = np.maximum(150, g.lognormal(np.log(1000), 0.5, 5000)).astype(int)
    a_kw = dict(n_a=24.0) if a_mode == 0 else dict(n_a=12.0, a_islands=3.0)
    rows = []
    for dd in range(first, last):
        t0 = time.time()
        s = simulate_tables(rl, al, bl, seed=seed_of(k, linked, dd),
                            n_genomes=n_gen, total_bases=depth * n_gen * 4e6, islands=5, singletons=12, clustered=True,
                            linked=linked, min_frac=0.8, n_truth_anchors=100_000, genome_length=(4_000_000, 4_000_000),
                            density_sigma=1.0, truth_stratified=True, positions=True,
                            density_cap=float(CAP) if CAP else None, **a_kw)
        t1 = time.time()
        f, ln = s.features, s.lengths
        r = net_linked_test(ln, f[f["class"] == "ARG"], f[f["class"] == "MGE"], strata=ln.genome.to_numpy(),
                            loci="positions")
        tr = s.truth.query("mode == 'end' and w == 2000").iloc[0]
        truth = tr.net_units if r.get("a_units_d") is not None else tr.net_strat
        rows.append(dict(setting=k, linked=linked, dataset=dd, seed=seed_of(k, linked, dd), depth=depth, hosts=n_gen, a_mode=a_mode,
                         molecules=len(ln), sim_s=t1 - t0, test_s=time.time() - t1, truth=truth,
                         truth_genes=tr.net_strat, truth_unstrat=tr.net,
                         **{kk: v for kk, v in r.items() if kk not in ("w",)}))
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
