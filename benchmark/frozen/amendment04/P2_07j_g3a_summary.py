#!/usr/bin/env python3
"""P2 step 2.7: G3a summary with the pass rules of amendment 04 draft v2 (development run 6 now; G3 after freezing).

Input: per-dataset tables written by P2_07_feasibility.py (RUN/s*.tsv). Per test (net_w at 2 kb, same molecule)
and setting, under no linkage: two-sided rejection at 0.05, enrichment-side rejection (one-sided p_enrich <= 0.025),
mean estimate; under linkage (0.03, 0.1): power, coverage of the community's truth by the net_w 95% interval
(truth = unit-level net at 2 kb from the simulator). Rules (draft v2, Part 3):
  R1 size: pooled two-sided rejection within [0.025, 0.06] and every setting <= 0.08, for each test (lower bound
     revised from 0.04: net_w is conservative, 3.3-3.7% in runs 6-7, which is valid; 0.025 still catches a test
     whose variance is overstated);
  R2 enrichment side (one-sided 2.5%): pooled <= 0.030 and every setting <= 0.050, for each test (revised from
     "every setting <= 0.04", which fails an exactly calibrated test in 25% of gates with 500 datasets per setting;
     the revised rule passes it in 99% and catches a setting with a true 6% in 80%);
  R3 bias: |mean net_w| <= 0.01 and |mean net_mol| <= 0.01 under no linkage, every setting;
  R4 coverage: the net_w interval covers the truth in >= 90% of datasets, every setting, under no linkage and
     under each linked share.
Power is reported, not gated.
usage: P2_07j_g3a_summary.py RUN_DIR OUT_TSV
"""
import glob
import sys

import numpy as np
import pandas as pd

PROF = {0: "soil", 1: "wastewater", 2: "gut"}


def main():
    run, out = sys.argv[1], sys.argv[2]
    t = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{run}/s*.tsv"))], ignore_index=True)
    t["cover_w"] = (t.net_w_lo <= t.truth) & (t.truth <= t.net_w_hi)
    rows = []
    for (k, lk), d in t.groupby(["setting", "linked"]):
        row = dict(setting=k, profile=PROF[k // 6], depth=d.depth.iloc[0], hosts=d.hosts.iloc[0],
                   A="islands" if d.a_mode.iloc[0] else "independent", linked=lk, n=len(d),
                   truth=d.truth.mean(), cover_w=d.cover_w.mean(), sec=(d.sim_s + d.test_s).mean())
        for name in ("net_w", "net_mol"):
            row[f"{name}_mean"] = d[name].mean()
            row[f"{name}_reject"] = (d[f"{name}_p"] <= 0.05).mean()
            row[f"{name}_enrich"] = (d[f"{name}_p_enrich"] <= 0.025).mean()
            row[f"{name}_eff"] = d[f"{name}_eff"].median()
        rows.append(row)
    s = pd.DataFrame(rows)
    s.to_csv(out, sep="\t", index=False)
    null, lnk = s[s.linked == 0], s[s.linked > 0]
    nt = t[t.linked == 0]
    verdict = {}
    for name in ("net_w", "net_mol"):
        pooled = float((nt[f"{name}_p"] <= 0.05).mean())
        verdict[f"R1 {name}"] = (0.025 <= pooled <= 0.06 and null[f"{name}_reject"].max() <= 0.08,
                                 f"pooled {pooled:.3f}, max setting {null[f'{name}_reject'].max():.3f}")
        pe = float((nt[f"{name}_p_enrich"] <= 0.025).mean())
        verdict[f"R2 {name}"] = (pe <= 0.03 and null[f"{name}_enrich"].max() <= 0.05,
                                 f"pooled {pe:.3f}, max setting {null[f'{name}_enrich'].max():.3f}")
        verdict[f"R3 {name}"] = (null[f"{name}_mean"].abs().max() <= 0.01, f"max |mean| {null[f'{name}_mean'].abs().max():.4f}")
    verdict["R4 net_w"] = (s.cover_w.min() >= 0.90, f"min coverage {s.cover_w.min():.3f}")
    for kk, (ok, txt) in verdict.items():
        print(f"{kk}: {'pass' if ok else 'FAIL'} ({txt})")
    for lk, d in lnk.groupby("linked"):
        print(f"power at linked {lk}: net_w {d.net_w_reject.min():.2f}-{d.net_w_reject.max():.2f}, "
              f"net_mol {d.net_mol_reject.min():.2f}-{d.net_mol_reject.max():.2f}")
    print(f"seconds per dataset: {t.sim_s.add(t.test_s).mean():.1f} mean, by depth "
          f"{t.groupby('depth').apply(lambda x: (x.sim_s + x.test_s).mean()).round(1).to_dict()}")


if __name__ == "__main__":
    main()
