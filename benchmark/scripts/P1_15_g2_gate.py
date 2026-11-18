#!/usr/bin/env python3
"""P1 step 1.15: the G2 verdict (amendment 03 v2, Part 2).

G2a (--g2a, P1_15_g2a.py output): cases with E(N2-ev) >= 10 count; passes if >= 75% of them have |ratio - 1| <= 0.05
and the median |ratio - 1| <= 0.05.
G2b (--g2b, P1_15_g2b.py outputs): per setting, datasets with E >= 5 count; passes if in every setting the two-sided
rejection rate of the conditional p value at 0.05 lies within [0.03, 0.06] and the one-sided (enrichment) rate is
<= 0.06. Reported: mid-p rates, rates without strata, of the registered length-only null and of N2-ev, the share of
loci on >= 2 reads.
usage: P1_15_g2_gate.py --g2a FILE --g2b FILE [FILE ...] --out PREFIX
"""
import argparse

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g2a", required=True)
    ap.add_argument("--g2b", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    t = pd.read_csv(a.g2a, sep="\t")
    c = t[t.counts.fillna(False).astype(bool)]
    dev = (c.ratio - 1).abs()
    g2a = dict(cases=len(t), counting=len(c), within_5pct=int((dev <= 0.05).sum()),
               share_within=float((dev <= 0.05).mean()) if len(c) else np.nan,
               median_abs_dev=float(dev.median()) if len(c) else np.nan)
    g2a["pass"] = bool(len(c) > 0 and g2a["share_within"] >= 0.75 and g2a["median_abs_dev"] <= 0.05)

    b = pd.concat([pd.read_csv(f, sep="\t") for f in a.g2b], ignore_index=True)
    rows = []
    for k, d in b.groupby("setting"):
        n = d[d.E >= 5]
        rows.append(dict(setting=int(k), datasets=len(d), counting=len(n), mean_E=n.E.mean(), mean_oe=n.oe.mean(),
                         reject_two=(n.p_two <= 0.05).mean(), reject_enr=(n.p_enr <= 0.05).mean(),
                         reject_two_mid=(n.p_two_mid <= 0.05).mean(), reject_enr_mid=(n.p_enr_mid <= 0.05).mean(),
                         nostrata_oe=n.nostrata_oe.mean(), nostrata_reject_two=(n.nostrata_p_two <= 0.05).mean(),
                         registered_oe=n.registered_oe.mean(), registered_reject_two=(n.registered_p_two <= 0.05).mean(),
                         n2ev_oe=n.n2ev_oe.mean(), n2ev_reject_two=(n.n2ev_p_two <= 0.05).mean(),
                         share_loci_2plus=n.share_loci_2plus.mean(),
                         conditioned_A=int((n.conditioned == "A").sum())))
    s = pd.DataFrame(rows)
    s["pass"] = (s.reject_two >= 0.03) & (s.reject_two <= 0.06) & (s.reject_enr <= 0.06)
    g2b_pass = bool(s["pass"].all() and len(s) == 6)
    s.to_csv(f"{a.out}_g2b_settings.tsv", sep="\t", index=False)
    verdict = dict(**{f"g2a_{k}": v for k, v in g2a.items()}, g2b_pass=g2b_pass,
                   g2b_settings_passing=int(s["pass"].sum()), g2_pass=bool(g2a["pass"] and g2b_pass))
    pd.DataFrame([verdict]).to_csv(f"{a.out}_verdict.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(s.round(4).to_string(index=False))
    for k, v in verdict.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
