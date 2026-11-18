#!/usr/bin/env python3
"""P1 step 1.15, gate G2a (amendment 03 v2, Part 2): analytic vs permutation expectation, with the excluded-volume
correction, on the 8 G1b runs (native reads, primary rule).

Features: `linkstat prep` with CARD, BacMet2 EXP and mobileOG (class MGE = mobility core), primary rule, written to
tmp/g2/<run>.primary.features.tsv.gz by scripts/P1_15_g2_driver.sh. Pairs ARG:MGE and ARG:MRG, after the same-ORF rule.
Per run and pair: E(N0-ev) (nulls.analytic.cooccurrence_test_ev) and E(N2-ev) (nulls.binned.binned_test with 1 kb bins
of the same exposure, roles as in N0-ev); ratio = E(N0-ev) / E(N2-ev); a case counts if E(N2-ev) >= 10.
Reported, not gated: the registered length-only N0 and N2, observed counts and O/E.
usage: P1_15_g2a.py OUT_TSV
"""
import sys

import pandas as pd

from linkstat.features import drop_same_orf
from linkstat.io import read_lengths
from linkstat.nulls.analytic import cooccurrence_test, cooccurrence_test_ev
from linkstat.nulls.binned import binned_test

L0 = "/mnt/rmrfc/workflow/iskander/linkstat"
PAIRS = [("ARG", "MGE"), ("ARG", "MRG")]


def main():
    out = sys.argv[1]
    runs = pd.read_csv(f"{L0}/public/screen/g1b_selected.tsv", sep="\t").run_accession
    rows = []
    for run in runs:
        ln = read_lengths(f"{L0}/public/{run}/{run}.lengths.tsv.gz")
        f = pd.read_csv(f"{L0}/tmp/g2/{run}.primary.features.tsv.gz", sep="\t", dtype={"molecule_id": str})
        L, ids = ln.length.to_numpy(), ln.molecule_id
        for a, b in PAIRS:
            fa, fb, n_dual = drop_same_orf(f[f["class"] == a], f[f["class"] == b])
            if fa.empty or fb.empty:
                rows.append(dict(run=run, pair=f"{a}:{b}", note="class empty"))
                continue
            r, cond_on, x = cooccurrence_test_ev(ln, fa, fb, tost_bound=None)
            ha, hb = ids.isin(fa.molecule_id).to_numpy(), ids.isin(fb.molecule_id).to_numpy()
            hx, hy = (hb, ha) if cond_on == "A" else (ha, hb)
            n2 = binned_test(L, hx, hy, bin_by=x)
            reg = cooccurrence_test(L, ha, hb, (fa.end - fa.start).to_numpy(float),
                                    (fb.end - fb.start).to_numpy(float), tost_bound=None)
            reg2 = binned_test(L, ha, hb)
            rows.append(dict(run=run, pair=f"{a}:{b}", n_dual=n_dual, conditioned=b if cond_on == "B" else a,
                             molecules_A=int(ha.sum()), molecules_B=int(hb.sum()), observed=r.observed,
                             E_N0ev=r.expected, E_N2ev=n2.expected, ratio=r.expected / n2.expected
                             if n2.expected > 0 else float("nan"), counts=n2.expected >= 10,
                             oe_N0ev=r.oe, oe_N2ev=n2.oe, E_N0_registered=reg.expected, oe_N0_registered=reg.oe,
                             E_N2_registered=reg2.expected, oe_N2_registered=reg2.oe))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
