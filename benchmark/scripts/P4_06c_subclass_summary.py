"""P4 step 4.6, part c (descriptive): summary tables from the sub-class linkage tables of P4_06b.

  subclass.mge_difference.tsv  per class: net linked fraction of all ARG to integration/excision genes minus that to
                               transfer genes, with the delete-one-project jackknife of the difference
  subclass.composition.tsv     per class: share of the ARG units in each mechanism and in each drug family
  subclass.influence.tsv       per sub-class and class: range of the delete-one-project estimates and the project
                               whose removal moves the estimate most
Also writes manuscript/supplement/Table_S11_subclasses.tsv (the class table with the influence columns).
Nothing here is a registered test.

usage: python -I P4_06c_subclass_summary.py TABLE_DIR      (TABLE_DIR: tables/P4/sub)
"""
import importlib.util
import sys

import numpy as np
import pandas as pd
from scipy import stats

L = "/mnt/rmrfc/workflow/iskander/linkstat"


def frozen_p4():
    s = importlib.util.spec_from_file_location("p4_03_confirm", f"{L}/frozen/amendment07/P4_03_confirm.py")
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def main():
    out = sys.argv[1]
    fz = frozen_p4()
    d = pd.read_csv(f"{out}/subclass.classes.tsv", sep="\t")
    envs = sorted(d.env.unique())

    def tab(sub, e):
        return pd.read_csv(f"{out}/linkage/{sub}/{e.replace(' ', '_')}.tsv", sep="\t").set_index("left_out")

    rows = []
    for e in envs:
        a, b = tab("mge_intexc", e), tab("mge_transfer", e)
        full = a.net_w["none"] - b.net_w["none"]
        reps = (a.net_w - b.net_w).drop("none").to_numpy(float)
        se, k = fz.jack_se(full, reps)
        t = stats.t.ppf(0.975, k - 1) if k >= 2 else np.nan
        rows.append(dict(env=e, net_w_intexc=a.net_w["none"], net_w_transfer=b.net_w["none"], difference=full, se_jack=se,
                         lo95=full - t * se, hi95=full + t * se, jackknife_reps=k, reps_positive=int((reps > 0).sum())))
    pd.DataFrame(rows).to_csv(f"{out}/subclass.mge_difference.tsv", sep="\t", index=False)

    comp = []
    for pre in ("mech_", "fam_"):
        u = d[d["sub"].str.startswith(pre)].pivot(index="sub", columns="env", values="units").fillna(0)
        comp.append((u / u.sum()).reset_index().melt(id_vars="sub", var_name="env", value_name="share_of_arg_units"))
    pd.concat(comp).to_csv(f"{out}/subclass.composition.tsv", sep="\t", index=False)

    inf = []
    for r in d.itertuples():
        g = tab(r.sub, r.env)
        if "net_w" not in g or not np.isfinite(r.net_w):
            continue
        rep = g.net_w.drop("none").dropna()
        if len(rep) < 2:
            continue
        i = (rep - r.net_w).abs().idxmax()
        inf.append(dict(sub=r.sub, env=r.env, loo_min=rep.min(), loo_max=rep.max(), most_influential_project=i, net_w_without_it=rep[i]))
    inf = pd.DataFrame(inf)
    inf.to_csv(f"{out}/subclass.influence.tsv", sep="\t", index=False)
    s11 = d.drop(columns="registered_all").merge(inf, on=["sub", "env"], how="left")
    s11.to_csv(f"{L}/manuscript/supplement/Table_S11_subclasses.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).round(4).to_string(index=False))
    print("cells", len(d), "| with an estimate", int(np.isfinite(d.net_w).sum()), "| of these under 100 sides", int((np.isfinite(d.net_w) & (d.sides_w < 100)).sum()))


if __name__ == "__main__":
    main()
