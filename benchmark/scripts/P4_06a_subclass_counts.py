"""P4 step 4.6, part a: counts of the registered ARG and MGE calls by sub-class, per environment class.

No linkage statistic. ARG calls get the resistance mechanism and the drug class of the CARD index (CARD 4.0.1, the
release of the search database); MGE calls get the major mobileOG category of the mobility core table.

usage: python -I P4_06a_subclass_counts.py OUT_DIR [--procs 12]     (OUT_DIR: tables/P4/sub)
"""
import argparse
from multiprocessing import Pool

import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
CARD = "/mnt/rmrfc/workflow/iskander/kmaarg/db/card/aro_index.tsv"


def card_table():
    c = pd.read_csv(CARD, sep="\t", dtype=str)
    c = c.rename(columns={"CARD Short Name": "gene", "Drug Class": "drug", "Resistance Mechanism": "mech", "AMR Gene Family": "family"})
    return c[["gene", "drug", "mech", "family"]].drop_duplicates("gene").set_index("gene")


def one(run):
    f = pd.read_csv(f"{L}/public/{run}/{run}.core.features.tsv.gz", sep="\t", dtype={"read_id": str}, usecols=["class", "gene"])
    a = f[f["class"] == "ARG"].gene.value_counts().rename("n").reset_index().assign(cls="ARG")
    m = f[f["class"] == "MGE"].gene.str.split("|").str[3].value_counts().rename("n").reset_index().assign(cls="MGE")
    return pd.concat([a, m]).assign(run=run)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=12)
    a = ap.parse_args()
    sel = pd.read_csv(f"{L}/tables/P4/confirm/p4_final.tsv", sep="\t", dtype=str)[["run_accession", "env", "project"]]
    with Pool(a.procs) as p:
        t = pd.concat(p.map(one, sel.run_accession, chunksize=4)).merge(sel, left_on="run", right_on="run_accession").drop(columns="run_accession")
    c = card_table()
    arg = t[t.cls == "ARG"].join(c, on="gene")
    print("ARG calls", int(arg.n.sum()), "| without a CARD index row", int(arg.n[arg.mech.isna()].sum()), sorted(arg.gene[arg.mech.isna()].unique())[:20], flush=True)
    arg.to_csv(f"{a.out}/arg_calls_by_gene.tsv", sep="\t", index=False)
    t[t.cls == "MGE"].to_csv(f"{a.out}/mge_calls_by_category.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
    print(t[t.cls == "MGE"].pivot_table(index="gene", columns="env", values="n", aggfunc="sum", fill_value=0).to_string())
    print(arg.pivot_table(index="mech", columns="env", values="n", aggfunc="sum", fill_value=0, margins=True).to_string())
    d = arg.pivot_table(index="drug", columns="env", values="n", aggfunc="sum", fill_value=0, margins=True).sort_values("All", ascending=False)
    print(d.head(30).to_string())
    print(arg.groupby("gene").n.sum().sort_values(ascending=False).head(25).to_string())


if __name__ == "__main__":
    main()
