"""P3 step 3.1 (amendment 05, Part A): the confirmatory public sample, chosen from archive metadata only.

Input: public/screen/ena_runs.tsv (ENA portal, fetched 2026-11-14; sha256 in the amendment). No read is accessed.
Eligible run: ONT; library source METAGENOMIC; strategy WGS (or OTHER for PRJNA1414918, as in amendment 02 A2.1);
selection RANDOM or unspecified; >= 1 Gbp; mean read length (bases / reads) >= 1,500 bp; not a mock, reference or
synthetic sample; not used in P0-P2 (the 9 P0 runs, the 8 G1b runs, the mock runs) and not from a unit used there.
Unit (independence key): wastewater site; gut patient; pig farm; manure sample; chicken; soil plot (PRJEB86780, title
without the time point) or sample accession.
Per unit, the run with the most bases represents it (ties: accession order).
Per project with quota k: the units are sorted by mean read length and cut into k consecutive groups of (nearly) equal
size; one unit is drawn from each group (numpy default_rng seeded with crc32 of the project accession), so the sample
spans the project's read lengths. The remaining units, in a seeded random order, are the replacement list.
Quotas (12 per environment class, 48 runs): wastewater PRJNA1207805 12; gut PRJNA929328 6, PRJEB108131 6; livestock
PRJNA1279599 4, PRJNA1414918 4, ERP183232 4; soil PRJEB58634 4, PRJEB86780 4, PRJNA1400213 4.
usage: python -I P3_01_select.py OUT_PREFIX      -> OUT_PREFIX.eligible.tsv, .selected.tsv, .replacements.tsv
"""
import re
import sys
import zlib

import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
QUOTA = {"PRJNA1207805": ("wastewater", 12), "PRJNA929328": ("human gut", 6), "PRJEB108131": ("human gut", 6),
         "PRJNA1279599": ("livestock", 4), "PRJNA1414918": ("livestock", 4), "ERP183232": ("livestock", 4),
         "PRJEB58634": ("soil", 4), "PRJEB86780": ("soil", 4), "PRJNA1400213": ("soil", 4)}
OTHER_OK = {"PRJNA1414918"}
USED_RUNS = {"SRR32022307", "SRR32022410", "SRR32022392", "SRR23278727", "SRR23278735", "SRR23278724", "SRR34068995",
             "SRR34068994", "ERR17169941",                                              # P0
             "SRR34068997", "SRR34068996", "SRR34068990", "SRR23278673", "SRR23278681", "SRR23278721", "SRR32022521",
             "SRR32022302"}                                                             # G1b / P2
MOCK = re.compile(r"zymo|mock|standard|reference|synthetic", re.I)


def main():
    out = sys.argv[1]
    r = pd.read_csv(f"{L}/public/screen/ena_runs.tsv", sep="\t", low_memory=False)
    g = pd.read_csv(f"{L}/tmp/prjna929328_samples.tsv", sep="\t")[["run_accession", "sample_alias"]]
    r = r.merge(g, on="run_accession", how="left")
    x = r[r.project.isin(QUOTA) & (r.instrument_platform == "OXFORD_NANOPORE")].copy()
    x["env"] = x.project.map(lambda p: QUOTA[p][0])
    x["mean_len"] = x.base_count / x.read_count
    title = x.sample_title.fillna("").astype(str)
    unit = x.sample_accession.astype(str)
    m = x.project == "PRJNA929328"
    unit[m] = "patient " + x.sample_alias[m].astype(str).str.split("_").str[0]
    m = x.project == "PRJNA1207805"
    unit[m] = "site " + title[m].str.extract(r"^[A-Za-z]+_([A-Za-z]+)\d+", expand=False).astype(str)
    m = x.project == "PRJEB108131"
    unit[m] = "patient " + title[m].str.extract(r"^(CAR-\d+)", expand=False).astype(str)
    m = x.project == "PRJEB86780"
    unit[m] = "plot " + title[m].str.replace(r"_T\d+$", "", regex=True)
    m = x.project.isin(["PRJNA1414918", "PRJNA1279599", "ERP183232"])
    unit[m] = title[m]
    x["unit"] = x.project + ": " + unit
    used_units = set(x.loc[x.run_accession.isin(USED_RUNS), "unit"]) | {
        "PRJNA929328: patient 15", "PRJNA929328: patient 16", "PRJNA929328: patient 17",
        "PRJNA1207805: site CW", "PRJNA1207805: site I", "PRJNA1207805: site TP"}
    ok = ((x.library_source == "METAGENOMIC")
          & ((x.library_strategy == "WGS") | ((x.library_strategy == "OTHER") & x.project.isin(OTHER_OK)))
          & x.library_selection.isin(["RANDOM", "unspecified"])
          & (x.base_count >= 1e9) & (x.mean_len >= 1500)
          & ~title.str.contains(MOCK) & ~x.run_accession.isin(USED_RUNS) & ~x.unit.isin(used_units)
          & ~unit.str.endswith("nan"))
    e = x[ok].sort_values(["project", "unit", "base_count", "run_accession"], ascending=[True, True, False, True])
    cols = ["env", "project", "run_accession", "unit", "sample_title", "instrument_model", "read_count", "base_count",
            "mean_len"]
    e[cols].to_csv(f"{out}.eligible.tsv", sep="\t", index=False)
    rep = e.drop_duplicates("unit")                      # the largest run of each unit
    sel, repl = [], []
    for p, (env, k) in QUOTA.items():
        d = rep[rep.project == p].sort_values(["mean_len", "run_accession"]).reset_index(drop=True)
        rng = np.random.default_rng(zlib.crc32(p.encode()))
        groups = np.array_split(np.arange(len(d)), k) if len(d) >= k else [np.array([i]) for i in range(len(d))]
        pick = [int(rng.choice(gr)) for gr in groups if len(gr)]
        s = d.iloc[pick].assign(length_group=range(1, len(pick) + 1))
        sel.append(s)
        rest = d.drop(index=pick)
        rest = rest.iloc[rng.permutation(len(rest))].assign(order=range(1, len(rest) + 1))
        repl.append(rest)
        print(f"{env:11s} {p:13s} eligible runs {int((e.project == p).sum()):4d}  units {len(d):4d}  quota {k:2d}  "
              f"selected {len(s):2d}  mean length {d.mean_len.min():.0f}-{d.mean_len.max():.0f}")
    pd.concat(sel)[cols + ["length_group"]].to_csv(f"{out}.selected.tsv", sep="\t", index=False)
    pd.concat(repl)[cols + ["order"]].to_csv(f"{out}.replacements.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
