"""P3 extension (amendment 06, Part A): the candidate runs of stage 2 and the order in which they are tried.

Input: public/screen/ena_ext_runs_2026-11-17.tsv (ENA portal search of 2026-11-17: ONT, library source METAGENOMIC,
>= 3 Gbp; run metadata only) and public/screen/p3x_projects.tsv (the project table of Part A: class, include, reason).
Eligible run: strategy WGS; selection RANDOM or unspecified; >= 3 Gbp; archive mean read length >= 6,000 bp; organism
name of one of the three classes; project not used in P0-P3 (any role); project included in the project table; for
PRJNA1398240 only the runs whose sample alias ends in "-raw-long" (the others are enrichments).
Unit: the BioSample, or the individual or site where the sample title names it (UNIT_RULES). The run with the longest
archive mean read length represents its unit.
Order (no random element): within a project, units by archive mean read length, longest first; within a class, projects
by the mean read length of their first unit, longest first. Try order = round 1 (first unit of every project in project
order), round 2 (second units), and so on. Ties by accession.
Output: OUT_PREFIX.eligible.tsv (all eligible runs) and OUT_PREFIX.order.tsv (one row per unit, with class, project,
project_rank, unit_rank, try_order).
usage: python -I P3X_01_candidates.py OUT_PREFIX
"""
import re
import sys

import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
NAMES = {
    "human gut": ["human gut metagenome", "human feces metagenome", "gut metagenome", "feces metagenome"],
    "livestock": ["bovine gut metagenome", "bovine metagenome", "rumen metagenome", "manure metagenome",
                  "chicken gut metagenome", "pig gut metagenome", "pig metagenome"],
    "wastewater": ["wastewater metagenome", "activated sludge metagenome", "sludge metagenome"],
}
UNIT_RULES = {  # project: regex on the sample title; group 1 is the unit
    "PRJEB88320": r"^(D\d+)_T\d+$",          # donor, two time points
    "PRJNA1134722": r"^(.*?)(?:\s+\d{4}-\d{2}-\d{2}.*)?$",  # treatment plant station, several dates
}


def main():
    out = sys.argv[1]
    d = pd.read_csv(f"{L}/public/screen/ena_ext_runs_2026-11-17.tsv", sep="\t", dtype=str)
    for c in ("read_count", "base_count"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["mean_len"] = d.base_count / d.read_count
    proj = pd.read_csv(f"{L}/public/screen/p3x_projects.tsv", sep="\t", dtype=str)
    inc = proj[proj.include == "yes"].set_index("project")["class"]
    used = pd.read_csv(f"{L}/public/screen/ena_runs.tsv", sep="\t", dtype=str)
    used = set(used.study_accession) | set(used.project)
    assert not (set(inc.index) & used), "an included project was used before"
    name2class = {n: c for c, ns in NAMES.items() for n in ns}
    e = d[(d.library_strategy == "WGS") & d.library_selection.isin(["RANDOM", "unspecified"]) & (d.base_count >= 3e9)
          & (d.mean_len >= 6000) & d.scientific_name.isin(name2class) & d.study_accession.isin(inc.index)].copy()
    e["class"] = e.study_accession.map(inc)
    assert (e.scientific_name.map(name2class) == e["class"]).all(), "organism name and project class disagree"
    alias = {}
    if (e.study_accession == "PRJNA1398240").any():  # sample alias from the frozen table of this project
        a = pd.read_csv(f"{L}/public/screen/p3x_PRJNA1398240_alias.tsv", sep="\t", dtype=str)
        alias = dict(zip(a.run_accession, a.sample_alias))
        keep = (e.study_accession != "PRJNA1398240") | e.run_accession.map(alias).fillna("").str.endswith("-raw-long")
        e = e[keep]

    def unit(r):
        rule = UNIT_RULES.get(r.study_accession)
        m = re.match(rule, str(r.sample_title)) if rule else None
        return f"{r.study_accession}:{m.group(1).strip()}" if m else f"{r.study_accession}:{r.sample_accession}"

    e["unit"] = [unit(r) for r in e.itertuples()]
    e = e.sort_values(["mean_len", "run_accession"], ascending=[False, True])
    cols = ["class", "study_accession", "run_accession", "unit", "sample_accession", "sample_title", "instrument_model",
            "read_count", "base_count", "mean_len", "first_public"]
    e[cols].rename(columns={"study_accession": "project"}).to_csv(f"{out}.eligible.tsv", sep="\t", index=False)
    u = e.drop_duplicates("unit").copy()
    u["unit_rank"] = u.groupby("study_accession").cumcount() + 1
    first = u[u.unit_rank == 1].sort_values(["class", "mean_len", "study_accession"], ascending=[True, False, True])
    first["project_rank"] = first.groupby("class").cumcount() + 1
    u = u.merge(first[["study_accession", "project_rank"]], on="study_accession")
    u = u.sort_values(["class", "unit_rank", "project_rank"]).reset_index(drop=True)
    u["try_order"] = u.groupby("class").cumcount() + 1
    u[cols + ["project_rank", "unit_rank", "try_order"]].rename(columns={"study_accession": "project"}).to_csv(
        f"{out}.order.tsv", sep="\t", index=False)
    g = u.groupby(["class", "study_accession"]).agg(units=("unit", "size"), longest=("mean_len", "max"),
                                                    project_rank=("project_rank", "first")).reset_index()
    print(g.sort_values(["class", "project_rank"]).round(0).to_string(index=False))
    print(u.groupby("class").agg(projects=("study_accession", "nunique"), units=("unit", "size")).to_string())


if __name__ == "__main__":
    main()
