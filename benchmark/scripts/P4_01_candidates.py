"""P4 step 4.1b: candidate projects of the Phase 4 sample, by rule, from archive metadata only.

Input: public/screen/ena_p4_runs_2026-11-18.tsv (ENA portal search: ONT, METAGENOMIC, >= 1 Gbp).
Eligible run: strategy WGS; selection RANDOM or unspecified; archive mean read length >= 2,000; its project has no
run that was fetched in P0-P3 (a run directory under public/, whatever its role).
Class by the organism name of the run (CLASS below). Names that do not name an environment ("metagenome",
"gut metagenome", "feces metagenome", ...) are "to read": their class comes from the project record in step 4.1c.
Output: public/screen/p4_candidates.runs.tsv (eligible runs with class_auto) and public/screen/p4_candidates.projects.tsv
(one row per project and class_auto: runs, samples, median and maximum mean read length, Gbp, title). Step 4.1c adds
class, include and reason by hand in public/screen/p4_projects.tsv; nothing here looks at reads.
usage: python -I P4_01_candidates.py
"""
import glob
import os

import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
CLASS = {
    "human gut": ["human gut metagenome", "human feces metagenome"],
    "livestock": ["bovine gut metagenome", "bovine metagenome", "rumen metagenome", "manure metagenome", "chicken gut metagenome",
                  "pig gut metagenome", "pig metagenome", "Sus scrofa", "Bos taurus", "sheep gut metagenome", "goat gut metagenome",
                  "chicken metagenome", "poultry metagenome", "livestock metagenome"],
    "wastewater": ["wastewater metagenome", "sewage metagenome", "hospital metagenome"],
    "activated sludge": ["activated sludge metagenome", "sludge metagenome", "anaerobic digester metagenome", "bioreactor metagenome",
                         "biogas fermenter metagenome"],
    "soil": ["soil metagenome", "rhizosphere metagenome", "permafrost metagenome", "peat metagenome", "compost metagenome",
             "soil crust metagenome"],
    "freshwater": ["lake water metagenome", "freshwater metagenome", "river metagenome", "aquatic metagenome", "groundwater metagenome",
                   "drinking water metagenome", "pond metagenome", "wetland metagenome"],
    "marine": ["seawater metagenome", "marine metagenome", "marine plankton metagenome", "estuary metagenome"],
    "sediment": ["sediment metagenome", "marine sediment metagenome", "freshwater sediment metagenome", "lake sediment metagenome"],
    "to read": ["metagenome", "gut metagenome", "feces metagenome", "human metagenome", "mixed culture metagenome", "biofilm metagenome",
                "landfill metagenome", "lagoon metagenome", "water metagenome", "environmental metagenome", "Homo sapiens"],
}
NAME = {n: c for c, names in CLASS.items() for n in names}


def main():
    d = pd.read_csv(f"{L}/public/screen/ena_p4_runs_2026-11-18.tsv", sep="\t", dtype=str)
    for c in ("read_count", "base_count"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["mean_len"] = d.base_count / d.read_count
    # projects with a run fetched in P0-P3
    fetched = {os.path.basename(x) for x in glob.glob(f"{L}/public/*") if os.path.isdir(x)}
    look = pd.concat([pd.read_csv(f"{L}/public/screen/{f}", sep="\t", dtype=str)[["run_accession", "study_accession"]]
                      for f in ("ena_runs.tsv", "ena_ext_runs_2026-11-17.tsv", "ena_p4_runs_2026-11-18.tsv")]).drop_duplicates()
    used = set(look.study_accession[look.run_accession.isin(fetched)])
    missing = sorted(r for r in fetched if r[:3] in ("SRR", "ERR", "DRR") and r not in set(look.run_accession))
    e = d[(d.library_strategy == "WGS") & d.library_selection.isin(["RANDOM", "unspecified"]) & (d.mean_len >= 2000)].copy()
    e["used_project"] = e.study_accession.isin(used)
    e["class_auto"] = e.scientific_name.map(NAME).fillna("out of scope")
    print(f"{len(d):,} runs; eligible by metadata {len(e):,} in {e.study_accession.nunique()} projects; "
          f"projects fetched before: {len(used)} (runs fetched before but in no screen: {len(missing)})")
    n = e[~e.used_project & (e.class_auto != "out of scope")]
    n.to_csv(f"{L}/public/screen/p4_candidates.runs.tsv", sep="\t", index=False)
    g = n.groupby(["class_auto", "study_accession"]).agg(
        runs=("run_accession", "size"), samples=("sample_accession", "nunique"), median_mean_len=("mean_len", "median"),
        max_mean_len=("mean_len", "max"), gbp=("base_count", lambda x: x.sum() / 1e9), names=("scientific_name", lambda x: "; ".join(sorted(set(x)))),
        instruments=("instrument_model", lambda x: "; ".join(sorted(set(x)))), first_public=("first_public", "min"),
        title=("study_title", "first")).reset_index().sort_values(["class_auto", "samples"], ascending=[True, False])
    g.to_csv(f"{L}/public/screen/p4_candidates.projects.tsv", sep="\t", index=False)
    s = g.groupby("class_auto").agg(projects=("study_accession", "size"), samples=("samples", "sum"),
                                    units_at_5_per_project=("samples", lambda x: int(x.clip(upper=5).sum())),
                                    projects_ge3_samples=("samples", lambda x: int((x >= 3).sum())))
    pd.set_option("display.width", 200)
    print(s.to_string())
    print("out of scope by name:", e[~e.used_project & (e.class_auto == "out of scope")].study_accession.nunique(), "projects")


if __name__ == "__main__":
    main()
