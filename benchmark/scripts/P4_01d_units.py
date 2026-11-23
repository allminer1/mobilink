"""P4 step 4.1d: runs, units and the selected sample of Phase 4, from archive metadata only.

Inputs: public/screen/p4_projects.tsv (hand-entered project table), public/screen/p4_candidates.runs.tsv (eligible runs)
and the run records in public/screen/p4_records/.
Per run: class from the project's run_rule (or the project's class), unit from its unit_rule (or the BioSample).
Per unit: the run with the longest archive mean read length represents the unit (ties by accession).
Per project and class: at most 5 units, spread over read length without a random element: units ranked by mean read
length, and the ranks round(i (n - 1) / 4), i = 0..4, are taken.
Class rule: a class is kept if at least 6 projects and 20 selected units remain.
Outputs: public/screen/p4_runs_classed.tsv (every eligible run of the included projects, with class, unit and the text
the rules were applied to), public/screen/p4_units.tsv (one row per unit, selected yes or no), public/screen/p4_selected.tsv
(the selected runs of the kept classes) and a summary on stdout.
usage: python -I P4_01d_units.py [--show PROJECT]
"""
import io
import re
import sys

import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
MAX_UNITS = 5


def main():
    show = sys.argv[2] if len(sys.argv) > 2 and sys.argv[1] == "--show" else None
    pt = pd.read_csv(f"{L}/public/screen/p4_projects.tsv", sep="\t", dtype=str).fillna("")
    pt = pt[pt.include == "yes"].set_index("project")
    runs = pd.read_csv(f"{L}/public/screen/p4_candidates.runs.tsv", sep="\t", dtype=str)
    runs["mean_len"] = pd.to_numeric(runs.mean_len)
    runs["base_count"] = pd.to_numeric(runs.base_count)
    runs = runs[runs.study_accession.isin(pt.index)].drop_duplicates("run_accession")
    out = []
    for p, g in runs.groupby("study_accession"):
        rec = pd.read_csv(io.StringIO(open(f"{L}/public/screen/p4_records/{p}.runs.tsv").read()), sep="\t", dtype=str).fillna("")
        g = g.merge(rec, on="run_accession", how="left", suffixes=("", "_rec")).fillna("")
        cols = [c for c in ("experiment_title", "sample_title", "isolation_source", "host", "sample_alias", "sample_description") if c in g]
        g["text"] = g[cols].agg(" | ".join, axis=1).str.lower()
        rule, urule = pt.run_rule[p], pt.unit_rule[p]
        if rule:
            items = [x.split("=>") for x in rule.split(" ;; ")]
            g["class"] = [next((c for rx, c in items if re.search(rx, t)), "") for t in g.text]
        else:
            g["class"] = pt["class"][p]
        if urule:
            g["unit"] = [(m.group(1) if (m := re.search(urule, t)) else "") for t in g.text]
            bad = g[(g["class"] != "") & (g.unit == "")]
            if len(bad):
                print(f"!! {p}: unit rule matches nothing in {len(bad)} runs, e.g. {bad.text.iloc[0][:140]}")
            g["unit"] = [u if u else s for u, s in zip(g.unit, g.sample_accession)]
        else:
            g["unit"] = g.sample_accession
        out.append(g[["run_accession", "study_accession", "class", "unit", "sample_accession", "mean_len", "base_count", "instrument_model", "first_public", "text"]])
    R = pd.concat(out, ignore_index=True)
    R.to_csv(f"{L}/public/screen/p4_runs_classed.tsv", sep="\t", index=False)
    if show:
        pd.set_option("display.width", 250, "display.max_colwidth", 150, "display.max_rows", 200)
        print(R[R.study_accession == show][["run_accession", "class", "unit", "mean_len", "text"]].to_string(index=False))
        return
    K = R[R["class"] != ""].sort_values(["mean_len", "run_accession"], ascending=[False, True])
    U = K.groupby(["study_accession", "class", "unit"], sort=False).head(1).copy()
    U["runs_of_unit"] = U.set_index(["study_accession", "class", "unit"]).index.map(K.groupby(["study_accession", "class", "unit"]).size())
    U["selected"] = "no"
    for (p, c), g in U.groupby(["study_accession", "class"]):
        g = g.sort_values(["mean_len", "run_accession"], ascending=[False, True])
        n = len(g)
        take = sorted({round(i * (n - 1) / (MAX_UNITS - 1)) for i in range(MAX_UNITS)}) if n > MAX_UNITS else range(n)
        U.loc[g.index[list(take)], "selected"] = "yes"
    U = U.sort_values(["class", "study_accession", "mean_len"], ascending=[True, True, False])
    U.drop(columns="text").to_csv(f"{L}/public/screen/p4_units.tsv", sep="\t", index=False)
    S = U[U.selected == "yes"]
    # class rule (plan 4.1): a class is kept if at least 6 projects and 20 selected units remain
    k = S.groupby("class").agg(projects=("study_accession", "nunique"), units=("unit", "size"))
    kept = sorted(k.index[(k.projects >= 6) & (k.units >= 20)])
    S.drop(columns="text")[S["class"].isin(kept)].to_csv(f"{L}/public/screen/p4_selected.tsv", sep="\t", index=False)
    print("classes kept:", kept, "; dropped:", sorted(set(k.index) - set(kept)), "; runs in the sample:", int(S["class"].isin(kept).sum()))
    pd.set_option("display.width", 220)
    print(pd.DataFrame({"projects": U.groupby("class").study_accession.nunique(), "units": U.groupby("class").size(),
                        "selected": S.groupby("class").size(), "median_mean_len": S.groupby("class").mean_len.median().round(0),
                        "selected_ge_5kb": S[S.mean_len >= 5000].groupby("class").size()}).fillna(0).astype(int).to_string())
    print("selected runs:", len(S), "in", S.study_accession.nunique(), "projects; Gbp to fetch at 3 Gbp per run:",
          round(S.base_count.clip(upper=3.3e9).sum() / 1e9))
    dropped = R[R["class"] == ""].groupby("study_accession").size()
    print("runs left out by a run rule:", dropped.to_dict())
    none = sorted(set(pt.index) - set(K.study_accession))
    print("included projects without a classed run:", none)


if __name__ == "__main__":
    main()
