"""P4 step 4.1c, first part: archive records of the candidate projects, for reading (metadata only).

For each project of public/screen/p4_candidates.projects.tsv: the ENA study record (title, description) and, for its
eligible runs, the experiment and sample fields that say how the library was made and what was sampled.
Output: public/screen/p4_records/<project>.study.tsv and <project>.runs.tsv; a digest for reading in
public/screen/p4_records.digest.tsv (one row per project).
usage: python -I P4_01b_records.py
"""
import io
import os
import time
import urllib.parse
import urllib.request

import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
O = f"{L}/public/screen/p4_records"
API = "https://www.ebi.ac.uk/ena/portal/api/search"
RUN_FIELDS = ("run_accession,experiment_title,library_name,library_construction_protocol,instrument_model,sample_alias,"
              "sample_description,host,isolation_source,environment_biome,country,collection_date")


def get(params, tries=3):
    url = API + "?" + urllib.parse.urlencode(params)
    for k in range(tries):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "metadata-screen"}), timeout=90).read().decode()
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(5 * (k + 1))
    return f"ERROR\t{err}"


def table(txt):
    return pd.read_csv(io.StringIO(txt), sep="\t", dtype=str) if txt.strip() and not txt.startswith("ERROR") else pd.DataFrame()


def main():
    os.makedirs(O, exist_ok=True)
    c = pd.read_csv(f"{L}/public/screen/p4_candidates.projects.tsv", sep="\t", dtype=str)
    runs = pd.read_csv(f"{L}/public/screen/p4_candidates.runs.tsv", sep="\t", dtype=str)
    rows = []
    for p in c.study_accession.unique():
        fs, fr = f"{O}/{p}.study.tsv", f"{O}/{p}.runs.tsv"
        if not os.path.exists(fs):
            open(fs, "w").write(get(dict(result="study", query=f'study_accession="{p}"', fields="study_accession,study_title,study_description,center_name", format="tsv", limit=0)))
            time.sleep(0.4)
        if not os.path.exists(fr):
            open(fr, "w").write(get(dict(result="read_run", query=f'study_accession="{p}" AND instrument_platform="OXFORD_NANOPORE"', fields=RUN_FIELDS, format="tsv", limit=0)))
            time.sleep(0.4)
        s, r = table(open(fs).read()), table(open(fr).read())
        if len(r):
            r = r[r.run_accession.isin(runs.run_accession)]
        first = lambda col: "; ".join(sorted(set(r[col].dropna().astype(str))))[:300] if len(r) and col in r else ""  # noqa: E731
        rows.append(dict(project=p, title=(s.study_title.iloc[0] if len(s) else ""), description=(str(s.study_description.iloc[0])[:900] if len(s) else ""),
                         center=(s.center_name.iloc[0] if len(s) else ""), eligible_runs=len(r),
                         protocol=first("library_construction_protocol"), experiment_titles=first("experiment_title"),
                         hosts=first("host"), isolation_source=first("isolation_source"), biome=first("environment_biome"),
                         country=first("country"), sample_descriptions=first("sample_description")))
    d = pd.DataFrame(rows).merge(c.groupby("study_accession").agg(class_auto=("class_auto", lambda x: "; ".join(sorted(set(x)))), samples=("samples", lambda x: sum(map(int, x))),
                                                                  names=("names", "first"), instruments=("instruments", "first")).reset_index(),
                                 left_on="project", right_on="study_accession").drop(columns="study_accession")
    d.to_csv(f"{L}/public/screen/p4_records.digest.tsv", sep="\t", index=False)
    print(len(d), "projects;", int((d.title == "").sum()), "without a study record;", int((d.eligible_runs == 0).sum()), "without run records")


if __name__ == "__main__":
    main()
