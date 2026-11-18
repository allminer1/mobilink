"""P3 step 3.2: the kit of each confirmatory run from its read ends (amendment 05 Part A, "kit and chemistry").

Input: public/<run>/<run>.adapters.tsv.gz (P2_05e in P3_02: ligation adapter 'lsk', rapid primer 'rad', native
barcodes 'nb:NBxx'; region start, end or interior) and the read count of <run>.lengths.tsv.gz.
Rule, fixed before the tables were read (2026-11-17): a pattern counts if it lies at the start of >= 5% of reads.
  native barcoding  a native barcode counts (the dominant barcode is given);
  rapid             the rapid primer counts and no native barcode does;
  ligation          the ligation adapter counts and neither of the above does;
  not determinable  nothing counts (reads trimmed before submission).
A run whose kit is not determinable from the reads keeps the kit and chemistry of its project (p3_project_kits.tsv)
when a publication or the record states one. It is replaced only if kit and chemistry are both unknown.
usage: python -I P3_02c_kits.py RUN_LIST_TSV OUT_TSV
"""
import os
import sys

import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
MIN = 0.05


def main():
    sel = pd.read_csv(sys.argv[1], sep="\t", dtype=str)
    kits = pd.read_csv(f"{L}/public/screen/p3_project_kits.tsv", sep="\t", dtype=str).set_index("project")
    rows = []
    for r in sel.itertuples():
        P = f"{L}/public/{r.run_accession}/{r.run_accession}"
        if not os.path.exists(f"{L}/public/{r.run_accession}/P3_DONE"):
            continue
        n = sum(1 for _ in pd.read_csv(f"{P}.lengths.tsv.gz", sep="\t", usecols=[0], chunksize=10**6, dtype=str)
                for _ in range(len(_)))
        a = pd.read_csv(f"{P}.adapters.tsv.gz", sep="\t", dtype={"read_id": str})
        a["kind"] = a.pattern.str.split(":").str[0]
        sh = {f"{k}_{g}": a[(a.kind == k) & (a.region == g)].read_id.nunique() / n
              for k in ("lsk", "rad", "nb") for g in ("start", "end", "interior")}
        nb = a[(a.kind == "nb") & (a.region == "start")].pattern.value_counts(normalize=True)
        if sh["nb_start"] >= MIN:
            kit = "native barcoding"
        elif sh["rad_start"] >= MIN:
            kit = "rapid"
        elif sh["lsk_start"] >= MIN:
            kit = "ligation"
        else:
            kit = "not determinable"
        k = kits.loc[r.project]
        stated_kit = "" if str(k.kit).startswith("from the read-end") else k.kit
        stated_chem = "" if str(k.chemistry).startswith("not stated") else k.chemistry
        rows.append(dict(run_accession=r.run_accession, env=r.env, project=r.project, reads=n, **sh,
                         top_barcode=nb.index[0] if len(nb) else "", top_barcode_share=float(nb.iloc[0]) if len(nb) else 0.0,
                         kit_from_reads=kit, kit_stated=stated_kit, chemistry_stated=stated_chem,
                         replace=(kit == "not determinable") and not stated_kit and not stated_chem))
    t = pd.DataFrame(rows)
    t.to_csv(sys.argv[2], sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(t[["run_accession", "project", "reads", "lsk_start", "rad_start", "nb_start", "lsk_interior", "nb_interior",
             "top_barcode", "top_barcode_share", "kit_from_reads", "replace"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
