"""P3 step 3.5, exploratory: screen for reads of the ONT lambda control strand (DCS, about 3.56 kb) in the public runs
of both stages and in the 8 deep soils.

The Kraken2 Standard-16 database labels control reads Escherichia coli (taxid 562; D78), so in the registered analyses
they would sit in the E. coli stratum. Screen (no alignment; the reads of the public runs were not kept):
  public runs   among the reads labelled 562: the share with a length of 3,400-3,620 bases ("window"), and the count
  all sets      window count against the mean of the two neighbouring windows of the same width (3,180-3,400 and
                3,620-3,840), over all reads
A run is flagged if >= 100 reads labelled 562 fall in the window and they are >= 30% of the reads labelled 562, or
(soils, no Kraken2 labels for all reads) if the window holds >= 1.5 times the neighbouring mean.
Output: tables/P3/P3_06c_control_strand.tsv.
"""
import numpy as np
import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
LO, HI = 3400, 3620


def window_ratio(x):
    w = HI - LO
    c = ((x >= LO) & (x < HI)).sum()
    nb = (((x >= LO - w) & (x < LO)).sum() + ((x >= HI) & (x < HI + w)).sum()) / 2
    return int(c), float(c / nb) if nb else np.nan


def main():
    runs = [(r, "stage 1", e) for r, e in pd.read_csv(f"{L}/tables/P3/confirm/p3_final.tsv", sep="\t")[["run_accession", "env"]].itertuples(index=False)]
    runs += [(r, "stage 2", e) for r, e in pd.read_csv(f"{L}/public/screen/p3x_final.tsv", sep="\t")[["run_accession", "env"]].itertuples(index=False)]
    rows = []
    for run, stage, env in runs:
        P = f"{L}/public/{run}/{run}"
        ln = pd.read_csv(f"{P}.lengths.tsv.gz", sep="\t", dtype={"read_id": str})
        k = pd.read_csv(f"{P}.kraken.tsv.gz", sep="\t", dtype={"read_id": str}, usecols=["read_id", "taxid"])
        e = ln[ln.read_id.isin(set(k.read_id[k.taxid == 562]))].length.to_numpy()
        c, ratio = window_ratio(ln.length.to_numpy())
        ce = int(((e >= LO) & (e < HI)).sum())
        rows.append(dict(set=run, group=stage, env=env, reads=len(ln), window_reads=c, window_ratio=ratio, ecoli_reads=len(e),
                         ecoli_window_reads=ce, ecoli_window_share=ce / len(e) if len(e) else np.nan,
                         window_share_of_all_reads=ce / len(ln),
                         flag=bool(ce >= 100 and len(e) and ce / len(e) >= 0.3)))
    for s in (278, 279, 280, 281, 304, 307, 308, 310):
        x = pd.read_csv(f"{L}/lengths/{s}.lengths.tsv.gz", sep="\t", usecols=["length"]).length.to_numpy()
        c, ratio = window_ratio(x)
        rows.append(dict(set=str(s), group="deep soil", env="soil", reads=len(x), window_reads=c, window_ratio=ratio,
                         flag=bool(ratio >= 1.5)))
    t = pd.DataFrame(rows)
    t.to_csv(f"{L}/tables/P3/P3_06c_control_strand.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)
    print(t[t.flag | (t.group == "deep soil") | (t.window_ratio >= 1.3)].round(4).to_string(index=False))
    print(t.groupby("group").agg(sets=("set", "size"), flagged=("flag", "sum"), max_ratio=("window_ratio", "max")).to_string())


if __name__ == "__main__":
    main()
