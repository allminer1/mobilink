"""P3 step 3.4a, part 1: the loci of the published pilot (Ridder, Front. Microbiol. 17:1994598) as linkstat input.

Reads ridder_study/ only (nothing is written there). Catalogue Rv1 (801 ARG, 349 MRG loci; 119 metaFlye MAGs and the
four unbinned contig sets), placed on contigs through the Bakta GFF coordinates.
Bakta renamed the contigs to contig_1..n in FASTA order, for the MAGs and for the unbinned sets. The pilot's
p9_build_master.py maps the MAG names back to the assembly names but keeps the Bakta names for the unbinned sets.
Here both are mapped back (checked: the Bakta sequence length equals the assembly contig length), and the pilot's
naming is kept in the column contig_pilot so that its tables can be reproduced.
Molecules: the contigs of a sample's dereplicated MAGs and of its unbinned set (the annotated universe).
Host stratum: the MAG; unbinned contigs are "unresolved".
usage: python -I P3_05a_pilot_loci.py OUT_DIR   -> OUT_DIR/<sample>.lengths.tsv, <sample>.features.tsv, placement.tsv
"""
import csv
import glob
import os
import re
import sys

import pandas as pd

R = "/mnt/rmrfc/workflow/iskander/ridder_study"
SAMPLES = ["278", "279", "280", "281"]


def fasta_names(path):
    with open(path) as f:
        return [x[1:].split()[0] for x in f if x.startswith(">")]


def gff(path):
    """(bakta contig lengths, locus_tag -> (bakta contig, start0, end))."""
    ln, loc = {}, {}
    with open(path) as f:
        for x in f:
            if x.startswith("##sequence-region"):
                p = x.split()
                ln[p[1]] = int(p[3])
            elif "\tCDS\t" in x:
                p = x.split("\t")
                m = re.search(r"locus_tag=([^;]+)", p[8])
                if m:
                    loc[m.group(1)] = (p[0], int(p[3]) - 1, int(p[4]))
    return ln, loc


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    stat = []
    for s in SAMPLES:
        alen = pd.read_csv(f"{R}/09_colocalization/contig_lengths_{s}.tsv", sep="\t", header=None,
                           names=["contig", "length"]).set_index("contig").length
        mol, place = {}, {}                        # contig -> host ; (set, protein_id) -> (contig, start, end, pilot)
        for fa in sorted(glob.glob(f"{R}/03_binning/drep/95ani/dereplicated_genomes/{s}_*.fa")):
            mag = os.path.basename(fa)[:-3]
            b2o = {f"contig_{i + 1}": n for i, n in enumerate(fasta_names(fa))}
            g = f"{R}/05_annotation/bakta_mags/{mag}/{mag}.gff3"
            if not os.path.exists(g):
                continue
            bl, loc = gff(g)
            assert all(alen[b2o[b]] == n for b, n in bl.items()), mag
            for c in b2o.values():
                assert c not in mol, c
                mol[c] = mag
            for lt, (b, st, en) in loc.items():
                place[("mags", f"{mag}__{lt}")] = (b2o[b], st, en, b2o[b])
        b2o = {f"contig_{i + 1}": n for i, n in enumerate(fasta_names(f"{R}/05_annotation/unbinned_contigs/{s}_unbinned.fasta"))}
        bl, loc = gff(f"{R}/05_annotation/bakta_unbinned/{s}/{s}.gff3")
        assert all(alen[b2o[b]] == n for b, n in bl.items()), s
        n_both = sum(c in mol for c in b2o.values())
        for c in b2o.values():
            mol.setdefault(c, "unresolved")
        for lt, (b, st, en) in loc.items():
            place[(f"unbinned_{s}", lt)] = (b2o[b], st, en, b)      # the pilot kept the Bakta name b
        rows = []
        for path, cl in ((f"{R}/06_resistome_mobilome/integrated/arg_catalog.tsv", "ARG"),
                         (f"{R}/06_resistome_mobilome/integrated/hmrg_catalog.tsv", "MRG")):
            with open(path) as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    if r["sample"] != s:
                        continue
                    p = place.get((r["set"], r["protein_id"]))
                    gene = (r.get("rgi_aro") or r.get("amrfp_gene") or r.get("bacmet_gene")
                            or r.get("amrfp_stress_gene") or "").strip()
                    if p is None:
                        stat.append(dict(sample=s, cls=cl, protein_id=r["protein_id"], placed=False))
                        continue
                    rows.append(dict(molecule_id=p[0], start=p[1], end=p[2], **{"class": cl}, gene=gene,
                                     set=r["set"], protein_id=r["protein_id"], contig_pilot=p[3]))
                    stat.append(dict(sample=s, cls=cl, protein_id=r["protein_id"], placed=True))
        f = pd.DataFrame(rows)
        f["host"] = f.molecule_id.map(mol)
        f.to_csv(f"{out}/{s}.features.tsv", sep="\t", index=False)
        ln = pd.DataFrame({"molecule_id": list(mol), "host": list(mol.values())})
        ln["length"] = ln.molecule_id.map(alen).astype(int)
        ln[["molecule_id", "length", "host"]].to_csv(f"{out}/{s}.lengths.tsv", sep="\t", index=False)
        # the pilot's master table, by its own contig names, against the corrected placement
        m = pd.read_csv(f"{R}/09_colocalization/master_contigs_{s}.tsv", sep="\t")
        mine = f.groupby(["contig_pilot", "class"]).size().unstack(fill_value=0).reindex(columns=["ARG", "MRG"], fill_value=0)
        j = m.set_index("contig")[["n_ARG", "n_HMRG"]].join(mine, how="outer").fillna(0).astype(int)
        same = bool(((j.n_ARG == j.ARG) & (j.n_HMRG == j.MRG)).all())
        true = f.groupby(["molecule_id", "class"]).size().unstack(fill_value=0).reindex(columns=["ARG", "MRG"], fill_value=0)
        print(f"{s}: molecules {len(ln)} ({(ln.host != 'unresolved').sum()} in {ln.host.nunique() - 1} MAGs; "
              f"{n_both} in both sets), assembly contigs {len(alen)}, Mb {ln.length.sum() / 1e6:.0f} of {alen.sum() / 1e6:.0f}; "
              f"loci ARG {int((f['class'] == 'ARG').sum())} MRG {int((f['class'] == 'MRG').sum())}; "
              f"master table reproduced by pilot names: {same}; "
              f"contigs with both, pilot names {int(((j.ARG > 0) & (j.MRG > 0)).sum())} "
              f"(master {int(m.coloc.sum())}), assembly names {int(((true.ARG > 0) & (true.MRG > 0)).sum())}; "
              f"resistance contigs pilot {len(j)} true {len(true)}")
    st = pd.DataFrame(stat)
    st.groupby(["sample", "cls"]).placed.agg(["sum", "size"]).to_csv(f"{out}/placement.tsv", sep="\t")
    print(st.groupby(["sample", "cls"]).placed.agg(["sum", "size"]).to_string())


if __name__ == "__main__":
    main()
