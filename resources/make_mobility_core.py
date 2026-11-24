#!/usr/bin/env python3
"""The MGE "mobility core": mobileOG-db entries that move DNA within or between cells.

The rule uses only the annotation fields of mobileOG-db (beatrix-1.6), never the hits of a dataset. An entry is in
the core if all of:
  1. its major category is integration/excision or transfer;
  2. its evidence is Manual or Homology (not Keyword Search);
  3. its minor category names none of: regulation, partitioning, competence, replication/recombination/repair
     (regulators, plasmid maintenance, natural transformation and housekeeping recombination have chromosomal
     homologs);
  4. its gene family (lower-case name without a trailing _N or -N) is not in the housekeeping list H below.

usage: make_mobility_core.py mobileOG-db-beatrix-1.6-All.csv mobility_core.tsv.gz
The shipped src/mobilink/data/mobility_core.tsv.gz (109,325 entries) is the output for beatrix-1.6.
"""
import sys

import pandas as pd

MAJOR = {"integration/excision", "transfer"}
EVIDENCE = {"Manual", "Homology"}
MINOR_EXCLUDE = ("regulation", "partitioning", "competence", "replication/recombination/repair")
H = sorted({
    "xerc", "xerd",                                                              # chromosome-dimer resolution
    "reca", "recx", "recb", "recc", "recd", "recf", "recj", "recn", "reco", "recq", "recr",
    "ruva", "ruvb", "ruvc", "uvra", "uvrb", "uvrc", "uvrd", "muts", "mutl", "nth", "rada",   # recombination, repair
    "dnaa", "dnab", "dnag", "pola", "liga", "gyra", "gyrb", "topa", "ssb",          # replication
    "para", "parb", "ftsk", "smc",                                                 # segregation, division
    "clpp", "clpx",                                                                # proteases
    "copr", "ompr", "phob",                                                        # response regulators
})


def main():
    src, out = sys.argv[1], sys.argv[2]
    c = pd.read_csv(src, dtype=str, keep_default_na=False,
                    usecols=["mobileOG Entry Name", "Name", "Major mobileOG Category", "Minor mobileOG Categories",
                             "Evidence", "Database", "mobileOG fasta Header"])
    p = pd.DataFrame({"sseqid": c["mobileOG fasta Header"].str.split(" ").str[0], "entry": c["mobileOG Entry Name"],
                      "gene": c["Name"], "major": c["Major mobileOG Category"],
                      "minor": c["Minor mobileOG Categories"].replace({"NA": ""}),
                      "evidence": c["Evidence"].replace({"Keyword Search": "Keyword"}), "source": c["Database"]})
    p["family"] = p.gene.str.lower().str.replace(r"[_\-]\d+$", "", regex=True)
    keep = (p.major.isin(MAJOR) & p.evidence.isin(EVIDENCE)
            & ~p.minor.fillna("").str.contains("|".join(MINOR_EXCLUDE), regex=True) & ~p.family.isin(H))
    mc = p[keep]
    mc[["sseqid", "entry", "gene", "family", "major", "minor", "evidence", "source"]].to_csv(out, sep="\t", index=False)
    print(f"{len(p):,} entries, {len(mc):,} in the mobility core -> {out}")


if __name__ == "__main__":
    main()
