#!/usr/bin/env python3
"""mobileOG-db lookup table for `mobilink prep --mobileog-lookup`.

DIAMOND reports subject ids up to the first space of the FASTA header. Many mobileOG-db headers contain a space
(462,024 of 775,257 in beatrix-1.6, e.g. 'Plasmid RefSeq', 'Keyword Search'), so the evidence field is cut off in
DIAMOND output. This table restores gene, category and evidence from the full header.

usage: make_mobileog_lookup.py mobileOG-db_beatrix-1.6.All.faa mobileog_all.lookup.tsv.gz
"""
import gzip
import sys

src, out = sys.argv[1], sys.argv[2]
seen = set()
with open(src) as fh, gzip.open(out, "wt") as fo:
    fo.write("sseqid\tgene\tcategory\tevidence\n")
    for ln in fh:
        if ln[0] != ">":
            continue
        h = ln[1:].rstrip("\n")
        key = h.split()[0]
        f = h.split("|")
        if key in seen:
            sys.exit(f"duplicate DIAMOND key {key}")
        seen.add(key)
        fo.write(f"{key}\t{f[1]}\t{f[3]}\t{f[6]}\n")
print(f"{len(seen)} keys written to {out}")
