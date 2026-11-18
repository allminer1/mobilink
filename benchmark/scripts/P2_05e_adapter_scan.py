"""P2 step 2.5: find ONT adapter and barcode sequences anywhere in reads (chimera flag "mid-read adapter").

Each read is scanned in 500 bp windows every 400 bp. In each window edlib (HW mode) reports the best hit of each
pattern within its edit-distance limit. Hits of one pattern within 50 bp are merged. A hit is 'start' if it begins
in the first 150 bp, 'end' if it ends in the last 150 bp, else 'interior'.

Patterns (Dorado sequences; tables/P2/mocks/ont_barcodes_dorado.tsv), both orientations:
  lsk   ligation adapter core CCTGTACTTCGTTCAGTTACGTATTGC, k = 5         (all kits with the LSK114 adapter)
  rad   rapid primer (50 bp), k = 10                                     (rapid and rapid-barcoding kits)
  nb    AAGGTTAA + NBxx + CAGCACCT (40 bp), k = 6, for the barcodes given (native barcoding)
Random 150 bp sequence never matched lsk (k = 5) or rad (k = 10) in 20,000 tries.

usage: python -I P2_05e_adapter_scan.py READS_FA BARCODES_TSV OUT_TSV_GZ PROCS [NBxx,NByy,...|all|none]
output: read_id  read_len  pattern  strand  start  end  ed  region
"""
import gzip
import multiprocessing as mp
import sys

import edlib

LSK = "CCTGTACTTCGTTCAGTTACGTATTGC"
RAD = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
WIN, STEP, EDGE, MERGE = 500, 400, 150, 50
COMP = str.maketrans("ACGTN", "TGCAN")
PATTERNS = []


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    """FASTA or FASTQ (plain or gzip); the read id is the first word of the header."""
    name, seq = None, []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        first = f.read(1)
        f.seek(0)
        if first == "@":
            for h in f:
                s = next(f).strip().upper()
                next(f)
                next(f)
                yield h[1:].split()[0], s
            return
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip().upper())
    if name is not None:
        yield name, "".join(seq)


def scan(rec):
    name, s = rec
    L = len(s)
    hits = []
    for pname, strand, q, k in PATTERNS:
        found = []
        o = 0
        while o < L:
            seg = s[o:o + WIN]
            if len(seg) >= len(q) // 2:
                r = edlib.align(q, seg, mode="HW", task="locations", k=k)
                if r["editDistance"] >= 0:
                    a, b = r["locations"][0]
                    found.append((o + a, o + b + 1, r["editDistance"]))
            if o + WIN >= L:
                break
            o += STEP
        found.sort()
        merged = []
        for a, b, d in found:
            if merged and a - merged[-1][0] <= MERGE:
                if d < merged[-1][2]:
                    merged[-1] = (a, b, d)
                continue
            merged.append((a, b, d))
        for a, b, d in merged:
            region = "start" if a < EDGE else "end" if b > L - EDGE else "interior"
            hits.append(f"{name}\t{L}\t{pname}\t{strand}\t{a}\t{b}\t{d}\t{region}\n")
    return hits


def init(patterns):
    PATTERNS.extend(patterns)


def main():
    fa, bc_tsv, out, procs = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    which = sys.argv[5] if len(sys.argv) > 5 else "none"
    pats = [("lsk", "+", LSK, 5), ("lsk", "-", rc(LSK), 5), ("rad", "+", RAD, 10), ("rad", "-", rc(RAD), 10)]
    if which != "none":
        nb = dict(line.split() for line in open(bc_tsv) if line.startswith("NB"))
        names = sorted(n for n in nb if int(n[2:]) <= 24) if which == "all" else which.split(",")
        for n in names:
            q = "AAGGTTAA" + nb[n] + "CAGCACCT"
            pats += [(f"nb:{n}", "+", q, 6), (f"nb:{n}", "-", rc(q), 6)]
    n_reads = n_hit = 0
    with gzip.open(out, "wt") as fo, mp.Pool(procs, initializer=init, initargs=(pats,)) as pool:
        fo.write("read_id\tread_len\tpattern\tstrand\tstart\tend\ted\tregion\n")
        for hits in pool.imap(scan, read_fasta(fa), chunksize=200):
            n_reads += 1
            if hits:
                n_hit += 1
                fo.writelines(hits)
    print(f"[2.5] adapter scan: {n_reads} reads, {n_hit} with a hit, {len(pats)} patterns")


if __name__ == "__main__":
    main()
