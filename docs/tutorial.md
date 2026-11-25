# Tutorial: one public sample from download to report

This page runs `mobilink` on a public nanopore metagenome and explains the report that comes out. Every command was
run as written on 2026-11-25, with the databases downloaded anew, and the numbers below are from that run. Use it to
check your installation: if your numbers match, the installation works.

| | |
|---|---|
| Sample | SRR33541077: human gut (neonate), nanopore R10.4.1, project PRJNA1260441 |
| Reads used | The first 150,000 reads: 403,468,987 bases, read N50 4.9 kb |
| Disk | 16 GB for the Kraken2 index, and 12 GB more for its archive until that is deleted; 0.4 GB for the gene databases, 0.4 GB for the reads, 45 MB of results |
| Memory | 16 GB (the Kraken2 index is held in memory) |
| Time of the analysis | 7 minutes with 8 threads (reads); 6 minutes for the optional assembly and 1 minute for its analysis |
| Time of the downloads | Depends on your connection. Ours took about 50 minutes for the databases and 60 for the reads |

## 1. Install

```bash
mamba create -n mobilink -c conda-forge -c bioconda mobilink diamond minimap2 kraken2 sra-tools
mamba activate mobilink
```

`sra-tools` is needed only to fetch the example reads.

`mobilink db` and `mobilink run` belong to version 1.1. Until that version is released, the conda package does not
have them; install the tools as above without `mobilink` and add the package from its development branch:

```bash
mamba create -n mobilink -c conda-forge -c bioconda python=3.12 pip diamond minimap2 kraken2 sra-tools
mamba activate mobilink
pip install git+https://github.com/allminer1/mobilink@wp3-run
```

## 2. Gene databases

```bash
mobilink db --out mobilink_db --threads 8
```

This fetches CARD, BacMet2, mobileOG-db and NCBI UniVec, checks them and builds the DIAMOND indexes. The last lines
should read:

```text
card: version 4.0.2, 6,059 sequences. ...
bacmet: version 2.0 EXP, 753 sequences. ...
mobileog: version beatrix-1.6, 775,257 sequences. ...
database directory ready: mobilink_db
```

- CARD is fetched as its current release, so your version can be later than 4.0.2. The other two are fixed versions.
- The mobileOG-db file (355 MB) comes from Zenodo and can be slow. A download that is cut is taken up again where it
  stopped. If it fails for good, the message names the file to download by hand and the option to give it with.

## 3. Kraken2 index for host labels

```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20241228.tar.gz
mkdir k2_standard_16gb && tar -xzf k2_standard_16gb_20241228.tar.gz -C k2_standard_16gb
rm k2_standard_16gb_20241228.tar.gz
```

This is the Standard-16 index of December 2024, the one used in the paper. Other indexes are listed on the
[Kraken2 index page](https://benlangmead.github.io/aws-indexes/k2). The directory must hold `ktaxonomy.tsv`.

## 4. Reads

```bash
fastq-dump -X 150000 --stdout SRR33541077 | gzip -1 > SRR33541077.fastq.gz
zcat SRR33541077.fastq.gz | awk 'NR%4==2{n++; b+=length($0)} END{print n, b}'
```

The second command must print `150000 403468987`. If it prints fewer reads, the download was cut: run it again.

## 5. Run

```bash
mobilink run SRR33541077.fastq.gz --db mobilink_db --kraken-db k2_standard_16gb --out results --threads 8
```

The nine steps are printed as they start (lengths, three DIAMOND searches, features, host labels, overlaps, chimera
flags, test, descriptive tables, report). At the end the result table and the notes are printed, and the report is
at `results/SRR33541077.report.html`. A run that is stopped goes on where it stopped when the command is given again.

## 6. What you should see

The first lines of the report:

> ARG and MGE are linked: 4.8% of the flanks of ARG units have a whole MGE gene within 2 kb, beyond what chance gives
> (95% interval -0.039 to 0.136); p = 0.032. The test and the interval use different methods; the interval is the
> more cautious one and includes 0 here, so treat the size as uncertain.

The result table (`results/SRR33541077.test.tsv`; the report shows the same numbers rounded):

| Pair | Net linked within 2 kb | 95% interval | p | Units | Observed / expected linked units |
|---|---|---|---|---|---|
| ARG–MGE | 0.048 | −0.039 to 0.136 | 0.032 | 1,738 | 132 / 34.1 |
| MRG–MGE | 0.037 | −0.343 to 0.416 | 0.054 | 2,400 | 110 / 47.0 |
| ARG–MRG | 0.005 | −0.019 to 0.028 | 0.980 | 1,787 | 16 / 19.6 |

The input box: 150,000 reads, 0.40 Gbp, N50 4.9 kb; 2,154 ARG calls, 3,289 MRG calls and 6,775 MGE calls; 4% of the
reads without a host label; 227 of 7,286 judged reads flagged as chimeras.

**How exact is "the same"?** The random numbers of the test have a fixed seed. We ran these reads twice, with gene
databases downloaded on two occasions, and the gene calls and the result table were identical to the last digit.
Your numbers will differ slightly if CARD has a later release than 4.0.2 or if you use another Kraken2 index.

## 7. Reading this report

**The headline.** Of the flanks of resistance-gene units, 4.8% more have a mobility gene within 2 kb than chance
would give. The test finds this (p = 0.032): 132 linked units were seen where random positions in the same hosts
give 34. The interval, however, runs from −0.04 to 0.14 and includes 0. The report says why in the column "Notes":
*few effective hosts (G_eff 2.4)*. Most of the resistance genes of this sample sit in two or three species. The
interval is made by leaving out one host at a time, and with so few hosts each one moves the estimate a lot. The
honest reading is: resistance genes are linked to mobility genes in this sample, in a few hosts, and the size of the
link is not well known. More reads of the same sample would not narrow this interval much; more samples would.

**Metals.** MRG–MGE has almost the same estimate (0.037) but a p value just above 0.05 and a far wider interval.
The report says "no linkage shown" and adds that this is not evidence of no linkage. ARG–MRG is the opposite case:
an estimate near 0 with a narrow interval (−0.02 to 0.03), so strong linkage between the two classes can be ruled
out here.

**The two notes on QnrB.** The report flags two quinolone resistance genes:

> QnrB56: 245 carrier reads, 82% of them within 3% of 2.67 kb, the pattern of a small plasmid sequenced whole. The
> reads do not look like a cloning vector (30% of the read matches UniVec): probably a natural small plasmid; kept.

Many reads of one length that carry one gene are what a small plasmid gives when each read is the whole plasmid.
That can be a laboratory plasmid that got into the library, which would have to be removed. `mobilink` therefore
aligns such reads to the vector collection UniVec. Here one vector record matches about 30% of each read, below the
rule for a laboratory plasmid (35% and at least two records), so the reads stay in. The table "Genes on reads of one
fixed length" shows the numbers. When a cluster is judged to be a laboratory plasmid, the report adds a result table
without its reads, and that table is the one to use.

**Chimeras.** Without the 227 flagged reads the ARG–MGE estimate is 0.047 instead of 0.048. Chimeric reads do not
make this result.

**What this depth can show.** For ARG–MGE the smallest net linked fraction that the test finds with 80% power at
0.4 Gbp is 0.040. To find 0.02 would take 1.6 Gbp of the same library, and 0.01 would take 6.3 Gbp. For ARG–MRG
the limit is 0.085, so the test alone could not have shown weak linkage of these two classes at this depth.

**Which genes.** The tables at the end are descriptive. In this sample the calls next to mobility genes include
*sul1*, *dfrA7* and *qacEdelta1* with an integrase as the commonest neighbour (16 of 22, 18 of 24 and 12 of 21
calls), and efflux genes such as *mdtM* next to a transposase (25 of 30). The shares in these tables are not
corrected for chance or for read length, so they describe this sample and should not be compared between samples.

## 8. The same sample as contigs

`mobilink run --contigs` takes an assembly in place of reads. To try it, assemble the same reads with metaFlye
(`mamba install -c bioconda flye`; 6 minutes with 16 threads and 2 GB of memory here) and run the tool on the
contigs (1 minute with 8 threads):

```bash
flye --nano-hq SRR33541077.fastq.gz --meta --threads 16 --out-dir flye
mobilink run flye/assembly.fasta --contigs --name SRR33541077_flye \
    --db mobilink_db --kraken-db k2_standard_16gb --out results_contigs --threads 8
```

Flye 2.9.6 gave 599 contigs with 22.4 Mbp (N50 165 kb, 16 of them circular). Another version of Flye can give a
slightly different assembly, so expect numbers close to the ones below, not identical ones.

| Pair | Net linked within 2 kb | 95% interval | p | Units | Observed / expected linked units |
|---|---|---|---|---|---|
| ARG–MGE | 0.066 | −0.097 to 0.229 | 0.057 | 95 | 16 / 4.7 |
| MRG–MGE | 0.023 | −0.051 to 0.096 | 0.087 | 120 | 9 / 3.3 |
| ARG–MRG | 0.002 | −0.008 to 0.011 | 0.274 | 97 | 2 / 1.8 |

Calls: 130 ARG, 196 MRG and 377 MGE.

**Reads and contigs compared.** The estimates agree: 0.048 from the reads and 0.066 from the contigs for ARG–MGE,
and near zero for ARG–MRG in both. What differs is how much the data can show. The reads held 1,738 ARG units,
because every read that covers a locus brings its own copy; the assembly holds each locus once, 95 units in all.
The report therefore says "no linkage shown" for the contigs (p = 0.057), and its section on depth gives 0.15 as
the smallest net linked fraction that the test would find in this assembly, against 0.04 for the reads. With a
small sample the reads are the stronger data. Contigs become useful when the assembly is large, or when the
question needs longer distances than the reads span.

**What the contig mode does differently.** There are no overlap and chimera steps. All gene units of one contig
count as one clump in the test, which the notes of the report state. The contigs are searched in windows of 30 kb,
so that every copy of a repeated gene on a long contig is found.

## 9. Next steps

- All output files are listed in [One command, from reads to a report](one_command.md), the columns of the result
  table in [The output table](output_columns.md).
- For several samples, run each one and combine them with `mobilink meta`
  ([Step by step on real data](real_data.md), section "Groups of runs").
- Before sequencing more, `mobilink power` projects the depth needed from a pilot run such as this one.
