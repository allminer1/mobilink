# From reads or contigs to a result

`mobilink run` does all of this in one command (see `one_command.md`). This page gives the single steps, for
pipelines that need them apart. These are the commands behind the public-run analyses of the paper, with the tool versions used there: DIAMOND 2.2.8,
Kraken2 with the Standard-16 database (2024-12-28), minimap2 2.27 and seqkit. Other versions should work, but the
thresholds were validated with these.

`T` is the number of threads. `reads.fa` is a FASTA or FASTQ file of nanopore reads (gzip is slower for DIAMOND; an
uncompressed FASTA was 2.6 times faster in our runs).

## 1. Databases (once)

```bash
diamond makedb --in protein_fasta_protein_homolog_model.fasta -d card_homolog        # CARD, protein homolog model
diamond makedb --in mobileOG-db_beatrix-1.6.All.faa -d mobileog_all                  # mobileOG-db beatrix-1.6
diamond makedb --in BacMet2_EXP_database.fasta -d bacmet2_exp                        # BacMet2, experimentally confirmed
python resources/make_mobileog_lookup.py mobileOG-db_beatrix-1.6.All.faa mobileog_all.lookup.tsv.gz
```

`mobilink db --out mobilink_db` does the same and also downloads the files. Two small tables that `mobilink prep`
needs come with the package: the BacMet2 class table (which entries are metal resistance genes) and the mobility core
(which mobileOG-db entries count as MGE). `mobilink db` copies them into the database directory; in the source tree
they are in `src/mobilink/data/`. See `resources/README.md`.

## 2. Read lengths

```bash
{ printf 'read_id\tlength\n'; seqkit fx2tab -n -i -l reads.fa; } | gzip > reads.lengths.tsv.gz
```

## 3. Gene calls

```bash
for db in card_homolog mobileog_all bacmet2_exp; do
  diamond blastx -q reads.fa -d $db -F 15 --range-culling --top 10 -p $T \
    -f 6 qseqid sseqid pident length qstart qend sstart send slen qlen evalue bitscore -o reads.$db.tsv
done

mobilink prep --lengths reads.lengths.tsv.gz \
  --card reads.card_homolog.tsv \
  --bacmet reads.bacmet2_exp.tsv --bacmet-classes mobilink_db/bacmet2_exp_classes.tsv \
  --mobileog reads.mobileog_all.tsv --mobileog-lookup mobileog_all.lookup.tsv.gz \
  --mobility-core mobilink_db/mobility_core.tsv.gz \
  --rule primary --out-prefix reads
```

- `-F 15 --range-culling` is DIAMOND's long-read mode. It tolerates frameshifts and reports several genes per read.
- The column list after `-f 6` is required in this order.
- `mobilink prep` writes `reads.primary.features.tsv.gz` and `reads.summary.tsv` (calls per class).
- The `primary` rule is the one used in the paper:

| Class | Database | Identity | Subject cover |
|---|---|---|---|
| ARG | CARD protein homolog model | ≥ 80% | ≥ 80% |
| MRG | BacMet2 experimentally confirmed, metal entries | ≥ 80% | ≥ 80% |
| MGE | mobileOG-db, mobility core | ≥ 60% | ≥ 80% |

Hits are collapsed to one call per ORF and database. An ORF that is called both as class A and as class B is removed
from that pair by `mobilink test`, so a single gene cannot link to itself. mobileOG-db calls outside the mobility core
get the class `MGE_other`; leave them out of the pairs unless you want them.

## 4. Hosts (reads)

```bash
kraken2 --db k2_standard_16gb --threads $T --confidence 0 --output reads.kraken.out reads.fa > /dev/null
awk -F'\t' -v OFS='\t' 'BEGIN{print "read_id","taxid","hit_kmers"}
  {h=0; n=split($5,p," "); for(i=1;i<=n;i++){split(p[i],q,":"); if(q[1]!="0" && q[1]!="A") h+=q[2]}
   print $2,$3,h}' reads.kraken.out | gzip > reads.kraken.tsv.gz
```

Classify all reads, not only the gene-carrying ones: the chance model places random anchors on the reads of each
host. Reads with fewer than 10 classified k-mers go to an "unresolved" stratum (`--min-hits`).

If the library contains the nanopore lambda control strand (DCS), remove those reads first. Kraken2 with the standard
databases labels them *Escherichia coli*, so align to the lambda genome instead:

```bash
minimap2 -x map-ont -c NC_001416.1.fa reads.fa | awk '$11 >= 200 {print $1}' | sort -u > lambda_reads.txt
```

and pass `--exclude lambda_reads.txt` to `mobilink test`.

## 5. Overlaps among the gene-carrying reads

```bash
zcat reads.primary.features.tsv.gz | awk -F'\t' 'NR>1{print $1}' | sort -u > carriers.txt
seqkit grep -f carriers.txt reads.fa > carriers.fa
minimap2 -x ava-ont -t $T carriers.fa carriers.fa | awk -F'\t' '$11 >= 1000' | cut -f1-12 | gzip > carriers.ava.paf.gz
```

Reads that cover the same genomic locus are not independent observations. mobilink uses these overlaps to join them
into clumps. Only reads with a class-A call are needed, so you can restrict `carriers.txt` to those if the file is
large.

## 6. Optional: chimeric reads

```bash
mobilink chimeras --paf carriers.ava.paf.gz --out chimeras.tsv
```

A chimera joins two unrelated molecules and can create a false link. On mock communities the default settings flagged
73–100% of reads that join two species and 0.6–2.2% of clean reads.

In real communities the rule flags more reads (a median of 12% of the judged reads, 4–39%, in 64 public runs),
because strains that differ by an insertion or a rearrangement also produce overlap breakpoints.
Some of these are true links. Treat the filtered result as a sensitivity analysis: run `mobilink test` with and
without `--exclude chimeras.tsv` and report both. The registered analyses of the paper use no filter.

## 7. Test

```bash
mobilink test --lengths reads.lengths.tsv.gz --features reads.primary.features.tsv.gz \
  --pairs ARG:MGE --pairs ARG:MRG --pairs MRG:MGE \
  --kraken reads.kraken.tsv.gz --taxonomy k2_standard_16gb/ktaxonomy.tsv \
  --loci carriers.ava.paf.gz --out results.tsv
```

`results.tsv` has one row per pair. `results.tsv.meta.json` records the version, the settings and the checksums of
the inputs. `docs/output_columns.md` explains the columns.

### Time and memory

Measured on 64 public runs of about 3 Gbp each (0.24–1.9 million reads), one thread:

| Command | CPU time, median (range) | Peak memory, median (range) | Grows with |
|---|---|---|---|
| `prep` | 0.2 min (0.05–6.3) | 1.2 GB (0.5–26) | DIAMOND hit rows: about 0.3 GB per million rows |
| `chimeras` | 0.08 min (0.02–7.2) | 0.5 GB (0.1–1.1) | overlap rows |
| `test` | 0.8 min (0.2–7.7) | 1.2 GB (0.3–2.1) | reads and gene calls |

The server was shared, so CPU times are approximate: `prep`, whose code did not change, took about twice as long in
an earlier measurement under a higher load. Version 1.0 removed the two large memory costs of version 0.5 (`chimeras` up to 104 GB, `test` up to 14 GB). The
results of `test` did not change; the flags of `chimeras` changed only by the fix named in the changelog (0.3% of
the flags):

- `mobilink test`: the descriptive O/E columns no longer build a matrix over all host strata. A soil run with 0.98
  million reads and 9,883 strata now needs 0.7 GB and 82 s of CPU time (before: 12.9 GB and 152 s). `--no-oe` still
  skips these columns; the estimates, intervals and p values are the same with and without them.
- `mobilink chimeras` reads the overlap file in one pass and keeps only small records on disk while it runs (44
  bytes per usable overlap, in the system's temporary directory or in `--tmp-dir`). The run with 445 million overlap
  rows now needs 0.65 GB and 7.5 CPU minutes (before: 104 GB and 21 CPU minutes). If the temporary directory of your
  system is held in memory (a `tmpfs`), give a directory on disk with `--tmp-dir`.

The direction of a pair matters. `ARG:MGE` asks what share of ARG are linked to an MGE; `MGE:ARG` asks what share of
MGE are linked to an ARG.

## Contigs

Annotate the contigs as in step 3 (for contigs longer than about 100 kb, cut them into overlapping windows for
DIAMOND, or use protein-level annotation and build the feature table yourself). Then:

```bash
mobilink test --lengths contigs.lengths.tsv --features contigs.features.tsv --pairs ARG:MGE \
  --strata-col bin --molecule-clumps --out contig_results.tsv
```

- The lengths table needs a column with the bin (genome) of each contig; give unbinned contigs one common label.
- `--molecule-clumps` treats all A genes of a contig as one observation. Without it, the same-contig test rejected
  up to 19% of simulated datasets with no linkage; with it, 0–3%.
- On contigs, read `net_mol`, `net_mol_p` and `mol_ratio` (same contig) first. `net_w` works as for reads.

## Several samples

Run each sample on its own. To pool samples, concatenate the tables and the overlap files (read names must be
unique), make the host label unique per sample (for example `sample|species`), and pass it with `--strata-col`.
Anchors then never leave their sample and host.
In Python, `mobilink.pairs.pool_runs` does this and keeps the overlaps of each sample apart.

## Groups of runs: environments, treatments, studies

To compare groups of runs, use `mobilink meta`. It needs a manifest with one row per run:

| Column | Content |
|---|---|
| `run` | A unique name |
| group column (default `group`) | What is compared, for example the environment |
| cluster column (default `cluster`) | The unit whose runs are not independent, for example the project or the site |
| `lengths`, `features` | The two tables of the run (paths; relative paths start at the manifest) |
| `kraken` (optional) | Kraken2 labels per read, as for `mobilink test`; needs `--taxonomy` |

```bash
mobilink meta --manifest runs.tsv --pair ARG:MGE --group-col environment --cluster-col project \
  --taxonomy k2_standard_16gb/ktaxonomy.tsv --out survey
```

- `survey.groups.tsv`: one row per group. `net_w` is the net linked fraction within 2 kb of the pooled runs of the
  group. The interval (`lo95`, `hi95`) comes from leaving out each cluster in turn. `clusters_effective` shows how
  evenly the class A units are spread over the clusters: a value near 1 means that one cluster carries the estimate.
  `net_w_min_without_one_cluster` and `net_w_max_without_one_cluster` show how far one cluster can move it.
- `survey.jackknife.tsv`: the estimate without each cluster.
- `survey.contrasts.tsv`: the difference between every two groups, with a Welch test on the jackknife standard
  errors and Holm-adjusted p values.

A group needs at least two clusters for an interval, and the interval is only as good as the number of clusters: with
fewer than about six it is very wide. Runs of one project share sampling site, extraction, library and basecalling,
so the project, not the run, is the unit of replication. The command gives estimates and contrasts only; for the
test of one run or one pooled group against chance, use `mobilink test`.

## Your own classes

Any feature table works. For example, to test whether integrons sit near toxin genes, write a table with
`molecule_id`, `start`, `end`, `class` holding `INT` and `TOX` rows and run `--pairs TOX:INT`.
