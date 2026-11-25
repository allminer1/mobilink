# mobilink

mobilink tests whether two gene classes sit close together on long reads or contigs more often than chance, and
measures how strongly they are linked with a value that does not depend on read length. It was built for
antibiotic resistance genes (ARG), metal resistance genes (MRG) and mobile genetic elements (MGE) in long-read
metagenomes, and it accepts any two classes of features.

## Where to start

| You have | Read |
|---|---|
| A FASTQ file or an assembly, and you want a result | [One command, from reads to a report](one_command.md) |
| A first run on public reads, with the expected numbers | [Tutorial on a public sample](tutorial.md) |
| Gene calls of your own, or other feature classes | [Step by step on real data](real_data.md) |
| No data yet | [Worked example with simulated reads](worked_example.md) |
| A result table to interpret | [The output table](output_columns.md) |
| Questions about the estimator and the test | [Method](method.md) |

## Install

```bash
mamba install -c conda-forge -c bioconda mobilink diamond minimap2 kraken2
```

or, from the source, `pip install .` in a clone of the repository. A `Dockerfile` in the repository builds an image
with mobilink, DIAMOND, minimap2 and Kraken2:

```bash
docker build -t mobilink .
docker run --rm -v "$PWD":/data mobilink mobilink db --out /data/mobilink_db
docker run --rm -v "$PWD":/data mobilink mobilink run /data/reads.fastq.gz --db /data/mobilink_db --out /data/results
```

The databases and the Kraken2 index are not part of the image.

## The measure in one paragraph

Each gene unit of the first class has two flanks. The net linked fraction is the share of flanks with a whole gene
of the second class within a window (2 kb by default), beyond what the gene densities of the same hosts give by
chance. It is 0 under chance and stays the same when reads get shorter. The usual share of reads that carry both
classes grows with read length and has no chance expectation; do not compare it between studies.
