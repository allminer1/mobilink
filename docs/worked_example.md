# Worked example

This example takes about one minute on a laptop. It uses simulated reads, so the true answer is known and you can
check that mobilink finds it. Every number below was produced with the commands shown (mobilink 0.5, default seeds).

## 1. Simulate a community

```bash
mobilink simulate --read-lengths examples/read_lengths.tsv.gz --n-genomes 60 --total-bases 1e9 --linked 0.3 \
  --out-prefix toy
```

```
80,435 reads, 26,244 features; true net linked (far-end distance, w = 2 kb): 0.1251
```

- `examples/read_lengths.tsv.gz` holds 20,000 read lengths drawn from a public human gut run (SRR29980985, read N50
  18.5 kb). The simulator draws read lengths from it.
- The community has 60 genomes. ARG are class A and MGE are class B. 30% of the ARG are placed next to an MGE; the
  rest are placed at random.
- Three files are written:
  - `toy.lengths.tsv.gz`: one row per read (`molecule_id`, `length`, `genome`);
  - `toy.features.tsv.gz`: one row per gene call on a read (`molecule_id`, `start`, `end`, `class`, `scov`);
  - `toy.truth.tsv`: the true net linked fraction, computed on the genomes.

The truth is not 0.30. Some of the ARG that were placed at random also land near an MGE, and the net linked fraction
removes exactly that chance share. Also, a pair counts only if the whole MGE lies within 2 kb. The true value for
this community is 0.125.

## 2. Run the test

```bash
mobilink test --lengths toy.lengths.tsv.gz --features toy.features.tsv.gz --pairs ARG:MGE --strata-col genome \
  --out toy.results.tsv
```

`--strata-col genome` names the host of each read. With real reads, the host comes from Kraken2 (`--kraken`).

| Column | Value | Meaning |
|---|---|---|
| `A_calls`, `B_calls` | 2,775; 23,457 | ARG and MGE calls |
| `net_w` | 0.132 | Net linked fraction within 2 kb (truth 0.125) |
| `net_w_lo`, `net_w_hi` | 0.103; 0.161 | 95% interval |
| `net_w_O`, `net_w_E` | 694; 59.2 | ARG sides with an MGE within 2 kb: observed, and expected by chance |
| `net_w_p` | 0.0001 | Test against chance (the smallest value that 9,999 draws can give) |
| `net_mol` | 0.316 | Net linked fraction for "an MGE anywhere on the same read" |
| `mol_ratio` | 2.15 | ARG carry an MGE on their read 2.15 times as often as chance predicts |
| `hosts_effective` | 10.9 | Effective number of hosts behind the estimate |
| `flags` | no locus clumps (each unit occurrence is its own clump) | See below |

The flag is expected here. The example gives no read overlaps, so each read counts as its own observation. With real
reads at more than about 1× coverage, give the overlaps with `--loci`; without them, p values are too small.

## 3. Cut the reads

```bash
mobilink oe --lengths toy.lengths.tsv.gz --features toy.features.tsv.gz --pairs ARG:MGE --strata-col genome \
  --fragment 9000 --fragment 4000 --out toy.cut.tsv
```

`--fragment N` cuts the reads at random points until their N50 is N, and estimates again.

| Read N50 | Naive share (`*_naive`) | Net linked fraction (`*_net`) |
|---|---|---|
| 18.5 kb (`intact_`) | 0.422 | 0.131 |
| 9 kb (`r9000_`) | 0.307 | 0.130 |
| 4 kb (`r4000_`) | 0.204 | 0.135 |

The naive share (ARG reads that also carry an MGE) halves. The net linked fraction does not move. This is why two
studies with different read lengths can compare net linked fractions but not naive shares.

## 4. How much data would you need?

```bash
mobilink power --lengths toy.lengths.tsv.gz --features toy.features.tsv.gz --a ARG --b MGE \
  --bases 5e8 --bases 1e9 --bases 5e9
```

```
 pilot     bases  molecules  expected  oe_enrichment  oe_depletion   sides_w  halfwidth
 False     5e+08  4.024e+04     288.6          1.147        0.8595      2018    0.01617
  True 9.994e+08  8.044e+04     576.8          1.103        0.9003      4033    0.01144
 False     1e+09  8.048e+04     577.2          1.103        0.9004      4036    0.01144
 False     5e+09  4.024e+05      2886          1.045        0.9553 2.018e+04   0.005114
bases for a net linked half-width of 0.05 at w = 2000: 5.23e+07
```

- `expected`: read-level co-occurrences expected by chance at that depth.
- `oe_enrichment`, `oe_depletion`: the observed-to-expected ratios that the read-level test can detect.
- `halfwidth`: expected half-width of the interval of the net linked fraction.

The gene densities and read lengths are taken from the pilot data, so the answer holds for a sample like the pilot.
This simulated community is rich in ARG. A soil sample with a read N50 of 5 kb can have fewer than one expected
event in 10 Gbp, and then no depth that is practical gives an estimate.

## Next

`docs/real_data.md` shows the same analysis on real reads and on contigs.
