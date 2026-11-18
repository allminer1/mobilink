# Amendment 06: a second stage for H1a and H2

**Status: approved by the user on 2026-11-17 (the reply "фззкщму", which is "approve" typed on a Russian keyboard
layout; D72) and frozen the same day, before any stage-2 run was fetched or probed.** The text is draft v1 with this
status block; nothing else was changed. Only run metadata of the archive and the experiment records of the candidate
projects had been read. Frozen with: this file's sha256 (`docs/preregistration_amendment_06.sha256`), package commit
`9819182`, the frozen files of amendment 05, and the checksummed scripts and inputs in `frozen/amendment06/`.

## Why

- Amendment 05 (stage 1) was run once on 48 runs. H1a and H2 came out "not determined": 3 runs counted and 6 are
  needed. A run counts only if its subsample N50 is ≥ 9 kb, and the stage-1 sample was drawn to span each project's
  read lengths, not to be long.
- The user decided on 2026-11-17 to register an extension ("Register an extension").
- **What is known now and could bias a new choice:** the stage-1 results, including that all 3 counting runs meet both
  criteria and that soil runs have too few ARG to count. Stage 2 is therefore built so that this knowledge cannot
  steer it: the sample is fixed by archive metadata and a read-length rule, the decision rule and the code are those of
  stage 1, and the primary verdict uses the new runs only.

## Scope

- **Tested:** H1a and H2 (under linkage), with the Part D rule of amendment 05, word for word.
- **Not tested again:** H1b, H3 and H4. Their stage-1 results stand. Stage-2 runs are chosen for long reads, so they
  are not added to the class results of stage 1.
- The stage-1 verdicts are not replaced. The paper reports both stages.

## Part A. Sample

- **Source:** `public/screen/ena_ext_runs_2026-11-17.tsv` (ENA portal search of 2026-11-17: ONT, library source
  METAGENOMIC, ≥ 3 Gbp; 7,198 runs; sha256 `1ef9fe38…`). Script `scripts/P3X_01_candidates.py`.
- **Eligible run:** strategy WGS; selection RANDOM or unspecified; ≥ 3 Gbp; archive mean read length ≥ 6,000 bp;
  organism name of one of the three classes below; project not used in P0–P3 in any role; not a mock, culture,
  enrichment or laboratory reactor; native DNA (no PCR or whole-genome amplification).
- **Classes:** human gut, livestock, wastewater, by the archive's organism name
  (human gut: "human gut", "human feces", "gut", "feces" metagenome of human samples; livestock: "bovine gut",
  "bovine", "rumen", "manure", "chicken gut", "pig gut", "pig" metagenome; wastewater: "wastewater", "activated
  sludge", "sludge" metagenome).
- **Soil is left out.** In stage 1 no soil run had the 30 ARG sides that the counting rule needs (the registered ARG
  class is too rare in soil reads), so soil runs would use the budget without being able to count. This choice uses a
  stage-1 fact; it concerns ARG counts, not the outcome of H1a or H2.
- **Projects:** `public/screen/p3x_projects.tsv` lists every project the rule finds (36), with the sample type and
  the kit as the record gives them. 25 are included and 11 excluded (3 gut projects of other animals; 2 livestock
  cultures; 6 laboratory reactors or consortia). In PRJNA1398240 only the runs with alias "-raw-long" are eligible;
  the others are enrichments.

| Class | Project | Sample type (record) | Units | Longest mean read length (kb) | Kit and chemistry (record) |
|---|---|---|---|---|---|
| Human gut | PRJNA1139951 | human gut metagenome | 5 | 12.5 | SQK-LSK114, kit 14 |
| Human gut | PRJNA820119 | human stool | 164 | 11.5 | SQK-LSK109 with native barcodes; R9.4 |
| Human gut | PRJNA1236482 | human gut metagenome (patients) | 13 | 9.4 | Native Barcoding 96 kit; R10.4.1 |
| Human gut | PRJDB8879 | human gut metagenome (diet study) | 1 | 8.6 | not in the record |
| Human gut | PRJEB98069 | human gut metagenome (ileal pouch) | 3 | 8.2 | not in the record |
| Human gut | PRJNA763692 | human stool | 5 | 8.2 | not in the record |
| Human gut | PRJEB88320 | stool of six healthy volunteers | 2 | 8.1 | not in the record |
| Human gut | PRJNA1000750 | human gut metagenome | 2 | 7.7 | not in the record |
| Human gut | PRJEB49168 | human stool | 4 | 6.5 | not in the record |
| Livestock | PRJNA1398240 | swine slurry ("-raw-long" runs) | 3 | 11.5 | not in the record |
| Livestock | PRJEB105185 | rumen of dairy cows | 9 | 9.1 | not in the record |
| Livestock | PRJNA1011201 | chicken faeces | 1 | 7.3 | SQK-RAD004 rapid kit |
| Livestock | PRJNA1162230 | rumen of dairy cows | 1 | 7.1 | ligation kit LSK112; R9.4.1 as stated |
| Livestock | PRJNA1466632 | rumen liquid | 2 | 6.9 | not in the record |
| Livestock | PRJEB83989 | rumen of beef cattle | 6 | 6.8 | not in the record; barcoded |
| Wastewater | PRJNA1134722 | sewage influent, one plant | 1 | 12.0 | SQK-LSK114; FLO-PRO114 |
| Wastewater | PRJEB83983 | activated sludge, plants worldwide | 75 | 11.9 | not in the record |
| Wastewater | PRJNA1067729 | activated sludge | 5 | 9.9 | not in the record |
| Wastewater | PRJEB67571 | treatment plant | 1 | 8.0 | R9.4.1 (file name) |
| Wastewater | PRJNA629478 | activated sludge of treatment plants | 20 | 7.9 | SQK-LSK109; R9 |
| Wastewater | PRJNA955434 | activated sludge and wastewater | 2 | 7.5 | not in the record |
| Wastewater | PRJDB12729 | activated sludge, landfill leachate | 1 | 7.3 | not in the record |
| Wastewater | PRJNA1028177 | wastewater | 1 | 7.2 | not in the record |
| Wastewater | PRJNA1437919 | activated sludge, municipal plant | 1 | 7.1 | not in the record |
| Wastewater | PRJNA1005436 | activated sludge | 1 | 7.1 | SQK-LSK112 |

- **Unit:** the BioSample; the donor in PRJEB88320 and the treatment plant in PRJNA1134722 (their samples are time
  points). The run with the longest archive mean read length represents its unit. 329 units: gut 199, livestock 22,
  wastewater 108.
- **Order, with no random element** (`public/screen/p3x.order.tsv`): within a project, units by archive mean read
  length, longest first; within a class, projects by the mean read length of their first unit, longest first. The try
  order is the first unit of every project in project order, then the second units, and so on.
- **Walk, per class, down the try order.** A candidate is fetched and annotated with the frozen stage-1 run script
  (`frozen/amendment05/P3_02_public_run.sh`: 3 Gbp subsample by seeded blocks, the registered annotation). It is
  **accepted** if:
  - its project has fewer than 2 accepted runs;
  - the subsample has ≥ 1 Gbp and N50 ≥ 9,000 bp;
  - the kit rule of amendment 05 Part A holds as applied in stage 1 (kit or chemistry known from the record, the
    linked publication or the read ends; a library found to be PCR- or whole-genome-amplified is not eligible).
- A candidate that is not accepted is listed with the reason, and no statistic is computed on it. A fetch that fails
  three times is listed and skipped.
- **Stop** at 6 accepted runs in the class, after 12 candidates tried in the class, or when the list ends. At most
  18 runs are accepted and at most 36 fetched.
- R9 and R10 runs both count, as in stage 1; the result is also given by chemistry.

## Part B. Per run

- As amendment 05 Part B, unchanged: registered classes, Kraken2 strata, unit loci from the carrier overlaps, the
  registered net linked test, and the frozen Part B driver for the naive share and the cut to N50 3 kb (20 draws).
- Code: `frozen/amendment05/P3_03_confirm.py` and package commit `9819182`, both unchanged. The package stays frozen
  until stage 2 has run.

## Part C. Decision rules

- **H1a (stage 2).** The rule of amendment 05 Part D on the accepted stage-2 runs alone. A run counts if its subsample
  N50 is ≥ 9 kb, at least 30 sides reach the window in both versions and the 3 kb half-width is ≤ 0.10.
  - Supported if at least two-thirds of the counting runs have naive share native / 3 kb ≥ 1.5.
  - Not supported if fewer than two-thirds do.
  - Fewer than 6 counting runs: not determined.
- **H2 under linkage (stage 2).** The same counting runs. Supported if at least two-thirds are stable
  (|native − 3 kb| ≤ max(0.02, 0.25 × |native|)) with overlapping intervals; otherwise as above.
- **Why stage 2 alone is the primary result:** the 3 counting runs of stage 1 are already known to meet both
  criteria. Adding them to the count would make the test easier to pass. They are shown next to it.
- **Reported next to the verdict, without a rule:** both stages together (3 + the stage-2 counting runs); the result
  by class and by chemistry; every accepted run that does not count, with the reason; every candidate that was not
  accepted.
- **Wording fixed now.** "Stage 1 (amendment 05): not determined, 3 of 3 counting runs met the criterion, 6 were
  needed. Stage 2 (amendment 06, runs chosen for N50 ≥ 9 kb after stage 1 was known): [verdict], k of n counting
  runs." A stage-2 verdict of "not supported" or "not determined" is reported in the same place and the same words.

## Part D. Order of work and freeze

1. The user approves this text.
2. Frozen before any fetch: this file's sha256; `p3x_projects.tsv`, `p3x.eligible.tsv`, `p3x.order.tsv`, the
   archive table and the PRJNA1398240 alias table; the selection, walk and analysis scripts; the package commit.
3. The walk runs. Read lengths, the kit check and the acceptance list are the only things looked at.
4. When the walk has ended in all three classes, the frozen analysis runs once on the accepted runs
   (`tables/P3/confirm_ext/`).
5. Deviations are logged in `docs/preregistration_amendment_06_deviations.md`.

## Feasibility (archive metadata and stage 1 only)

- In stage 1 the subsample N50 was 1.2 to 4.0 times the archive mean read length (median 1.56). A mean of 6 kb then
  gives N50 ≥ 9 kb in about half of the runs, and a mean of 7.5 kb in most. Units with a mean ≥ 7.5 kb: gut 116
  (8 projects), wastewater 54 (5 projects), livestock 10 (2 projects).
- Livestock is the weak class: 22 units in 6 projects, most of them rumen, where ARG may be too few for a run to
  count. The verdict needs 6 counting runs over all classes, not per class.
- Cost: up to 36 runs at about 15 min each with three in parallel, then about 1.5 h of analysis.
