# Phase 2 report: uncertainty layers

2026-11-16. Phase 2 is complete. This report summarises the steps and lists the decisions needed before Phase 3.
Details are in the step reports in `reports/P2/` and in decisions D35–D52 of `docs/analysis_plan.md`.

## In short

- **The read-level test changed, and it is calibrated.** Phase 1's molecule-level O/E test turned small errors in the
  chance model into false positives once counts were large. Phase 2 replaced it with the **net linked test**:
  - it compares each resistance-gene unit with random anchors on the same molecules, within host strata;
  - reads that cover the same locus count as one clump.

  Amendment 04 registered it (approved and frozen 2026-11-16). Gate G3 passed: size 5.0% and 4.4% in 16,200 simulated
  datasets, interval coverage ≥ 0.91 in every setting, and false-discovery share 0.031 and 0.029.
- **Real-data artefacts are small.**
  - Chimeras make few false links: 0.67% of same-molecule links in 7 public mock communities, none in the gut mock
    with our kit.
  - Kraken2 host strata are 99.6–99.9% pure.
  - The result does not hinge on one annotation threshold: no sign changes under 22 annotation variants.
- **The tool exists:** `linkstat` v0.2.0.dev0 (local commits only), with `linkstat test` running the registered test
  and `linkstat chimeras` for chimera flags.
- **First look at public runs (development data):** ARG:MGE linkage is positive in all 8 G1b runs. ARG:MRG is about 0.
  Where one host dominates a sample, the intervals are honestly wide.

## Steps

| Step | What | Result | Report |
|---|---|---|---|
| 2.0 | Re-fetch of the 8 G1b runs: carrier overlaps, AMRFinderPlus, Kraken2 on all reads | done for all 8 | interim report |
| 2.1 | Units (clustered classes) | N2-ev-u removes most of the G2 bias; kept as descriptive O/E | interim report |
| 2.2 | Simulator extensions (islands, host heterogeneity, plasmids, chimeras) | done | interim report |
| 2.3 | Locus clusters from overlaps | purity 0.87 at the density rule; used as unit loci | interim report |
| 2.4 | Locus-aware inference | molecule-level O/E fragile at high counts; led to the net linked test (D39) | interim report |
| 2.5 | Chimeras and mock communities (V4) | 7 runs, 5 labs; mock gate passes in all 7 (D48 reading) | `P2_05_mocks_report.md` |
| 2.6 | Annotation sensitivity (V6) | no sign reversal in 22 variants; registered thresholds on the plateau | `P2_06_annot_sens_report.md` |
| 2.7 | Amendment 04 and gate G3 | approved and frozen; G3a and G3b pass | `G3_report.md` |
| 2.8 | v0.2 | `linkstat test`, `linkstat chimeras`, reporting fixes; the shared test seed does not matter | D52 |

## What Phase 2 showed

1. **Host co-variation and locus pseudo-replication are the two big threats.** A-rich genomes are B-rich (association
   without linkage), and reads that cover one locus repeat one observation. The net linked test handles both: anchors
   within host strata, and loci as clumps. Without clumps the test rejected 41% under no linkage.
2. **Chimeras vary more than 100-fold between runs (0.01–5.5% of reads), but rarely create links.** Even the run with
   5.5% chimeric reads had only 1.0% false links. The chimera filter as first specified cost far more true links than
   it removed false ones; v0.2 uses a stricter rule.
3. **The size of the linkage depends on the class definitions.** ARG:MGE rises with the ARG identity threshold from
   40% to 70% and is flat from 70% to 90% (registered: 80%). The broader D3 MGE class gives larger values than the
   mobility core. The paper should give the class definitions with every estimate.
4. **One dominant host means little information.** In SRR32022521, 82% of resistance-gene units sit in E. coli, on 223
   loci covered by up to 724 reads, so there are about 1–2 effective hosts. The test reports this, and the intervals
   are wide.

## Open items carried into Phase 3

- The development runs looked conservative (net_w 3.3–3.7%) while G3a is at nominal size (5.0%). No cause was found:
  - not the code (exact reproduction);
  - not the seeds;
  - not dependence between datasets;
  - not the shared test seed (D52).

  The paper does not claim the test is conservative.
- Fold-back reads (0.07–1.7%) count a locus twice on one read; detection needs the reads and is left to the P3 data
  passes.

## Decisions for you

1. **Go to Phase 3?** Phase 3 analyses our 8 deep soils, the low-coverage samples and 30–50 public samples, and builds
   the package for release (docs, Bioconda, benchmark).
2. **The mock gate rule (D48).** The frozen text gates the mocks "only if they carry ≥ 10 true A–B pairs" but gives no
   threshold. My reading: in each mock run, the false share of read-level links after chimera filtering is ≤ 5% at its
   upper 95% bound. All 7 runs pass, also before filtering. Accept, or set another rule?
3. **Chimera filter in Phase 3.** The registered test has no chimera filter. Proposal: the primary analysis without a
   filter (as registered) and the filtered one as a sensitivity analysis.
4. **Read-to-bin rule for our soils** (open since D46; needed before the soils are tested). The amendment's default is
   the bin of the contig the read maps to (minimap2, MAPQ ≥ 20), else "unresolved". Alternative: Kraken2 for the
   soils too, as for public runs. Kraken2 would leave many soil reads unresolved, because soil taxa are poorly
   represented in the database.
5. **Before the soils (open since the handoff):**
   - permission to reuse the Pervomaisk samples and the 36 industrial samples in this paper;
   - whether the 36 include the 8 deep samples, so that no sample is counted twice.

## Files

- Reports: `reports/P2/G3_report.md`, `reports/P2/P2_05_mocks_report.md`, `reports/P2/P2_06_annot_sens_report.md`,
  `reports/P2/interim_2026-11-15.md`.
- Tables: `tables/G3/`, `tables/P2/P2_05_mocks.tsv`, `tables/P2/P2_05_soil_adapters.tsv`,
  `tables/P2/P2_06_annot_sens*.tsv`, `tables/P2/P2_08_seed_check.tsv`.
- Package: `pkg/` at `cf0db83` (v0.2.0.dev0, 100 tests).
