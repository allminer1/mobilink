# Amendment 05: deviations and rule applications after the freeze

The frozen text (`preregistration_amendment_05.md`) is not edited. Every departure from it, and every application of
its hand-checked rules, is listed here with its date and with what was known at the time.

## Deviation 1 (2026-11-17, user decision D66): the CARMA gut runs are replaced from the other gut project

- **What the frozen sample has:** 6 gut runs from PRJEB108131 (CARMA trial) and 6 from PRJNA929328.
- **What was found:** the paper that now cites PRJEB108131 (PMID 42878177, doi 10.1080/19490976.2026.2746814) states
  the Rapid PCR Barcoding Kit SQK-RPB004 (14 PCR cycles) on R9.4.1 flow cells. The libraries are PCR-amplified. The
  pre-registration allows only shotgun native DNA, and amendment 05 dropped PRJNA1088077 for the same reason. The
  archive metadata (selection "RANDOM") did not show it; the paper was not indexed under the accession on 2026-11-16.
- **Why this is a deviation:** the frozen replacement rule lists four reasons and draws the replacement from the same
  project. Every run of PRJEB108131 has the same library type, so the rule cannot be applied as written.
- **Decision (user, 2026-11-17, "okay i agree" to the recommended option):** the 6 runs are excluded. They are
  replaced by the next 6 unused patients on the frozen replacement list of PRJNA929328, in its seeded order:
  SRR23278684, SRR23278671, SRR23278685, SRR23278715, SRR23278728, SRR23278753. The gut class keeps 12 runs, all from
  one project and all R9.4. The frozen replacement rule applies to these runs as to any other.
- **What was known when it was decided:** no linkage result of any confirmatory run had been computed or looked at.
  Known: read counts, subsample sizes, N50 values and the read-end adapter shares.
- **Reported in addition:** the 6 CARMA runs are analysed separately as a non-registered sensitivity analysis
  (per-run estimates, and the pooled gut result with them added). It changes no decision.
- **Files:** `public/screen/p3_manual_exclusions.tsv`, `public/screen/p3_runlist.tsv` (role "replaces ... deviation 1").

## Rule application 1 (2026-11-17, D63 and D65): kit and chemistry

- **Reading of "the archive record":** it includes the experiment's design description and the read names of the
  submitted files. Fixed before any linkage result existed.
- **Read-end rule (fixed before the tables were read):** a pattern counts if it lies at the start of at least 5% of
  the reads (`scripts/P3_02c_kits.py`).
- **Result by project:** `public/screen/p3_project_kits.tsv`, `p3_kit_probe.tsv`, `p3_record_design.tsv`.
- **Excluded by the frozen rule ("kit and chemistry both undeterminable"):** ERR17169882 and ERR17169874 of ERP183232
  (reads trimmed; nothing in the record; no publication). Replaced from the project's frozen list.
- **Kept:** the four trimmed runs of PRJNA1414918, because the record states their kit ("Rapid Barcoding kit").

## Rule application 2 (2026-11-17, D67): the replacement runs and the final list

- **Kit rule on the 8 replacements** (`public/screen/p3_kit_probe.tsv`): all pass. Chicken ERR17169904 shows the
  ligation adapter at the start of 81% of reads and ERR17169864 a native barcode at 82%. The 6 gut replacements are
  of PRJNA929328, whose kit and chemistry are published.
- **Soil projects:** all 12 runs pass (kit or chemistry stated for each project). Within PRJEB86780 one run shows
  native barcodes and one the ligation adapter at the read starts; the chemistry is the same (R10.4.1), so both stay.
- **Final list:** `public/screen/p3_final.tsv`, 48 runs, 12 per class; 56 runs fetched, 8 excluded and replaced, no
  failed fetch and no run below 1 Gbp or N50 2 kb. The driver ended after round 2 with 0 further replacements.
- **Frozen analysis:** run once on this list (2026-11-17 08:47–10:18 UTC); the freeze check passed before the run.

## Additions outside Part E (2026-11-17, D69)

- `scripts/P3_04_sensitivity.py` also runs the six PCR-amplified gut runs of deviation 1 (pooled with the registered
  gut runs, pooled alone, each alone), as agreed in D66. No decision uses them.
- `scripts/P3_03b_descriptive.py` gives the random-effects summary named in Part C and the Part B columns for all
  runs, not only the counting ones. The registered verdicts are those of `tables/P3/confirm/verdict.txt`.

## Note for later draws

- The frozen replacement list of PRJNA929328 holds SRR23278737 (sample alias "Zymo", order 15), a mock control that
  the title-based screen did not catch. If the list is ever drawn that far, this run is skipped as ineligible
  ("not a mock, reference or synthetic sample").
