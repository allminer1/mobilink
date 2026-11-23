# Amendment 07: deviations and notes

Kept from the start of the data pass (2026-11-18). A deviation is a departure from the frozen text; a note records how
a rule of the frozen text was applied where the text leaves a detail open. Entries are dated and are written before the
result they could affect is known.

## Deviations

None so far.

## Notes

1. **2026-11-18 09:50 UTC. Amplified libraries known from the records.** The NCBI experiment records
   (`public/screen/p4_record_kits.tsv`) state a PCR barcoding kit for PRJNA1231276 (marine; SQK-RPB114.24) and
   PRJNA942915 (livestock; RPB004). By Part A ("amplified libraries are not eligible") all their runs are rejected.
   They are listed in `public/screen/p4_manual_exclusions.tsv`. The frozen walk applies this list to fetched runs
   only, so the runs were fetched before they were rejected. No statistic of these runs is computed.
2. **2026-11-18 09:50 UTC. Runs that the archive does not serve yet: one retry.** The frozen run script fetches from
   NCBI. Eight selected runs were released at ENA on 2026-11-05/06 or in October 2026 and NCBI holds no record of
   them yet (esearch count 0 on 2026-11-18): the 5 runs of PRJEB104244 (soil), the 2 of PRJEB125553 (soil) and the 1
   of PRJEB123233 (human gut). Their fetch fails within minutes. Part A says "a run is rejected if the fetch fails"
   and does not say how often a fetch is tried. Rule, set now, before any run of the sample has been analysed:
   - a selected run whose fetch fails while NCBI holds no record of its accession is tried once more, when the first
     pass over the sample has ended (its `P3_FAILED` marker is removed and the driver is started again);
   - if it fails again, it is rejected for good;
   - replacements of such a run that the walk tried in the meantime, and that failed for the same reason, are not
     tried again;
   - a fetch that fails although NCBI holds the record is rejected at once, as the text says.
   The first failure of this kind: ERR16963751, 09:46 UTC ("empty block").
3. **2026-11-19 04:15 UTC. A fetch that failed although NCBI holds the record.** ERR15289757 (PRJEB88618, soil;
   "empty block") is rejected for good, as note 2 says; the project has no other soil unit.
4. **2026-11-19 04:15 UTC. Runs with megabase-long reads are kept; one exploratory analysis is added.** In three soil
   runs of PRJEB74630, reads longer than 200 kb hold 47%, 40% and 10% of the bases (longest read 3.3 Mb; mean read
   length 2.2–3.8 kb), so their N50 (179 kb and 130 kb in two of them) does not describe the library. Such reads are
   most likely not single DNA molecules. The frozen text has no rule for this, so the runs stay in the sample and in
   every registered analysis. Seen in a scan of read lengths only (93 accepted runs so far; 4 runs have more than 1%
   of their bases in reads over 200 kb, all soil); no linkage statistic of any run has been computed. Added now, to be
   reported as exploratory: the class estimates without every run in which reads over 200 kb hold more than 1% of the
   bases. Where a table gives the read N50 of these runs, it is marked.
5. **2026-11-21 14:05 UTC. Outcome of the retry of note 2.** The retry ran at 08:19 UTC on 2026-11-21, when the first pass had
   ended (`logs/P4/P4_retry_DONE`). NCBI then held a record of ERR18104478 and ERR18104482 and none of the other
   six runs. All 8 fetches failed again ("empty block") and the 8 runs are rejected for good. The walk then tried
   every other unit of PRJEB104244 (58 units in all); each fetch failed in the same way. Three projects therefore
   have no accepted run: PRJEB104244 and PRJEB125553 (soil) and PRJEB123233 (human gut). No rule was changed.
6. **2026-11-21 14:05 UTC. Replacement runs in the four projects without a known kit.** After the 5 selected runs of PRJEB108219,
   PRJNA1469454 (human gut), PRJEB120797 (livestock) and PRJNA1444478 (wastewater) were rejected by the kit rule,
   the walk fetched the remaining units of these projects (4, 2, 1 and 4 runs). None has a kit pattern at the start
   of 5% of its reads (highest share 0.0008), and record and publication state neither kit nor chemistry, so these
   11 runs are rejected by the same rule (`public/screen/p4_manual_exclusions.tsv`, 52 runs). The four projects have
   no unit left. No rule was changed.
7. **2026-11-21 14:05 UTC. The sample of the confirmatory analysis.** No run is pending. `public/screen/p4_final.tsv` (sha256
   0f0d9b08d757…): 241 runs of 86 projects, all of them selected runs (no replacement was accepted): human gut 64
   runs / 19 projects, soil 63 / 22, freshwater 33 / 15, livestock 26 / 8, activated sludge 21 / 12, marine 18 / 10,
   wastewater 16 / 9. All 7 classes have at least 6 projects. Rejected: 52 kit rule (19 amplified, 33 kit and
   chemistry unknown, of which 11 were replacements), 62 fetch failed. Kit from the read ends of the 241 runs
   (`public/screen/p4_kit_probe.tsv`): ligation 76, native barcoding 61, rapid 10, not determinable 94 (kit or
   chemistry of these known from the record or the publication). No linkage statistic of any run has been computed.
