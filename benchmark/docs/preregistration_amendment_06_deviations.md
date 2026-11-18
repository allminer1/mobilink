# Amendment 06: deviations and notes

2026-11-18. The frozen analysis of stage 2 was run once (`scripts/P3X_03_run.sh`, 2026-11-17 23:42 UTC).

## Deviations from the rules

None. The sample rule, the walk, the kit rule, the per-run analysis and the decision rule were applied as frozen.

## Notes

1. **Time.** The feasibility note expected about 15 min per run. The runs took 0.7 to 4 h each, because the server
   was oversubscribed by another project. The walk ran from 2026-11-17 11:46 to 23:42 UTC. No rule depends on time.
2. **Driver restarts.** The fetch driver was started three times (11:46, 22:25 and 23:42 UTC). Each restart followed
   a kit-rule rejection that was entered after the previous pass had ended, as the driver's header describes. Every
   run that was fetched is a run the one-by-one walk reaches (25 tried, 25 listed in `p3x_status.tsv`).
3. **Kit rule, judgement calls.**
   - ERR12160882: chemistry R9.4.1 taken from the name of the submitted file in the ENA record, as the frozen project
     table already lists it.
   - DRR333437: kit and flow cell from the methods of the linked publication (PMID 41731555).
   - SRR24197343: rejected; the ligation adapter is at 1.7% of read starts, below the 5% threshold of the rule (D63).
4. **SRR29786052.** The record states that the sample was enriched for bacteria with capture particles before DNA
   extraction. The eligibility rule excludes enrichment cultures and amplified libraries; this is neither, and the
   project is in the frozen list. The run was kept; the count without it is 10 of 10.
5. **Control strand.** One accepted run (SRR33225186) holds reads of the ONT lambda control strand, which the Kraken2
   database of the run script assigns to *Escherichia coli* (found on 2026-11-17 in step 3.4d, D78). The frozen
   analysis was not changed.
6. **Lines of the frozen script outside the scope of stage 2.** `verdict.txt` prints an H1b line ("not supported",
   p = 0.0535) and H3 lines. Amendment 06 does not test these again; they are reported as descriptive.
