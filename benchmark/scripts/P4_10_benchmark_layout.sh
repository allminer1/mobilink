#!/usr/bin/env bash
# Phase 5: adds Phase 4 (amendment 07, the survey of 241 public runs, and the exploratory analyses of the paper)
# to the benchmark folder of the package repository (pkg/benchmark), which scripts/P3_06_benchmark_layout.sh built.
# Files are copied unchanged and the project layout is kept, so frozen/amendment07/SHA256SUMS still verifies.
# Not copied: reads and per-run tables (public/<run>/), databases, the raw archive records
# (public/screen/p4_records/), the third-party table of Go and Shrestha (external/) and the read-level copy made
# from it (tables/P4/reanalysis/their_reads.tsv.gz), and two large or path-bearing check files of tables/P4/v1check.
# Rerunnable: only the Phase 4 parts are replaced. README.md is written by hand.
set -euo pipefail
L=/mnt/rmrfc/workflow/iskander/linkstat; B=$L/pkg/benchmark; cd "$L"
# registration text, checksum, deviations
cp docs/preregistration_amendment_07.md docs/preregistration_amendment_07.sha256 docs/preregistration_amendment_07_deviations.md "$B/docs/"
# frozen folder
rm -rf "$B/frozen/amendment07"; mkdir -p "$B/frozen/amendment07"
find frozen/amendment07 -maxdepth 1 -type f -exec cp {} "$B/frozen/amendment07/" \;
# scripts: those of the checksum list, the drivers, and the exploratory analyses
cp scripts/P4_01_candidates.py scripts/P4_01b_records.py scripts/P4_01c_project_table.py scripts/P4_01d_units.py \
   scripts/P4_02_walk.py scripts/P4_02_driver.sh scripts/P4_02b_record_kits.py scripts/P4_02c_kits.py \
   scripts/P4_02r_retry_once.sh scripts/P4_02_design_dev.py scripts/P4_03_confirm.py scripts/P4_03_run.sh \
   scripts/P4_freeze_check.sh scripts/P4_04_sens_samples.py scripts/P4_04b_sens_variants.py \
   scripts/P4_06a_subclass_counts.py scripts/P4_06b_subclass_linkage.py scripts/P4_06c_subclass_summary.py \
   scripts/P4_07_driver.sh scripts/P4_07a_their_calls.py scripts/P4_07b_driver.sh scripts/P4_07c_estimates.py \
   scripts/P4_07d_kit_slopes.py scripts/P4_07e_compare.py scripts/P4_08_figures.py scripts/P4_09_release.sh \
   scripts/P4_09_release_one.sh scripts/P4_09_release_table.py scripts/P4_09_supp_tables.py \
   scripts/P4_10_benchmark_layout.sh "$B/scripts/"
# selection tables
cp public/screen/ena_p4_runs_2026-11-18.tsv public/screen/p4_*.tsv public/screen/p4_07_order.txt public/screen/p4_07b_order.txt "$B/public/screen/"
# result tables
rm -rf "$B/tables/P4"; mkdir -p "$B"/tables/P4/{confirm,sens,sub,reanalysis,kits,v1check}
cp tables/P4/P4_02_design_dev.*.tsv "$B/tables/P4/"
for d in confirm sens sub kits; do find "tables/P4/$d" -maxdepth 1 -type f -exec cp {} "$B/tables/P4/$d/" \; ; done
find tables/P4/reanalysis -maxdepth 1 -type f ! -name 'their_reads.tsv.gz' -exec cp {} "$B/tables/P4/reanalysis/" \;
cp tables/P4/v1check/P4_09_*.tsv tables/P4/v1check/meta_ARG_MGE.{groups,jackknife,contrasts}.tsv "$B/tables/P4/v1check/"
# reports
mkdir -p "$B/reports/P4"; cp reports/P4/*.md "$B/reports/P4/"
find "$B" -name __pycache__ -prune -exec rm -rf {} +
(cd "$B/frozen/amendment07" && sha256sum -c --quiet SHA256SUMS) && echo "amendment 07 checksums verify inside the benchmark folder"
(cd "$B" && sha256sum -c --quiet docs/preregistration_amendment_07.sha256) && echo "amendment 07 text verifies"
du -sh "$B"; find "$B" -type f | wc -l
