#!/usr/bin/env bash
# P4 step 4.9: the commands of package version 1.0 on one public run of Phase 3 (run on 2026-11-22 from the
# branch v1.0-dev in the worktree pkg_v1, commit d3d3c61; the branch is now merged into pkg, so the paths say pkg),
# as scripts/P3_06_runtime_one.sh did for version 0.5: prep (primary rule), chimeras, test (three registered pairs,
# Kraken2 species strata, unit loci, default seed, 9,999 draws). Output tmp/p49/<run>/.
# Nothing under tables/P3/confirm*, tables/P4/confirm or public/ is written; no frozen analysis is rerun.
# usage: P4_09_release_one.sh RUN
set -uo pipefail
L=/mnt/rmrfc/workflow/iskander/linkstat; D=$L/tmp/p49/bin; cd "$L"
run=$1; P=public/$run/$run; O=tmp/p49/$run; mkdir -p "$O/chim_tmp"
[[ -s $O/DONE ]] && exit 0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/usr/bin/time -v -o "$O/prep.time" nice -n 15 "$D/mobilink" prep --lengths "$P.lengths.tsv.gz" --card "$P.card.all.tsv.gz" \
  --bacmet "$P.bacmet.all.tsv.gz" --bacmet-classes pkg/resources/bacmet2_exp_classes.tsv \
  --mobileog "$P.mobileog.all.tsv.gz" --mobileog-lookup db/mobileog_all.lookup.tsv.gz \
  --mobility-core pkg/resources/mobility_core.tsv.gz --rule primary --out-prefix "$O/prep" > "$O/prep.log" 2>&1 || exit 1
/usr/bin/time -v -o "$O/chimeras.time" nice -n 15 "$D/mobilink" chimeras --paf "$P.ava.paf.gz" --tmp-dir "$O/chim_tmp" \
  --out "$O/chimeras.tsv" > "$O/chimeras.log" 2>&1 || exit 2
/usr/bin/time -v -o "$O/test.time" nice -n 15 "$D/mobilink" test --lengths "$P.lengths.tsv.gz" \
  --features "$O/prep.primary.features.tsv.gz" --pairs ARG:MGE --pairs ARG:MRG --pairs MRG:MGE \
  --kraken "$P.kraken.tsv.gz" --taxonomy db/kraken2/k2_standard_16gb_20241228/ktaxonomy.tsv --loci "$P.ava.paf.gz" \
  --out "$O/results.tsv" > "$O/test.log" 2>&1 || exit 3
rm -f "$O/chimeras.tsv"; rmdir "$O/chim_tmp" 2>/dev/null
date -Is > "$O/DONE"
