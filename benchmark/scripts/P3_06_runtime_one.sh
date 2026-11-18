#!/usr/bin/env bash
# P3 step 3.5 (V8): CPU time and peak memory of the released commands on one public run, from the tables in
# public/<run>/ (no reads are fetched). prep (primary rule), chimeras (overlap breakpoints) and test (the three
# registered pairs, Kraken2 species strata, unit loci from the carrier overlaps; default seed and 9,999 draws, as the
# frozen analyses). Output tmp/p35/<run>/ : *.time (GNU time -v), results.tsv, prep.primary.features.tsv.gz.
# The frozen analyses are not rerun: nothing under tables/P3/confirm* or public/ is written.
# usage: P3_06_runtime_one.sh RUN
set -uo pipefail
L=/mnt/rmrfc/workflow/iskander/linkstat; D=$L/envs/linkstat_dev/bin; V=/mnt/rmrfc/workflow/iskander/vko8soil; cd "$L"
run=$1; P=public/$run/$run; O=tmp/p35/$run; mkdir -p "$O"
[[ -s $O/DONE ]] && exit 0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/usr/bin/time -v -o "$O/prep.time" nice -n 15 "$D/mobilink" prep --lengths "$P.lengths.tsv.gz" --card "$P.card.all.tsv.gz" \
  --bacmet "$P.bacmet.all.tsv.gz" --bacmet-classes pkg/resources/bacmet2_exp_classes.tsv \
  --mobileog "$P.mobileog.all.tsv.gz" --mobileog-lookup db/mobileog_all.lookup.tsv.gz \
  --mobility-core pkg/resources/mobility_core.tsv.gz --rule primary --out-prefix "$O/prep" > "$O/prep.log" 2>&1 || exit 1
/usr/bin/time -v -o "$O/chimeras.time" nice -n 15 "$D/mobilink" chimeras --paf "$P.ava.paf.gz" --out "$O/chimeras.tsv" \
  > "$O/chimeras.log" 2>&1 || exit 2
/usr/bin/time -v -o "$O/test.time" nice -n 15 "$D/mobilink" test --lengths "$P.lengths.tsv.gz" \
  --features "$O/prep.primary.features.tsv.gz" --pairs ARG:MGE --pairs ARG:MRG --pairs MRG:MGE \
  --kraken "$P.kraken.tsv.gz" --taxonomy db/kraken2/k2_standard_16gb_20241228/ktaxonomy.tsv --loci "$P.ava.paf.gz" \
  --out "$O/results.tsv" > "$O/test.log" 2>&1 || exit 3
rm -f "$O/chimeras.tsv"
date -Is > "$O/DONE"
