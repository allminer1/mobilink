#!/usr/bin/env bash
# P3 step 3.5, exploratory sensitivity (outside the registration; no frozen table is touched): the registered per-run
# test of the runs with a control-strand signal (scripts/P3_06c_control_strand.py), repeated without the suspected
# control reads. Suspected control read: Kraken2 taxid 562 and length 3,400-3,619 (the reads were not kept, so no
# alignment is possible). Same features, strata, overlaps, seed and draws as the registered per-run analysis.
# Output tmp/p35c/<run>/results.tsv; table by scripts/P3_06c_control_table.py.
L=/mnt/rmrfc/workflow/iskander/linkstat; D=$L/envs/linkstat_dev/bin; cd "$L"
export OMP_NUM_THREADS=1
one() {
  run=$1; P=public/$run/$run; O=tmp/p35c/$run; mkdir -p "$O"
  [[ -s $O/results.tsv ]] && return 0
  join -t $'\t' <(zcat "$P.kraken.tsv.gz" | awk -F'\t' 'NR>1 && $3==562{print $1}' | sort) \
                <(zcat "$P.lengths.tsv.gz" | awk -F'\t' 'NR>1 && $2>=3400 && $2<3620{print $1}' | sort) > "$O/control_suspects.txt"
  nice -n 15 "$D/mobilink" test --lengths "$P.lengths.tsv.gz" --features "$P.core.features.tsv.gz" \
    --pairs ARG:MGE --pairs ARG:MRG --pairs MRG:MGE --kraken "$P.kraken.tsv.gz" \
    --taxonomy db/kraken2/k2_standard_16gb_20241228/ktaxonomy.tsv --loci "$P.ava.paf.gz" \
    --exclude "$O/control_suspects.txt" --no-oe --out "$O/results.tsv" > "$O/test.log" 2>&1
}
export -f one; export L D
printf '%s\n' SRR23278716 SRR23278685 SRR23278728 SRR23278753 ERR14825711 ERR14825617 SRR33225186 | xargs -P 4 -I{} bash -c 'one {}'
date -Is > logs/P3/P3_06c_DONE
