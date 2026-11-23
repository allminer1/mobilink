#!/usr/bin/env bash
# P4 (amendment 07, Parts B to D): the analysis script on the accepted runs (public/screen/p4_final.tsv). Run once,
# after the walk has ended and the kit rule has been applied. Output tables/P4/confirm/ (verdict.txt holds the verdicts).
L=/mnt/rmrfc/workflow/iskander/linkstat; G=$L/frozen/amendment07; D=$L/envs/linkstat_dev/bin; cd "$L"
source scripts/P4_freeze_check.sh
[[ -s logs/P4/P4_02_DONE ]] || { echo "[4.5] the walk has not ended"; exit 1; }
[[ -z "$("$D/python" -I scripts/P4_02_walk.py)" ]] || { echo "[4.5] the walk has pending runs"; exit 1; }
[[ -e tables/P4/confirm ]] && { echo "[4.5] tables/P4/confirm exists; refusing to mix runs"; exit 1; }
mkdir -p tables/P4/confirm
cp public/screen/p4_final.tsv public/screen/p4_status.tsv tables/P4/confirm/
OMP_NUM_THREADS=1 nice -n 10 "$D/python" -I "$G/P4_03_confirm.py" public/screen/p4_final.tsv tables/P4/confirm --procs 12 \
  > logs/P4/P4_03_confirm.log 2>&1 || { echo "[4.5] analysis failed"; exit 1; }
date -Is > logs/P4/P4_03_DONE
