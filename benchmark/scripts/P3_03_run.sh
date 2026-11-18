#!/usr/bin/env bash
# P3 step 3.3 driver: the confirmatory analysis of amendment 05 on public/screen/p3_final.tsv, with the frozen script.
L=/mnt/rmrfc/workflow/iskander/linkstat; F=$L/frozen/amendment05; D=$L/envs/linkstat_dev/bin; cd "$L"
source scripts/P3_freeze_check.sh
[[ -e tables/P3/confirm ]] && { echo "[3.3] tables/P3/confirm exists; refusing to mix runs"; exit 1; }
mkdir -p tables/P3/confirm
cp public/screen/p3_final.tsv public/screen/p3_status.tsv tables/P3/confirm/
OMP_NUM_THREADS=1 nice -n 10 "$D/python" -I "$F/P3_03_confirm.py" public/screen/p3_final.tsv tables/P3/confirm --procs 12 \
  > logs/P3/P3_03_confirm.log 2>&1 || { echo "[3.3] analysis failed"; exit 1; }
date -Is > logs/P3/P3_03_DONE
