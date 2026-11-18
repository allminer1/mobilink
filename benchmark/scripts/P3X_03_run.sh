#!/usr/bin/env bash
# P3 extension (amendment 06, Parts B and C): the frozen stage-1 analysis script, unchanged, on the accepted stage-2
# runs (public/screen/p3x_final.tsv). Run once, after the walk has ended and the kit rule has been applied.
# Output tables/P3/confirm_ext/ (verdict.txt holds the stage-2 verdicts of H1a and H2; its H1b and H3 lines are
# descriptive and no H4 line is written, because there is no soil run).
L=/mnt/rmrfc/workflow/iskander/linkstat; F=$L/frozen/amendment05; D=$L/envs/linkstat_dev/bin; cd "$L"
source scripts/P3X_freeze_check.sh
[[ -s logs/P3X/P3X_02_DONE ]] || { echo "[3.3x] the walk has not ended"; exit 1; }
[[ -z "$("$D/python" -I scripts/P3X_02_walk.py)" ]] || { echo "[3.3x] the walk has pending runs"; exit 1; }
[[ -e tables/P3/confirm_ext ]] && { echo "[3.3x] tables/P3/confirm_ext exists; refusing to mix runs"; exit 1; }
mkdir -p tables/P3/confirm_ext
cp public/screen/p3x_final.tsv public/screen/p3x_status.tsv tables/P3/confirm_ext/
OMP_NUM_THREADS=1 nice -n 10 "$D/python" -I "$F/P3_03_confirm.py" public/screen/p3x_final.tsv tables/P3/confirm_ext --procs 12 \
  > logs/P3X/P3X_03_confirm.log 2>&1 || { echo "[3.3x] analysis failed"; exit 1; }
date -Is > logs/P3X/P3X_03_DONE
