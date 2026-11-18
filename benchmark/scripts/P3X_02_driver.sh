#!/usr/bin/env bash
# P3 extension (amendment 06): fetch and annotate the stage-2 candidates in the order of the walk, 3 at a time, with the
# frozen stage-1 run script. Each round fetches the pending runs of the walk; it ends when the walk has no pending run.
# The kit rule is then applied by hand (public/screen/p3x_manual_exclusions.tsv); if a run is excluded, start this
# driver again. Marker logs/P3X/P3X_02_DONE (written each time the walk has no pending run).
L=/mnt/rmrfc/workflow/iskander/linkstat; F=$L/frozen/amendment05; D=$L/envs/linkstat_dev/bin; cd "$L"
source scripts/P3X_freeze_check.sh
mkdir -p logs/P3X logs/P3/runs; rm -f logs/P3X/P3X_02_DONE
for round in $(seq 1 20); do
  "$D/python" -I scripts/P3X_02_walk.py > tmp/p3x_pending.txt
  [[ -s tmp/p3x_pending.txt ]] || break
  echo "[3.2x] round $round: $(wc -l < tmp/p3x_pending.txt) runs $(date -Is)"
  xargs -P 3 -L 1 bash -c 'T=12 '"$F"'/P3_02_public_run.sh "$0" "$1" "$2" > logs/P3/runs/$0.out 2>&1 || { rm -rf tmp/p3/$0; mkdir -p public/$0; tail -1 logs/P3/runs/$0.out > public/$0/P3_FAILED; }' < tmp/p3x_pending.txt
  echo "[3.2x] round $round done $(date -Is)"
done
"$D/python" -I scripts/P3X_02_walk.py > /dev/null
cut -f2,7 public/screen/p3x_status.tsv | sort | uniq -c
date -Is > logs/P3X/P3X_02_DONE
