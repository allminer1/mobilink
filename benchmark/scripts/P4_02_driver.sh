#!/usr/bin/env bash
# P4 (amendment 07, Parts A and B): fetch and annotate the runs of the sample, 4 at a time, with the frozen run script
# of amendment 05. Each round fetches the pending runs of the walk (the sample, then replacements of rejected runs); it
# ends when the walk has no pending run. The kit rule is then applied by hand (public/screen/p4_manual_exclusions.tsv);
# if a run is excluded, start this driver again. Marker logs/P4/P4_02_DONE (written each time the walk has no pending run).
L=/mnt/rmrfc/workflow/iskander/linkstat; F=$L/frozen/amendment05; D=$L/envs/linkstat_dev/bin; cd "$L"
source scripts/P4_freeze_check.sh
mkdir -p logs/P4 logs/P3/runs tmp; rm -f logs/P4/P4_02_DONE
for round in $(seq 1 30); do
  "$D/python" -I scripts/P4_02_walk.py > tmp/p4_pending.txt
  [[ -s tmp/p4_pending.txt ]] || break
  echo "[4.4] round $round: $(wc -l < tmp/p4_pending.txt) runs $(date -Is)"
  xargs -P 4 -L 1 bash -c 'T=12 '"$F"'/P3_02_public_run.sh "$0" "$1" "$2" > logs/P3/runs/$0.out 2>&1 || { rm -rf tmp/p3/$0; mkdir -p public/$0; tail -1 logs/P3/runs/$0.out > public/$0/P3_FAILED; }' < tmp/p4_pending.txt
  echo "[4.4] round $round done $(date -Is)"
done
"$D/python" -I scripts/P4_02_walk.py > /dev/null
cut -f2,7 public/screen/p4_status.tsv | sort | uniq -c
date -Is > logs/P4/P4_02_DONE
