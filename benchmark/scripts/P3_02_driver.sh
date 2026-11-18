#!/usr/bin/env bash
# P3 step 3.2 driver: fetch and annotate the confirmatory runs of amendment 05 (frozen), 3 at a time, with the automatic
# replacements of Part A (failed fetch, subsample < 1 Gbp, N50 < 2 kb), up to 6 rounds. Uses the frozen copies of the
# run script and the check script. Marker logs/P3/P3_02_DONE.
L=/mnt/rmrfc/workflow/iskander/linkstat; F=$L/frozen/amendment05; D=$L/envs/linkstat_dev/bin; cd "$L"
source scripts/P3_freeze_check.sh
mkdir -p logs/P3/runs
for round in 1 2 3 4 5 6; do
  "$D/python" -I "$F/P3_02b_check.py" > /dev/null
  tail -n +2 public/screen/p3_runlist.tsv | awk -F'\t' '{print $1, $5, $6}' | while read r n b; do
    [[ -s public/$r/P3_DONE || -e public/$r/P3_FAILED ]] || echo "$r $n $b"; done > tmp/p3_pending.txt
  [[ -s tmp/p3_pending.txt ]] || break
  echo "[3.2] round $round: $(wc -l < tmp/p3_pending.txt) runs $(date -Is)"
  xargs -P 3 -L 1 bash -c 'T=12 '"$F"'/P3_02_public_run.sh "$0" "$1" "$2" > logs/P3/runs/$0.out 2>&1 || { rm -rf tmp/p3/$0; tail -1 logs/P3/runs/$0.out > public/$0/P3_FAILED; }' < tmp/p3_pending.txt
  added=$("$D/python" -I "$F/P3_02b_check.py")
  echo "[3.2] round $round done: $added replacement(s) added $(date -Is)"
  [[ "$added" == "0" ]] && break
done
"$D/python" -I "$F/P3_02b_check.py" > /dev/null
date -Is > logs/P3/P3_02_DONE
