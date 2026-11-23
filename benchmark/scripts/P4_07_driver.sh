#!/usr/bin/env bash
# P4 step 4.7 (exploratory, amendment 07 Part G): fetch and annotate the 40 runs of Go and Shrestha 2026 with the frozen
# run script of amendment 05, unchanged, 3 at a time (smallest run first). No eligibility or kit rule is applied: the
# sample is theirs. Runs: public/screen/p4_07_runs.tsv; order and counts: public/screen/p4_07_order.txt.
# A failed run gets public/<run>/P3_FAILED and is tried once more in a second round. Marker logs/P4/P4_07_DONE.
L=/mnt/rmrfc/workflow/iskander/linkstat; F=$L/frozen/amendment05; cd "$L"
mkdir -p logs/P4 logs/P3/runs tmp; rm -f logs/P4/P4_07_DONE
for round in 1 2; do
  : > tmp/p4_07_pending.txt
  while read -r run n b; do [[ -s public/$run/P3_DONE ]] || { rm -f public/$run/P3_FAILED; echo "$run $n $b" >> tmp/p4_07_pending.txt; }; done < public/screen/p4_07_order.txt
  [[ -s tmp/p4_07_pending.txt ]] || break
  echo "[4.7] round $round: $(wc -l < tmp/p4_07_pending.txt) runs $(date -Is)"
  xargs -P 3 -L 1 bash -c 'T=12 '"$F"'/P3_02_public_run.sh "$0" "$1" "$2" > logs/P3/runs/$0.out 2>&1 || { rm -rf tmp/p3/$0; mkdir -p public/$0; tail -1 logs/P3/runs/$0.out > public/$0/P3_FAILED; }' < tmp/p4_07_pending.txt
  echo "[4.7] round $round done $(date -Is)"
done
echo "done: $(ls public/*/P3_DONE | grep -c -F -f <(cut -d' ' -f1 public/screen/p4_07_order.txt)) of $(wc -l < public/screen/p4_07_order.txt)"
date -Is > logs/P4/P4_07_DONE
