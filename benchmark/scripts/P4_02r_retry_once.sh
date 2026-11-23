#!/usr/bin/env bash
# P4, note 2 of docs/preregistration_amendment_07_deviations.md: the one retry of the selected runs whose fetch failed
# while NCBI held no record of them. Waits until the driver log shows that round 1 has ended, then removes the
# P3_FAILED marker of those 8 runs, once, and records what NCBI says about each at that moment. The running driver
# picks them up at its next round. A second failure is final: this script never runs twice (marker logs/P4/P4_retry_DONE).
L=/mnt/rmrfc/workflow/iskander/linkstat; cd "$L"
[[ -e logs/P4/P4_retry_DONE ]] && { echo "retry already done"; exit 0; }
until grep -q "round 1 done" logs/P4/P4_02_driver.out; do sleep 60; done
{
  echo "retry at $(date -Is)"
  for r in ERR16963716 ERR16963732 ERR16963748 ERR16963751 ERR16963754 ERR18104478 ERR18104482 ERR18016248; do
    n=$(curl -s -m 30 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$r%5BAccession%5D" | grep -o '<Count>[0-9]*</Count>' | head -1)
    if [[ -e public/$r/P3_FAILED ]]; then rm -f "public/$r/P3_FAILED"; echo "$r marker removed; NCBI $n"; else echo "$r had no marker; NCBI $n"; fi
    sleep 1
  done
} > logs/P4/P4_retry_DONE
