#!/usr/bin/env bash
# P3 step 3.5 (V8): the released commands on every accepted public run of stage 1 (tables/P3/confirm/p3_final.tsv)
# and stage 2 (public/screen/p3x_final.tsv), 4 runs at a time, one thread each, nice 15.
# Then scripts/P3_06_runtime_table.py -> tables/P3/P3_06_runtime.tsv and tables/P3/P3_06_release_check.tsv.
L=/mnt/rmrfc/workflow/iskander/linkstat; D=$L/envs/linkstat_dev/bin; cd "$L"
mkdir -p tmp/p35 logs/P3
{ tail -n +2 tables/P3/confirm/p3_final.tsv | cut -f1; tail -n +2 public/screen/p3x_final.tsv | cut -f1; } > tmp/p35/runs.txt
echo "[3.5] $(wc -l < tmp/p35/runs.txt) runs $(date -Is)"
xargs -P 4 -L 1 bash -c 'scripts/P3_06_runtime_one.sh "$0" || echo "[3.5] $0 failed at step $?"' < tmp/p35/runs.txt
echo "[3.5] runs done $(date -Is)"
"$D/python" -I scripts/P3_06_runtime_table.py && date -Is > logs/P3/P3_06_DONE
