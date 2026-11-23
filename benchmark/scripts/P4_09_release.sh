#!/usr/bin/env bash
# P4 step 4.9: release check and run time of package version 1.0 on the 64 accepted public runs of Phase 3
# (the run list of scripts/P3_06_runtime.sh), 6 runs at a time, one thread each, nice 15.
# Then scripts/P4_09_release_table.py -> tables/P4/v1check/P4_09_{runtime,runtime.summary,release_check}.tsv.
L=/mnt/rmrfc/workflow/iskander/linkstat; cd "$L"
mkdir -p tmp/p49 logs/P4 tables/P4/v1check
{ tail -n +2 tables/P3/confirm/p3_final.tsv | cut -f1; tail -n +2 public/screen/p3x_final.tsv | cut -f1; } > tmp/p49/runs.txt
echo "[4.9] $(wc -l < tmp/p49/runs.txt) runs, branch commit $(git -C pkg rev-parse --short HEAD) $(date -Is)"
xargs -P 6 -L 1 bash -c 'scripts/P4_09_release_one.sh "$0" || echo "[4.9] $0 failed at step $?"' < tmp/p49/runs.txt
echo "[4.9] runs done $(date -Is)"
envs/linkstat_dev/bin/python -I scripts/P4_09_release_table.py && date -Is > logs/P4/P4_09_release_DONE
