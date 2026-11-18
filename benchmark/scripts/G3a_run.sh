#!/usr/bin/env bash
# Gate G3 (amendment 04 Part 3, frozen 2026-11-16, D46). Checks, before any G3 seed is used: the sha256 of
# docs/preregistration_amendment_04.md, the frozen scripts, read-length profiles and this driver
# (frozen/amendment04/SHA256SUMS), the package commit (frozen/amendment04/PACKAGE_COMMIT) and a clean package source.
# G3a: 18 settings x {no linkage: 500, linked 0.03: 200, linked 0.1: 200} datasets with the frozen
# P2_07_feasibility.py (net_linked_test with unit-locus clumps from positions), G3 seeds 90,000,000-90,017,899 (seed
# scheme 2: 90,000,000 + 1,000 setting + {0, 500, 700}[linked] + dataset), density cap 2.5, 2,000 hosts at 0.1 x;
# 32 processes at nice 15 (the process count does not change any result). Then the G3a summary (rules R1-R4) and the
# G3b simulation summary. Outputs tmp/G3a/s*.tsv and tables/G3/; logs logs/G3/; markers logs/G3a_DONE (datasets)
# and logs/G3_DONE (summaries). Refuses to start if tmp/G3a exists (no mixing of runs).
set -uo pipefail
L=/mnt/rmrfc/workflow/iskander/linkstat; E=$L/envs/linkstat_dev/bin; F=$L/frozen/amendment04; cd "$L"
[[ -s docs/preregistration_amendment_04.sha256 ]] || { echo "[G3] amendment 04 is not frozen"; exit 1; }
[[ "$(sha256sum docs/preregistration_amendment_04.md | cut -c1-64)" == "$(cut -c1-64 docs/preregistration_amendment_04.sha256)" ]] \
  || { echo "[G3] amendment 04 text differs from its frozen sha256"; exit 1; }
(cd "$F" && sha256sum -c --quiet SHA256SUMS) || { echo "[G3] frozen scripts or inputs changed"; exit 1; }
want=$(cat "$F/PACKAGE_COMMIT"); got=$(git -C pkg rev-parse HEAD)
[[ "$want" == "$got" ]] || { echo "[G3] package at $got, frozen commit is $want"; exit 1; }
[[ -z "$(git -C pkg status --porcelain -- src)" ]] || { echo "[G3] package source has uncommitted changes"; exit 1; }
[[ -e tmp/G3a ]] && { echo "[G3] tmp/G3a exists; refusing to mix runs"; exit 1; }
mkdir -p tmp/G3a logs/G3 tables/G3
echo "[G3] checks passed $(date -Is); package $got"
export OMP_NUM_THREADS=1 P2_SEED0=90000000 P2_SEED_SCHEME=2 P2_HOSTS01=2000 P2_DENSITY_CAP=2.5
{ for k in $(seq 0 17); do for c in $(seq 0 19); do echo "$k 0 $((c*25)) $(((c+1)*25))"; done
  for lk in 0.03 0.1; do for c in $(seq 0 7); do echo "$k $lk $((c*25)) $(((c+1)*25))"; done; done; done; } | xargs -P 32 -L 1 bash -c \
  'nice -n 15 '"$E"'/python -I '"$F"'/P2_07_feasibility.py $0 $1 $2 $3 tmp/G3a/s$0_l$1_$2.tsv > logs/G3/G3a_s$0_l$1_$2.log 2>&1 || echo "$0 $1 $2" >> logs/G3/G3a_failed'
date -Is > logs/G3a_DONE
[[ -s logs/G3/G3a_failed ]] && { echo "[G3] G3a chunks failed (logs/G3/G3a_failed); no summary"; exit 1; }
"$E/python" -I "$F/P2_07j_g3a_summary.py" tmp/G3a tables/G3/G3a_settings.tsv > tables/G3/G3a_verdict.txt 2> logs/G3/G3a_summary.log
"$E/python" -I "$F/G3b_summary.py" tmp/G3a tables/G3/G3b_settings.tsv > tables/G3/G3b_verdict.txt 2> logs/G3/G3b_summary.log
cat tables/G3/G3a_verdict.txt tables/G3/G3b_verdict.txt
date -Is > logs/G3_DONE
