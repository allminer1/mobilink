#!/usr/bin/env bash
# P1 step 1.15: run gate G2 (amendment 03 v2) with the frozen scripts in frozen/amendment03/.
# Checks the frozen checksums and the package commit, builds the G1b feature tables (prep; CARD + BacMet2 EXP +
# mobileOG with the mobility core), runs G2a, runs G2b (6 settings x 4,000 datasets, 16 processes at nice 15:
# each setting split into chunks), then the gate. Marker logs/P1_15_G2_DONE.
set -uo pipefail
L=/mnt/rmrfc/workflow/iskander/linkstat; E=$L/envs/linkstat_dev/bin; F=$L/frozen/amendment03; cd "$L"
(cd "$F" && sha256sum -c --quiet SHA256SUMS) || { echo "[G2] frozen scripts changed"; exit 1; }
want=$(cat "$F/PACKAGE_COMMIT"); got=$(git -C pkg rev-parse HEAD)
[[ "$want" == "$got" ]] || { echo "[G2] package at $got, frozen commit is $want"; exit 1; }
[[ -z "$(git -C pkg status --porcelain -- src)" ]] || { echo "[G2] package source has uncommitted changes"; exit 1; }
mkdir -p tmp/g2 tables/G2
export OMP_NUM_THREADS=1
BC=/mnt/rmrfc/workflow/iskander/vko8soil/db/bacmet/bacmet_exp_classes.tsv
for run in $(awk -F'\t' 'NR>1{print $3}' public/screen/g1b_selected.tsv); do
  P=public/$run/$run
  [[ -s tmp/g2/$run.primary.features.tsv.gz ]] || nice -n 15 "$E/linkstat" prep --lengths $P.lengths.tsv.gz \
    --card $P.card.all.tsv.gz --bacmet $P.bacmet.all.tsv.gz --bacmet-classes $BC --mobileog $P.mobileog.all.tsv.gz \
    --mobileog-lookup db/mobileog_all.lookup.tsv.gz --mobility-core frozen/partB/mobility_core.tsv.gz \
    --rule primary --out-prefix tmp/g2/$run > logs/G2_prep_$run.log 2>&1
done
nice -n 15 "$E/python" -I "$F/P1_15_g2a.py" tables/G2/g2a.tsv > logs/G2_g2a.log 2>&1
# G2b: 6 settings x 4,000 datasets in chunks; the heavy settings (3-5, 800 genomes) get more chunks
jobs=()
for k in 0 1 2; do for c in 0 1; do jobs+=("$k $((c*2000)) $(((c+1)*2000))"); done; done
for k in 3 4 5; do for c in 0 1 2 3; do jobs+=("$k $((c*1000)) $(((c+1)*1000))"); done; done
printf '%s\n' "${jobs[@]}" | xargs -P 16 -L 1 bash -c \
  'nice -n 15 '"$E"'/python -I '"$F"'/P1_15_g2b.py $0 $1 $2 tables/G2/g2b_s$0_$1.tsv > logs/G2_g2b_s$0_$1.log 2>&1'
"$E/python" -I "$F/P1_15_g2_gate.py" --g2a tables/G2/g2a.tsv --g2b tables/G2/g2b_s*.tsv --out tables/G2/G2 \
  > logs/G2_gate.log 2>&1
date -Is > logs/P1_15_G2_DONE
