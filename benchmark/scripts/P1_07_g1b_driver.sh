#!/usr/bin/env bash
# P1 step 1.7: G1b data and gate (amendment 02). Requires public/screen/g1b_selected.tsv and the frozen Part B code.
# Per selected run: 3 Gbp random-block subsample, DIAMOND and features (P0_07_public_run.sh, as in P0), then a block
# check (reads per block = spots per block; the D14 lesson) written to public/<run>/block_check.tsv.
# Then the frozen Part B driver (mobility core = primary; D3 = secondary) and the frozen gate. Marker logs/P1_07_DONE.
set -uo pipefail
L=/mnt/rmrfc/workflow/iskander/linkstat; E=$L/envs/linkstat_p0/bin; cd "$L"
(cd frozen/partB && sha256sum -c --quiet ../partB.sha256) || { echo "[P1_07] frozen Part B files changed"; exit 1; }
echo "ef49bad8611769a17cdf37ac36245c2096404df79d80e7c74e25d50c7d962ec8  frozen/partB_gate_P1_08_g1b_gate.py" | sha256sum -c --quiet \
  || { echo "[P1_07] frozen gate script changed"; exit 1; }
runs=$(awk -F'\t' 'NR>1{print $3}' public/screen/g1b_selected.tsv)
for run in $runs; do
  spots=$(awk -F'\t' -v r="$run" '$1==r{print $2}' public/screen/g1b_blockprobe.tsv)
  mlen=$(awk -F'\t' -v r="$run" '$1==r{printf "%d", $5}' public/screen/g1b_blockprobe.tsv)
  scripts/P0_07_public_run.sh "$run" "$spots" "$mlen" 3000000000 >> logs/P1_07_driver.log 2>&1 \
    || { echo "[P1_07] $run FAILED $(date -Is)" >> logs/P1_07_driver.log; continue; }
  "$E/python" -I - "$run" <<'PY' >> logs/P1_07_driver.log 2>&1
import sys, numpy as np, pandas as pd
run = sys.argv[1]; L = "/mnt/rmrfc/workflow/iskander/linkstat"
b = pd.read_csv(f"{L}/public/{run}/blocks.tsv", sep="\t")
spot = pd.read_csv(f"{L}/public/{run}/{run}.lengths.tsv.gz", sep="\t", usecols=["read_id"]).read_id.str.rsplit(".", n=1).str[1].astype(np.int64)
blk = np.searchsorted(b.first_spot.to_numpy(), spot.to_numpy(), "right") - 1
b["expected"] = b.last_spot - b.first_spot + 1
b["got"] = np.bincount(blk[blk >= 0], minlength=len(b))
b.to_csv(f"{L}/public/{run}/block_check.tsv", sep="\t", index=False)
short = b[b.got < b.expected]
print(f"[P1_07] {run} block check: {len(short)} of {len(b)} blocks short; reads {b.got.sum()} of {b.expected.sum()}")
PY
done
ok=$(for r in $runs; do [[ -s public/$r/DONE ]] && echo "$r"; done)
nice -n 15 "$E/python" -I frozen/partB/P1_06_g1b.py --runs $ok --b-class MC --out tables/P1_07_g1b_MC.tsv > logs/P1_07_g1b_MC.log 2>&1
nice -n 15 "$E/python" -I frozen/partB/P1_06_g1b.py --runs $ok --b-class D3 --out tables/P1_07_g1b_D3.tsv > logs/P1_07_g1b_D3.log 2>&1
"$E/python" -I frozen/partB_gate_P1_08_g1b_gate.py --table tables/P1_07_g1b_MC.tsv --selected public/screen/g1b_selected.tsv \
  --out tables/P1_08_g1b_gate.tsv > logs/P1_08_g1b_gate.log 2>&1
date -Is > logs/P1_07_DONE
