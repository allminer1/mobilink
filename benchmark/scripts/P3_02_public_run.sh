#!/usr/bin/env bash
# P3 step 3.2 (amendment 05 Part A and B): one confirmatory public run -> tables; the reads are deleted at the end.
#   subsample: 3 Gbp as random blocks of consecutive spots (48 blocks, or 8 for runs held only at ENA: ERR/DRR), one block
#     per stratum of the run, seeded by the accession (as P0_07); the whole run if it has < 3.3 Gbp;
#   lengths and N50; DIAMOND blastx (CARD homolog, mobileOG-db, BacMet2 EXP) and P0_02c features (as P0_07);
#   registered classes (P2_05i: ARG, MRG primary; MGE mobility core);
#   Kraken2 Standard-16 on all reads at confidence 0 and 0.1 (as P2_00_kraken_v5);
#   carriers as in P2 step 2.0 -> minimap2 ava-ont (block >= 1,000);
#   adapter scan of all reads (ligation adapter, rapid primer, NB01-24 with flanks; P2_05e) for the kit and for the
#     chimera sensitivity analysis.
# Outputs public/<run>/ (layout of the earlier public runs, so the frozen Part B driver reads them). Marker P3_DONE.
# usage: P3_02_public_run.sh RUN N_SPOTS BASE_COUNT          (T = threads, default 12)
set -euo pipefail
L=/mnt/rmrfc/workflow/iskander/linkstat; S=$L/envs/linkstat_sra/bin; E=$L/envs/linkstat_p0/bin; D=$L/envs/linkstat_dev/bin; cd "$L"
K=$L/envs/linkstat_kraken/bin/kraken2; KDB=$L/db/kraken2/k2_standard_16gb_20241228
run=$1; n=$2; bases=$3; T=${T:-12}; target=3000000000
P=public/$run; W=tmp/p3/$run; mkdir -p "$P" "$W" logs/P3/runs
[[ -s $P/P3_DONE ]] && { echo "[3.2] $run done earlier"; exit 0; }
fa=$W/$run.fa
if [[ ! -s $fa ]]; then
  [[ -s $P/blocks.tsv ]] || python3 - "$run" "$n" "$bases" "$target" > "$P/blocks.tsv" <<'PY'
import random, sys, zlib
run, n, bases, target = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
B = 8 if run[:3] in ("ERR", "DRR") else 48
print("block\tfirst_spot\tlast_spot")
if bases < 1.1 * target:
    print(f"0\t1\t{n}")
else:
    k = max(200, int(target / (bases / n) / B))
    rng = random.Random(zlib.crc32(run.encode()))
    step = n // B
    for i in range(B):
        lo = i * step + 1; hi = max(lo, (i + 1) * step - k + 1)
        s = rng.randint(lo, hi); print(f"{i}\t{s}\t{min(n, s + k - 1)}")
PY
  t0=$(date +%s)
  tail -n +2 "$P/blocks.tsv" | xargs -P 8 -L 1 bash -c \
    'for t in 1 2 3; do timeout 14400 "'"$S"'/fastq-dump" -N "$1" -X "$2" --fasta 0 --stdout "'"$run"'" > "'"$W"'/b$0.fa" 2>/dev/null && break; sleep 60; done'
  for f in "$W"/b*.fa; do [[ -s $f ]] || { echo "[3.2] $run: empty block $f"; exit 2; }; done
  [[ $(ls "$W"/b*.fa | wc -l) -eq $(( $(wc -l < "$P/blocks.tsv") - 1 )) ]] || { echo "[3.2] $run: blocks missing"; exit 2; }
  cat $(ls "$W"/b*.fa | sort -V) > "$fa.tmp" && mv "$fa.tmp" "$fa" && rm -f "$W"/b*.fa
  echo -e "fetch_seconds\t$(( $(date +%s) - t0 ))" > "$P/fetch.tsv"
fi
"$E/seqkit" stats -a -T "$fa" > "$P/subsample_stats.tsv"
[[ -s $P/$run.lengths.tsv.gz ]] || { printf 'read_id\tlength\n'; "$E/seqkit" fx2tab -n -i -l "$fa"; } | "$E/pigz" -p 2 > "$P/$run.lengths.tsv.gz"
for pair in card:$L/db/card_homolog mobileog:$L/db/mobileog_all bacmet:/mnt/rmrfc/workflow/iskander/vko8soil/db/bacmet/bacmet2_exp; do
  tag=${pair%%:*}; db=${pair#*:}; out=$P/$run.$tag.all.tsv.gz
  [[ -s $out ]] && continue
  nice -n 15 "$E/diamond" blastx -q "$fa" -d "$db" -F 15 --range-culling --top 10 -p "$T" \
    -f 6 qseqid sseqid pident length qstart qend sstart send slen qlen evalue bitscore --tmpdir tmp \
    -o "$W/$tag.all.tsv" > "logs/P3/runs/$run.$tag.log" 2>&1
  "$E/pigz" -p 4 -c "$W/$tag.all.tsv" > "$out.tmp" && mv "$out.tmp" "$out" && rm "$W/$tag.all.tsv"
done
[[ -s $P/$run.primary.features.tsv.gz ]] || nice -n 15 "$E/python" -I scripts/P0_02c_features.py "$run" \
  --arg "$P/$run.card.all.tsv.gz" --mrg "$P/$run.bacmet.all.tsv.gz" --mge "$P/$run.mobileog.all.tsv.gz" \
  --lengths "$P/$run.lengths.tsv.gz" --out-dir "$P" > "logs/P3/runs/$run.features.log" 2>&1
[[ -s $P/$run.core.features.tsv.gz ]] || "$D/python" -I scripts/P2_05i_core_features.py "$P/$run.primary.features.tsv.gz" \
  "$P/$run.mobileog.all.tsv.gz" "$P/$run.core.features.tsv.gz" > "logs/P3/runs/$run.core.log" 2>&1
if [[ ! -s $P/$run.kraken.tsv.gz ]]; then
  nr=$(grep -c '^>' "$fa")
  for c in 0 0.1; do
    nice -n 15 "$K" --db "$KDB" --threads "$T" --confidence $c --report "$W/kraken_c$c.report" --output "$W/kraken_c$c.out" \
      "$fa" 2> "logs/P3/runs/$run.kraken_c$c.log"
    [[ "$(wc -l < "$W/kraken_c$c.out")" -eq "$nr" ]] || { echo "[3.2] $run: Kraken2 (confidence $c) rows differ from reads"; exit 3; }
  done
  cmp -s <(cut -f2 "$W/kraken_c0.out") <(cut -f2 "$W/kraken_c0.1.out") || { echo "[3.2] $run: Kraken2 read order differs"; exit 3; }
  paste <(cut -f1-3,5 "$W/kraken_c0.out") <(cut -f3 "$W/kraken_c0.1.out") | awk -F'\t' -v OFS='\t' '
    BEGIN{print "read_id","status","taxid","hit_kmers","taxon_kmers","taxid_c01"}
    {h=0; x=0; n=split($4, p, " "); for(i=1;i<=n;i++){split(p[i], q, ":"); if(q[1]!="0" && q[1]!="A"){h+=q[2]; if(q[1]==$3) x+=q[2]}}
     print $2,$1,$3,h,x,$5}' | "$E/pigz" -p 4 > "$W/kraken.tsv.gz"
  mv "$W/kraken_c0.report" "$P/$run.kraken_c0.report"; mv "$W/kraken_c0.1.report" "$P/$run.kraken_c01.report"
  mv "$W/kraken.tsv.gz" "$P/$run.kraken.tsv.gz"
fi
{ for r in primary relaxed strict mgeperm; do zcat "$P/$run.$r.features.tsv.gz" | awk -F'\t' -v r=$r 'NR>1{print $1"\t"r}'; done
  zcat "$P/$run.card.all.tsv.gz" | awk -F'\t' '$3>=40{print $1"\tcard40"}'
  zcat "$P/$run.bacmet.all.tsv.gz" | awk -F'\t' '$3>=40{print $1"\tbacmet40"}'; } | sort -u > "$W/carriers.tsv"
cut -f1 "$W/carriers.tsv" | sort -u > "$W/ids.txt"
"$E/seqkit" grep -f "$W/ids.txt" "$fa" > "$W/carriers.fa"
[[ -s $P/$run.p2_carriers.tsv.gz ]] || "$E/pigz" -p 4 -c "$W/carriers.tsv" > "$P/$run.p2_carriers.tsv.gz"
if [[ ! -s $P/$run.ava.paf.gz ]]; then
  nice -n 15 minimap2 -x ava-ont -t "$T" "$W/carriers.fa" "$W/carriers.fa" 2> "logs/P3/runs/$run.ava.log" \
    | awk -F'\t' -v OFS='\t' '$11>=1000{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | "$E/pigz" -p 4 > "$W/ava.paf.gz"
  mv "$W/ava.paf.gz" "$P/$run.ava.paf.gz"
fi
[[ -s $P/$run.adapters.tsv.gz ]] || nice -n 15 "$D/python" -I scripts/P2_05e_adapter_scan.py "$fa" \
  tables/P2/mocks/ont_barcodes_dorado.tsv "$P/$run.adapters.tsv.gz" "$T" all > "logs/P3/runs/$run.adapters.log" 2>&1
echo -e "minimap2\t$(minimap2 --version)\nkraken2\tk2_standard_16gb_20241228\tconfidence 0 and 0.1, all reads\ndiamond\t$("$E/diamond" version | awk '{print $3}')" > "$P/p3_versions.tsv"
rm -rf "$W"
date -Is > "$P/P3_DONE"; echo "[3.2] $run done $(date -Is)"
