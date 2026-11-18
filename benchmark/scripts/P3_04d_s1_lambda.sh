#!/usr/bin/env bash
# P3 step 3.4d: find the reads of the ONT lambda control strand (DCS) in the low-coverage libraries of our first soil
# study by alignment to the lambda genome (db/lambda/NC_001416.1.fa, NCBI NC_001416.1).
# Reason: the Kraken2 database of the run script (standard, 16 GB) assigns these reads to Escherichia coli, not to
# Lambdavirus, so the Kraken2 taxon cannot be used to remove them (D78).
# Output public/<library>/<library>.lambda.paf.gz (minimap2 map-ont, all alignments to lambda). The rule that makes a
# read a control read is applied in scripts/P3_04d_s1_test.py.
set -euo pipefail
L=/mnt/rmrfc/workflow/iskander/linkstat; E=$L/envs/linkstat_p0/bin; cd "$L"
RAW=/mnt/rmrfc/data/oldpools/serverdata0/microbiome/raw; T=${T:-4}
awk -F'\t' 'NR>1 && $6=="analysed"{print $1"\t"$2}' tables/P3/P3_04d_s1_samples.tsv | while IFS=$'\t' read -r name file; do
  [[ -s public/$name/P3_DONE ]] || { echo "[3.4d lambda] $name is not annotated yet"; continue; }
  [[ -s public/$name/$name.lambda.paf.gz ]] && continue
  nice -n 15 "$E/minimap2" -x map-ont -c -t "$T" db/lambda/NC_001416.1.fa "$RAW/$file" 2> /dev/null < /dev/null \
    | cut -f1-12 | gzip > "public/$name/$name.lambda.paf.gz.tmp"
  mv "public/$name/$name.lambda.paf.gz.tmp" "public/$name/$name.lambda.paf.gz"
  echo "[3.4d lambda] $name $(zcat public/$name/$name.lambda.paf.gz | cut -f1 | sort -u | wc -l) reads with an alignment $(date -Is)"
done
n=$(ls public/EKR1_*/*.lambda.paf.gz 2>/dev/null | wc -l)
echo "[3.4d lambda] $n of 27 libraries $(date -Is)"
[[ $n -eq 27 ]] && date -Is > logs/P3/P3_04d_lambda_DONE
