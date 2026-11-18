#!/usr/bin/env bash
# Checks that amendment 06 is frozen and nothing frozen has changed: the amendment 05 check (text, scripts, package
# commit, clean package source), the amendment 06 text, and the checksummed files of frozen/amendment06/SHA256SUMS.
# Sourced by the P3X drivers.
L=/mnt/rmrfc/workflow/iskander/linkstat; G=$L/frozen/amendment06
source "$L/scripts/P3_freeze_check.sh"
[[ -s $L/docs/preregistration_amendment_06.sha256 ]] || { echo "[P3X] amendment 06 is not frozen"; exit 1; }
[[ "$(sha256sum $L/docs/preregistration_amendment_06.md | cut -c1-64)" == "$(cut -c1-64 $L/docs/preregistration_amendment_06.sha256)" ]] || { echo "[P3X] amendment 06 text changed"; exit 1; }
(cd "$G" && sha256sum -c --quiet SHA256SUMS) || { echo "[P3X] a frozen file of amendment 06 changed"; exit 1; }
