#!/usr/bin/env bash
# Checks that amendment 07 is frozen and nothing frozen has changed: the amendment text, the checksummed files of
# frozen/amendment07/SHA256SUMS (scripts, selection tables, the frozen scripts of amendments 02 and 05 that are reused),
# and the package at the commit of frozen/amendment07/PACKAGE_COMMIT with a clean source tree. Sourced by the P4 drivers.
L=/mnt/rmrfc/workflow/iskander/linkstat; G=$L/frozen/amendment07
[[ -s $L/docs/preregistration_amendment_07.sha256 ]] || { echo "[P4] amendment 07 is not frozen"; exit 1; }
[[ "$(sha256sum $L/docs/preregistration_amendment_07.md | cut -c1-64)" == "$(cut -c1-64 $L/docs/preregistration_amendment_07.sha256)" ]] || { echo "[P4] amendment 07 text changed"; exit 1; }
(cd "$G" && sha256sum -c --quiet SHA256SUMS) || { echo "[P4] a frozen file of amendment 07 changed"; exit 1; }
[[ "$(git -C $L/pkg rev-parse HEAD)" == "$(cat $G/PACKAGE_COMMIT)" ]] || { echo "[P4] the package is not at the frozen commit"; exit 1; }
[[ -z "$(git -C $L/pkg status --porcelain -- src)" ]] || { echo "[P4] the package source has uncommitted changes"; exit 1; }
