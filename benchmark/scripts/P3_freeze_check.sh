#!/usr/bin/env bash
# Checks that amendment 05 is frozen and nothing frozen has changed: the amendment text, the checksummed scripts and
# inputs (frozen/amendment05/SHA256SUMS), the package commit and a clean package source. Sourced by the P3 drivers.
L=/mnt/rmrfc/workflow/iskander/linkstat; F=$L/frozen/amendment05
[[ -s $L/docs/preregistration_amendment_05.sha256 ]] || { echo "[P3] amendment 05 is not frozen"; exit 1; }
[[ "$(sha256sum $L/docs/preregistration_amendment_05.md | cut -c1-64)" == "$(cut -c1-64 $L/docs/preregistration_amendment_05.sha256)" ]] || { echo "[P3] amendment 05 text changed"; exit 1; }
(cd "$F" && sha256sum -c --quiet SHA256SUMS) || { echo "[P3] a frozen file changed"; exit 1; }
[[ "$(git -C $L/pkg rev-parse HEAD)" == "$(cat $F/PACKAGE_COMMIT)" ]] || { echo "[P3] package commit differs from the frozen one"; exit 1; }
[[ -z "$(git -C $L/pkg status --porcelain -- src)" ]] || { echo "[P3] package source has uncommitted changes"; exit 1; }
