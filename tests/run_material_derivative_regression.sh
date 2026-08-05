#!/bin/sh
set -eu

: "${WM_PROJECT_DIR:?source OpenFOAM.com v2412 before running this test}"

repo_dir=$(CDPATH='' cd -- "$(dirname "$0")/.." && pwd)
test_dir="$repo_dir/tests/openfoam/materialDerivativeRegression"
case_dir=$(mktemp -d /tmp/curvature-material-derivative.XXXXXX)

cleanup()
{
    if [ "${CURVATURE_MATERIAL_DERIVATIVE_KEEP:-0}" = 1 ]; then
        echo "Material-derivative case retained at $case_dir"
    else
        rm -rf -- "$case_dir"
    fi
}
trap cleanup EXIT HUP INT TERM

(cd "$repo_dir/src/curvatureCorrectionModels" && wmake libso)
(cd "$test_dir" && wmake)

cp -a "$test_dir/case/." "$case_dir/"
blockMesh -case "$case_dir" > "$case_dir/log.blockMesh" 2>&1

curvatureMaterialDerivativeTest -case "$case_dir" \
    > "$case_dir/log.transient" 2>&1

foamDictionary "$case_dir/system/fvSchemes" \
    -entry ddtSchemes.default -set steadyState >/dev/null
curvatureMaterialDerivativeTest -case "$case_dir" -expectSteady \
    > "$case_dir/log.steady" 2>&1

grep -q "Temporal contribution: pass" "$case_dir/log.transient"
grep -q "Convective contribution: pass" "$case_dir/log.transient"
grep -q "SST correction: pass" "$case_dir/log.transient"
grep -q "SA correction: pass" "$case_dir/log.transient"
grep -q "Steady temporal contribution: pass" "$case_dir/log.steady"

echo "Material-derivative regression passed for transient/steady SST-RC and SA-RC"
