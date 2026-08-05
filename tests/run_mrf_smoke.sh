#!/bin/sh
set -eu

: "${WM_PROJECT_DIR:?source OpenFOAM.com v2412 before running this test}"

repo_dir=$(CDPATH='' cd -- "$(dirname "$0")/.." && pwd)
template_dir="$WM_PROJECT_DIR/tutorials/incompressible/simpleFoam/rotatingCylinders"

if [ ! -d "$template_dir" ]; then
    echo "ERROR: rotatingCylinders tutorial not found at $template_dir" >&2
    exit 2
fi

case_dir=$(mktemp -d /tmp/curvature-mrf-smoke.XXXXXX)
cleanup()
{
    if [ "${CURVATURE_MRF_KEEP:-0}" = 1 ]; then
        echo "MRF smoke case retained at $case_dir"
    else
        rm -rf -- "$case_dir"
    fi
}
trap cleanup EXIT HUP INT TERM

cp -a "$template_dir/." "$case_dir/"
foamDictionary "$case_dir/system/controlDict" \
    -entry libs -set '("libcurvatureCorrectionModels.so")' >/dev/null
foamDictionary "$case_dir/system/controlDict" -entry endTime -set 2 >/dev/null
foamDictionary "$case_dir/system/controlDict" -entry writeInterval -set 1 >/dev/null
foamDictionary "$case_dir/system/controlDict" -entry functions -remove >/dev/null

cp "$repo_dir/tutorials/incompressible/simpleFoam/curvedChannel/constant/turbulenceProperties" \
    "$case_dir/constant/turbulenceProperties"
foamDictionary "$case_dir/constant/turbulenceProperties" \
    -entry RAS/kOmegaSSTCCCoeffs/writeCurvatureDiagnostics -set on >/dev/null
cp "$repo_dir/tutorials/incompressible/simpleFoam/curvedChannel/0.orig/k" "$case_dir/0/k"
cp "$repo_dir/tutorials/incompressible/simpleFoam/curvedChannel/0.orig/omega" "$case_dir/0/omega"
cp "$repo_dir/tutorials/incompressible/simpleFoam/curvedChannel/0.orig/nut" "$case_dir/0/nut"

sed -i \
    -e 's/    inlet { type fixedValue; value uniform 0.375; }/    innerWall { type kqRWallFunction; value uniform 0.375; }/' \
    -e 's/    outlet { type inletOutlet; inletValue uniform 0.375; value uniform 0.375; }/    outerWall { type kqRWallFunction; value uniform 0.375; }/' \
    -e '/    walls { type kqRWallFunction; value uniform 0.375; }/d' \
    "$case_dir/0/k"
sed -i \
    -e 's/    inlet { type fixedValue; value uniform 100; }/    innerWall { type omegaWallFunction; value uniform 100; }/' \
    -e 's/    outlet { type inletOutlet; inletValue uniform 100; value uniform 100; }/    outerWall { type omegaWallFunction; value uniform 100; }/' \
    -e '/    walls { type omegaWallFunction; value uniform 100; }/d' \
    "$case_dir/0/omega"
sed -i \
    -e 's/    inlet { type calculated; value uniform 0; }/    innerWall { type nutkWallFunction; value uniform 0; }/' \
    -e 's/    outlet { type calculated; value uniform 0; }/    outerWall { type nutkWallFunction; value uniform 0; }/' \
    -e '/    walls { type nutkWallFunction; value uniform 0; }/d' \
    "$case_dir/0/nut"

cp "$repo_dir/tutorials/incompressible/simpleFoam/curvedChannel/system/fvSchemes" \
    "$case_dir/system/fvSchemes"
cp "$repo_dir/tutorials/incompressible/simpleFoam/curvedChannel/system/fvSolution" \
    "$case_dir/system/fvSolution"
foamDictionary "$case_dir/system/fvSolution" -entry SIMPLE.pRefCell -set 0 >/dev/null
foamDictionary "$case_dir/system/fvSolution" -entry SIMPLE.pRefValue -set 0 >/dev/null

blockMesh -case "$case_dir" > "$case_dir/log.blockMesh" 2>&1
simpleFoam -case "$case_dir" > "$case_dir/log.simpleFoam" 2>&1

mpi_procs=${CURVATURE_MRF_MPI_PROCS:-2}
cp "$repo_dir/tutorials/incompressible/simpleFoam/smitsCurve/baseCase/system/decomposeParDict" \
    "$case_dir/system/decomposeParDict"
foamDictionary "$case_dir/system/decomposeParDict" \
    -entry numberOfSubdomains -set "$mpi_procs" >/dev/null
foamDictionary "$case_dir/system/controlDict" -entry startFrom -set latestTime >/dev/null
foamDictionary "$case_dir/system/controlDict" -entry endTime -set 3 >/dev/null
decomposePar -case "$case_dir" -latestTime -force > "$case_dir/log.decomposePar" 2>&1
mpirun -np "$mpi_procs" simpleFoam -case "$case_dir" -parallel \
    > "$case_dir/log.simpleFoam-mpi" 2>&1
reconstructPar -case "$case_dir" -latestTime > "$case_dir/log.reconstructPar" 2>&1

postProcess -case "$case_dir" -latestTime -func "fieldMinMax(OmegaRot)" \
    > "$case_dir/log.OmegaRot" 2>&1

grep -q "min(mag(OmegaRot)) = 100" "$case_dir/log.OmegaRot"
grep -q "max(mag(OmegaRot)) = 100" "$case_dir/log.OmegaRot"

diagnostic_fields="curvatureIndicator productionBase productionCorrected productionRatio fr1Clipped"
for field in $diagnostic_fields; do
    test -s "$case_dir/3/$field"
done
if grep -Eiq "(^|[^[:alpha:]])(nan|inf)([^[:alpha:]]|$)" \
    "$case_dir/log.simpleFoam" "$case_dir/log.simpleFoam-mpi" \
    "$case_dir/3/fr1" "$case_dir/3/OmegaRot" \
    $(for field in $diagnostic_fields; do printf "%s " "$case_dir/3/$field"; done); then
    echo "ERROR: non-finite value found in SST MRF diagnostics" >&2
    exit 1
fi

echo "Restarting the same MRF case with SpalartAllmarasCC"
cp "$case_dir/3/k" "$case_dir/3/nuTilda"
foamDictionary "$case_dir/3/nuTilda" -entry FoamFile/object -set nuTilda >/dev/null
foamDictionary "$case_dir/3/nuTilda" -entry dimensions -set "[0 2 -1 0 0 0 0]" >/dev/null
foamDictionary "$case_dir/3/nuTilda" -entry boundaryField/innerWall \
    -set "{ type fixedValue; value uniform 0; }" >/dev/null
foamDictionary "$case_dir/3/nuTilda" -entry boundaryField/outerWall \
    -set "{ type fixedValue; value uniform 0; }" >/dev/null
foamDictionary "$case_dir/constant/turbulenceProperties" \
    -entry RAS/RASModel -set SpalartAllmarasCC >/dev/null
foamDictionary "$case_dir/constant/turbulenceProperties" \
    -entry RAS/SpalartAllmarasCCCoeffs \
    -set "{ curvatureCorrection on; writeCurvatureDiagnostics on; Cr1 1; Cr2 12; Cr3 1; }" >/dev/null
foamDictionary "$case_dir/system/fvSchemes" \
    -entry "divSchemes/div(phi,nuTilda)" -set "Gauss upwind" >/dev/null
foamDictionary "$case_dir/system/fvSolution" -entry solvers/nuTilda \
    -set "{ solver smoothSolver; smoother GaussSeidel; tolerance 1e-7; relTol 0.1; }" >/dev/null
foamDictionary "$case_dir/system/fvSolution" \
    -entry relaxationFactors/equations/nuTilda -set 0.7 >/dev/null
foamDictionary "$case_dir/system/controlDict" -entry endTime -set 4 >/dev/null
simpleFoam -case "$case_dir" > "$case_dir/log.simpleFoam-sa" 2>&1
for field in fr1 OmegaRot $diagnostic_fields; do
    test -s "$case_dir/4/$field"
done
if grep -Eiq "(^|[^[:alpha:]])(nan|inf)([^[:alpha:]]|$)" \
    "$case_dir/log.simpleFoam-sa" "$case_dir/4/fr1" "$case_dir/4/OmegaRot" \
    $(for field in $diagnostic_fields; do printf "%s " "$case_dir/4/$field"; done); then
    echo "ERROR: non-finite value found in SA MRF diagnostics" >&2
    exit 1
fi

echo "MRF smoke test passed for SST serial/MPI and SA serial diagnostics: OmegaRot=100 rad/s"
