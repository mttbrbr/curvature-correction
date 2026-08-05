#!/bin/sh
set -eu

: "${WM_PROJECT_DIR:?source OpenFOAM.com v2412 before running this test}"

repo_dir=$(CDPATH='' cd -- "$(dirname "$0")/.." && pwd)
test_dir="$repo_dir/tests/openfoam/mrfMappingRegression"
case_dir=$(mktemp -d /tmp/curvature-mrf-mapping.XXXXXX)
mpi_procs=${CURVATURE_MRF_MAPPING_MPI_PROCS:-2}

cleanup()
{
    if [ "${CURVATURE_MRF_MAPPING_KEEP:-0}" = 1 ]; then
        echo "MRF mapping case retained at $case_dir"
    else
        rm -rf -- "$case_dir"
    fi
}
trap cleanup EXIT HUP INT TERM

if [ "$mpi_procs" -ne 2 ]; then
    echo "ERROR: this four-cell regression requires exactly two MPI ranks" >&2
    exit 2
fi

(cd "$repo_dir/src/curvatureCorrectionModels" && wmake libso)
(cd "$test_dir" && wmake)

cp -a "$test_dir/case/." "$case_dir/"
blockMesh -case "$case_dir" > "$case_dir/log.blockMesh" 2>&1
topoSet -case "$case_dir" > "$case_dir/log.topoSet" 2>&1

set_scenario()
{
    scenario=$1
    rm -f "$case_dir/constant/MRFProperties" \
        "$case_dir/constant/MRFPropertiesSecond"

    if [ "$scenario" != no-mrf ]; then
        cp "$test_dir/scenarios/$scenario/MRFProperties" \
            "$case_dir/constant/MRFProperties"
    fi
    if [ "$scenario" = duplicate ]; then
        cp "$test_dir/scenarios/duplicate/MRFPropertiesSecond" \
            "$case_dir/constant/MRFPropertiesSecond"
    fi
}

run_success()
{
    scenario=$1
    cli_scenario=${2:-$scenario}
    set_scenario "$scenario"
    curvatureMRFMappingTest -case "$case_dir" -scenario "$cli_scenario" \
        > "$case_dir/log.$scenario.serial" 2>&1
    grep -q "MRF scenario $cli_scenario: pass" \
        "$case_dir/log.$scenario.serial"
}

run_parallel()
{
    scenario=$1
    set_scenario "$scenario"
    decomposePar -case "$case_dir" -force \
        > "$case_dir/log.$scenario.decomposePar" 2>&1
    mpirun -np "$mpi_procs" curvatureMRFMappingTest \
        -case "$case_dir" -scenario "$scenario" -parallel \
        > "$case_dir/log.$scenario.parallel" 2>&1
    grep -q "MRF scenario $scenario: pass" \
        "$case_dir/log.$scenario.parallel"
}

run_failure()
{
    scenario=$1
    cli_scenario=$2
    expected=$3
    set_scenario "$scenario"

    if curvatureMRFMappingTest -case "$case_dir" -scenario "$cli_scenario" \
        > "$case_dir/log.$scenario.expected-failure" 2>&1; then
        echo "ERROR: $cli_scenario unexpectedly succeeded" >&2
        exit 1
    fi

    grep -q "$expected" "$case_dir/log.$scenario.expected-failure"
}

run_success no-mrf
run_success mapping
run_success time
run_parallel mapping
run_parallel time
run_failure overlap overlap "belongs to overlapping active MRF zones"
run_failure duplicate duplicate-registry \
    "Expected at most one IOMRFZoneList in mesh registry"

echo "MRF mapping regression passed: empty, inactive, disjoint, timed, overlap, duplicate, serial and MPI"
