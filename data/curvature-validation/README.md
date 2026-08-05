# Generated validation data

This directory is intentionally empty in Git. Validation tables may reproduce
or interpolate external reference series, so they are generated locally or in
CI and uploaded as run artifacts instead of being redistributed in releases.

Fetch authorised benchmark inputs with `tools/fetch_benchmark_data.py`, run the
relevant OpenFOAM cases, then use `tools/collect_curvature_validation_data.py`.
Source licences and checksums are recorded in `data/benchmark-sources.json`.
