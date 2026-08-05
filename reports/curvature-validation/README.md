# Generated validation reports

Reports, metrics, and plots are generated locally or in CI and are not committed.
This avoids silently relicensing embedded external reference series. Use
`tools/generate_curvature_validation_report.py --strict-acceptance` for a release
gate: critical checks that are failed or not verified make the command fail.

Generated report directories should be published as CI artifacts together with
the manifest revision and workflow run identifier.
