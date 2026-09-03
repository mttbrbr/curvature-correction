# Generated validation reports

Reviewed validation evidence is stored under one immutable campaign directory
per run. Schema v2 records provenance and computes acceptance states from
observations and preregistered thresholds. Failed or unverified critical checks
make `--strict-acceptance` fail.

Raw third-party reference profiles are not committed. Each accepted campaign
contains Markdown, HTML, aggregate JSON, permitted plots, and a copy plus
SHA-256 of the private archive manifest. The private workflow opens a pull
request; it never pushes evidence directly to `main`.
