# Changelog

All notable changes follow semantic versioning.

## Unreleased

- Add scientifically explicit runtime names `kOmegaSST2003RC` and
  `SpalartAllmarasNoFt2RC` with RC enabled by default.
- Retain deprecated `kOmegaSSTCC` and `SpalartAllmarasCC` aliases through 1.x.
- Make MRF support explicitly experimental.
- Introduce computed schema-v2 validation gates and a private-runner workflow.
- Apply GPL-3.0-or-later repository-wide and add OPENFOAM® trade-mark notices.
- Remove runnable tutorial cases from the public repository.

The first tagged release from this line will be `v0.1.0`. Version `v1.0.0` is
reserved for a campaign satisfying every critical gate in
`docs/validation-protocol.md`.
