# Third-party notices

## OpenFOAM

The adapter targets OpenFOAM.com v2412 and links against GPL-covered OpenFOAM
libraries. `SpalartAllmarasCC` and `kOmegaSSTCC` contain code derived from the
corresponding OpenFOAM implementations. Their source files retain upstream
copyright notices and are distributed under `GPL-3.0-or-later`.

OpenFOAM itself is a build/runtime dependency and is not bundled. Official
licensing information is available at <https://www.openfoam.com/licensing>.

This offering is not approved or endorsed by OpenCFD Limited, producer and
distributor of the OpenFOAM software via www.openfoam.com, and owner of the
OPENFOAM® and OpenCFD® trade marks.

OPENFOAM® is a registered trade mark of OpenCFD Limited, producer and
distributor of the OpenFOAM software via <https://www.openfoam.com/>.

## Publications

The implementation uses model equations described by Shur et al.
(<https://doi.org/10.2514/2.1058>) and Smirnov and Menter
(<https://doi.org/10.1115/1.3070573>). The publications are cited, not copied or
redistributed.

## TMBWG validation data

The optional validation tooling references the TMBWG Turbulence Modeling
Resource repository at revision
`ddf1840b6e928d1de60d920cafe431ed476f36a1`, distributed under CC0-1.0:
<https://github.com/TMBWG/turbmodels>.

External inputs are not committed. Their URLs, revisions, checksums and status
are recorded in `data/benchmark-sources.json`. The fetcher only materializes
assets whose redistribution status is recorded as `downloadable`.

The Smits experimental series and medium/fine grid archive remain on hold
because file-level redistribution rights were not verified. The optional
Appelbaum--Kloker--Wenzel dataset is CC BY 4.0 and would require attribution,
version and modification notices before adoption.

Runnable OpenFOAM tutorial definitions and automation are distributed under
the repository GPL-3.0-or-later licence. Downloaded benchmark inputs and
generated simulation results are excluded from version control.

## Xiao rotating-pipe reference

`data/reference-models.json` contains transcriptions of Table 1 and equations
3.2/3.4 from M. Xiao, A. Ceci, P. Orlandi and S. Pirozzoli, Journal of Fluid
Mechanics 996 (2024), A24, <https://doi.org/10.1017/jfm.2024.811>. These reference
values retain their CC BY 4.0 licence
(<https://creativecommons.org/licenses/by/4.0/>); source, authors and
transcription details are recorded in the file. No raw DNS data is bundled.
