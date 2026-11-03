# External TMBWG inputs

Benchmark bytes are intentionally not versioned here. The quick Smits grid and
CC0 CFL3D reference curves are fetched from TMBWG revision
`ddf1840b6e928d1de60d920cafe431ed476f36a1` and verified against
`data/benchmark-sources.json`:

```bash
python3 tools/fetch_benchmark_data.py --case smits
```

The command stores files under ignored `.cache/curvature-benchmarks/` and creates
ignored symlinks at the paths expected by the tutorial. The experimental Smits
series and NASA medium/fine archive are on licence hold and are never fetched.
For `Allrun-gridStudy`, supply lawfully obtained grids with
`SMITS_MEDIUM_GRID` and `SMITS_FINE_GRID`. See `THIRD_PARTY.md`.
