# External TMBWG inputs

Benchmark bytes are intentionally not versioned here. The 2DZP grids and CC0
CFL3D reference curves are fetched from TMBWG revision
`ddf1840b6e928d1de60d920cafe431ed476f36a1` and verified against
`data/benchmark-sources.json`:

```bash
python3 tools/fetch_benchmark_data.py --case flatplate
```

The command stores files under ignored `.cache/curvature-benchmarks/` and creates
ignored symlinks at the paths expected by the tutorial. See `THIRD_PARTY.md` for
licence and provenance details.
