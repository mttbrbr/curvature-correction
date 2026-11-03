# External TMBWG inputs

Benchmark bytes are intentionally not versioned here. The NACA 0012 grids and
CC0 CFL3D SA/SA-RC curves are fetched from TMBWG revision
`ddf1840b6e928d1de60d920cafe431ed476f36a1` and verified against
`data/benchmark-sources.json`:

```bash
python3 tools/fetch_benchmark_data.py --case naca0012
```

The command stores files under ignored `.cache/curvature-benchmarks/` and creates
ignored symlinks at the paths expected by the tutorial. Experimental and
digitised paper data remain excluded. See `THIRD_PARTY.md`.
