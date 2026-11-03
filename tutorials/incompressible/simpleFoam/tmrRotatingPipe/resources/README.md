# Rotating-pipe external resources

The Zaets experimental profiles are fetched from the TMBWG repository pinned
at revision `ddf1840b6e928d1de60d920cafe431ed476f36a1` under CC0-1.0. They are
checksum-verified and materialised here as ignored symlinks:

```bash
python3 tools/fetch_benchmark_data.py --case rotatingpipe
```

No external grid is redistributed. The tutorial uses a project-authored
full-cylinder O-grid derived from the dimensions documented by TMR.
