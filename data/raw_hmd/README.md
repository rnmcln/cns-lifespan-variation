# Human Mortality Database source files

This directory is intentionally empty in the public repository.

The analysis uses five 5x1 series per country (`Deaths`, `Exposures`, `bltper`, `fltper`,
`mltper`) for 28 populations, 140 files in total. Human Mortality Database terms of use
prohibit redistribution, so the files are not included here.

To obtain them, register free of charge at https://www.mortality.org and, in the working
directory described in the top-level README, create `.hmdenv` containing

```
HMD_USER=<your registered email>
HMD_PASS=<your password>
```

and run `python3 fetch_hmd.py`, which writes the files into `raw/` in the working directory.

`data/derived/hmd_file_manifest.csv` records the series vintage, file size and
full SHA-256 digest of every file used, so a retrieval can be checked against the versions
that produced the reported results. Vintages range from November 2021 to April 2026, and
results are vintage-dependent.
