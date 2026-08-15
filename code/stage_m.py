"""Stage M: file manifest with complete digests.

Rebuilds the Human Mortality Database file manifest with full 64-character
SHA-256 digests for every retrieved file, including the Deaths_5x1 series used to weight
the open interval, which the earlier manifest omitted.

Writes additional_file_1_hmd_manifest.csv.
"""
import pandas as pd, numpy as np, hashlib, glob, os, re

MAPc = {"AT": "AUT", "BE": "BEL", "BG": "BGR", "HR": "HRV", "CZ": "CZE", "DK": "DNK",
        "EE": "EST", "FI": "FIN", "FR": "FRATNP", "DE": "DEUTNP", "EL": "GRC",
        "HU": "HUN", "IS": "ISL", "IE": "IRL", "IT": "ITA", "LV": "LVA", "LT": "LTU",
        "LU": "LUX", "NL": "NLD", "NO": "NOR", "PL": "POL", "PT": "PRT", "SK": "SVK",
        "SI": "SVN", "ES": "ESP", "SE": "SWE", "CH": "CHE", "UK": "GBR_NP"}
INV = {v: k for k, v in MAPc.items()}

# ---------------------------------------------------------------- 1. manifest
rows = []
for p in sorted(glob.glob("raw/*_5x1.txt")):
    b = open(p, "rb").read()
    rows.append(dict(file=os.path.basename(p), bytes=len(b),
                     sha256=hashlib.sha256(b).hexdigest()))
dg = pd.DataFrame(rows)
if dg.empty:
    raise SystemExit("no Human Mortality Database files in raw/; run fetch_hmd.py first")
vin = pd.read_csv("hmd_vintages.csv")
dg["hmd_code"] = dg.file.str.replace(r"_(Deaths|Exposures|bltper|fltper|mltper)_5x1\.txt$",
                                     "", regex=True)
dg["series"] = dg.file.str.extract(r"_(Deaths|Exposures|bltper|fltper|mltper)_5x1\.txt$")
dg["geo"] = dg.hmd_code.map(INV)
dg["series_key"] = dg.series + "_5x1"
vint = vin.set_index("hmd_code").series_vintage.to_dict()
dg["vintage"] = dg.hmd_code.map(vint)
M = dg[["geo", "hmd_code", "series_key", "bytes", "vintage", "sha256"]].rename(
    columns={"series_key": "file", "vintage": "series_vintage"}).sort_values(
    ["geo", "file"]).reset_index(drop=True)
assert M.sha256.str.len().eq(64).all(), "digest length"
per = M.groupby("hmd_code").size()
incomplete = per[per != 5]
if len(incomplete):
    raise SystemExit("incomplete series for %s; expected five files per country, run "
                     "fetch_hmd.py" % ", ".join(incomplete.index))
if len(per) != 28:
    raise SystemExit("found %d of 28 populations in raw/; run fetch_hmd.py" % len(per))
M.to_csv("hmd_file_manifest.csv", index=False)

print("manifest: %d files across %d populations, all digests 64 characters"
      % (len(M), len(per)))
