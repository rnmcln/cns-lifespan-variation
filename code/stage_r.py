"""Stage R. Eurostat disclosure control within the analysed frame.

stage_f.py counts confidentiality-flagged C70 to C72 age cells across the whole Eurostat
extract, which includes countries outside the analysis (Cyprus, Liechtenstein, Malta) and
years outside 2012 to 2018. This stage restricts the count to the 28 countries and the
primary period actually analysed, so that the figures reported for the frame describe the
frame. Runs from the shipped files alone.

Writes DISCLOSURE.json.
"""
import json

import pandas as pd

E = pd.read_csv("eurostat_full.csv")
R = pd.read_csv("results_countryyear.csv")
FRAME = sorted(R.geo.unique())
YEARS = (2012, 2018)

f = E[(E.which == "cns") & E.flag.notna() & (E.flag != "") & E.deaths.isna()]
f = f[f.geo.isin(FRAME) & f.time.between(*YEARS)]

YOUNG = {"0", "1-4", "5-9", "10-14", "15-19"}
out = {"frame_countries": len(FRAME), "frame_years": list(YEARS)}
for sx in ("T", "M", "F"):
    s = f[f.sex == sx]
    out["flagged_" + sx] = int(len(s))
    out["flagged_countries_" + sx] = sorted(s.geo.unique())
    out["flagged_below20_pct_" + sx] = round(100 * float(s.k.isin(YOUNG).mean()), 1) if len(s) else None
out["flagged_countries_sexspec_extra"] = sorted(
    (set(out["flagged_countries_M"]) | set(out["flagged_countries_F"]))
    - set(out["flagged_countries_T"]))

json.dump(out, open("DISCLOSURE.json", "w"), indent=1)
print(json.dumps(out, indent=1))
