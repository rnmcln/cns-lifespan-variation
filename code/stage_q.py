"""Stage Q. Rank agreement between reference life tables, and a literally common open interval.

Two questions the earlier stages left open.

1. How far do the national and common references disagree about the ordering of countries?
   Stage I counted the countries moving two or more places, which depends on an arbitrary
   cut. This stage reports rank correlations, the distribution of absolute displacement and
   the number of pairwise country orderings that reverse, none of which requires a cut.

2. The GBD 2023 reference used in the primary analysis is common at every age except the
   open interval, where the four sub-interval values are averaged with each country-year's
   own distribution of all-cause deaths above 95 years. The reference is therefore not
   literally common there. This stage recomputes the comparison with one fixed open-interval
   value applied to every country-year, under two choices of that value, and reports how far
   the results move.

Inputs : gbd/tmrlt_gbd2023.csv, hmd_all.csv, eurostat_full.csv, results_countryyear.csv,
         official_reference_country.csv, raw/*_Deaths_5x1.txt
Outputs: RANKS.json, rank_agreement.csv, open_interval_common.csv, hmd_open_interval_deaths.csv
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

from engine import ESP, KEYS, build_lt, e_within  # noqa: F401

OPEN_COLS = ["d95_99", "d100_104", "d105_109", "d110p"]
SUB_AGES = np.array([97.5, 102.5, 107.5, 110.0])

O = pd.read_csv("gbd/tmrlt_gbd2023.csv")
oa = O["Age"].values.astype(float)
oe = O["Life Expectancy"].values
e_star = lambda a: np.interp(a, oa, oe)          # noqa: E731
SUB_E = e_star(SUB_AGES)
STARTS = np.array([0, 1] + list(range(5, 95, 5)), float)
WID = np.array([1.0, 4.0] + [5.0] * 18)

# ----------------------------------------------------------------- 1. rank agreement
g = pd.read_csv("official_reference_country.csv").set_index("geo")


def agreement(a, b, label_a, label_b):
    """Rank agreement between two rate vectors, ranked from highest rate downwards."""
    ra = pd.Series(a).rank(ascending=False)
    rb = pd.Series(b).rank(ascending=False)
    d = (rb - ra).abs()
    rho = stats.spearmanr(a, b).statistic
    tau = stats.kendalltau(a, b)
    n = len(a)
    n_pairs = n * (n - 1) // 2
    # discordant pairs: orderings that reverse between the two references
    disc = 0
    av, bv = np.asarray(a, float), np.asarray(b, float)
    for i in range(n):
        for j in range(i + 1, n):
            if (av[i] - av[j]) * (bv[i] - bv[j]) < 0:
                disc += 1
    return {
        "comparison": "%s vs %s" % (label_a, label_b),
        "n": n,
        "spearman_rho": round(float(rho), 3),
        "kendall_tau_b": round(float(tau.statistic), 3),
        "median_abs_shift": float(d.median()),
        "p90_abs_shift": float(d.quantile(0.9)),
        "max_abs_shift": int(d.max()),
        "n_shift_ge2": int(d.ge(2).sum()),
        "n_pairs": n_pairs,
        "n_pairs_reversed": disc,
        "pct_pairs_reversed": round(100.0 * disc / n_pairs, 1),
    }


rows = [
    agreement(g.asyr_nat, g.asyr_off, "national", "GBD 2023"),
    agreement(g.asyr_nat, g.asyr_bp, "national", "minimum-mortality"),
    agreement(g.asyr_off, g.asyr_bp, "GBD 2023", "minimum-mortality"),
]
RA = pd.DataFrame(rows)
RA.to_csv("rank_agreement.csv", index=False)

# ------------------------------------------- 2. literally common open interval
H = pd.read_csv("hmd_all.csv")
E = pd.read_csv("eurostat_full.csv")
Rr = pd.read_csv("results_countryyear.csv")
MAPc = {"AT": "AUT", "BE": "BEL", "BG": "BGR", "HR": "HRV", "CZ": "CZE", "DK": "DNK",
        "EE": "EST", "FI": "FIN", "FR": "FRATNP", "DE": "DEUTNP", "EL": "GRC", "HU": "HUN",
        "IS": "ISL", "IE": "IRL", "IT": "ITA", "LV": "LVA", "LT": "LTU", "LU": "LUX",
        "NL": "NLD", "NO": "NOR", "PL": "POL", "PT": "PRT", "SK": "SVK", "SI": "SVN",
        "ES": "ESP", "SE": "SWE", "CH": "CHE", "UK": "GBR_NP"}
OPEN_AGES = ["95-99", "100-104", "105-109", "110+"]


def read_hmd(path):
    L = open(path).read().split("\n")
    h = [i for i, l in enumerate(L) if l.strip().startswith("Year")][0]
    cols = L[h].split()
    df = pd.DataFrame([l.split() for l in L[h + 1:] if l.strip()], columns=cols)
    df["Year"] = df.Year.astype(int)
    for c in cols[2:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


# All-cause deaths above 95 years, by country-year. Derived here from the same raw series the
# rest of the pipeline uses, so the shipped copy is an output rather than an input.
wrows = []
for geo, cc in sorted(MAPc.items()):
    d = read_hmd("raw/%s_Deaths_5x1.txt" % cc)
    d = d[d.Age.isin(OPEN_AGES)]
    for yr, sgrp in d.groupby("Year"):
        if yr < 2012 or yr > 2021:
            continue
        v = sgrp.set_index("Age")["Total"].reindex(OPEN_AGES).fillna(0).values.astype(float)
        wrows.append([geo, int(yr)] + list(v))
W = pd.DataFrame(wrows, columns=["geo", "Year"] + OPEN_COLS)
W.to_csv("hmd_open_interval_deaths.csv", index=False)

cns = E[E.which == "cns"].rename(columns={"time": "Year"})
cns = cns[cns.k.notna()]
piv = cns.pivot_table(index=["geo", "Year", "sex", "k"], values="deaths", dropna=False).reset_index()
ok = (Rr[(Rr.sex == "T") & (Rr.year.between(2012, 2018))][["geo", "sex", "year"]]
      .rename(columns={"year": "Year"}))
M = (H.merge(piv, on=["geo", "Year", "sex", "k"], how="inner")
      .merge(ok, on=["geo", "sex", "Year"])
      .sort_values(["geo", "Year", "ord"]))

wf = W[W.Year.between(2012, 2018)].set_index(["geo", "Year"])[OPEN_COLS]
OWc = {k: (v.values / v.values.sum() if v.values.sum() > 0 else np.array([1.0, 0, 0, 0]))
       for k, v in wf.iterrows()}

# pooled weighting: one distribution of deaths above 95, formed over the whole frame
pooled = wf.sum().values
pooled = pooled / pooled.sum()
E_POOLED = float((pooled * SUB_E).sum())
E_FIXED = float(e_star(97.5))            # the single sub-interval midpoint

rows = []
for (geo, yr), s in M.groupby(["geo", "Year"]):
    s = s.sort_values("ord")
    d = np.nan_to_num(s.deaths.values)
    ax = s.ax.values
    expo = s.expo.values
    if d.sum() == 0:
        continue
    base = e_star(STARTS + np.clip(ax[:20] / WID, 0, 1) * WID)

    def summarise(e_open):
        eo = np.concatenate([base, [e_open]])
        ypd = float((d * eo).sum() / d.sum())
        asyr = float((np.divide(d * eo, expo, out=np.zeros_like(d), where=expo > 0) * ESP).sum() * 1e5)
        return ypd, asyr

    e_cy = float((OWc[(geo, yr)] * SUB_E).sum()) if (geo, yr) in OWc else E_POOLED
    y_cy, a_cy = summarise(e_cy)
    y_po, a_po = summarise(E_POOLED)
    y_fx, a_fx = summarise(E_FIXED)
    rows.append(dict(geo=geo, year=yr, deaths=d.sum(), d_open=d[20],
                     e_open_countryyear=e_cy, e_open_pooled=E_POOLED, e_open_fixed=E_FIXED,
                     ypd_countryyear=y_cy, ypd_pooled=y_po, ypd_fixed=y_fx,
                     asyr_countryyear=a_cy, asyr_pooled=a_po, asyr_fixed=a_fx))

D = pd.DataFrame(rows)
D.to_csv("open_interval_common.csv", index=False)

cty = D.groupby("geo").apply(lambda z: pd.Series({
    "ypd_countryyear": (z.ypd_countryyear * z.deaths).sum() / z.deaths.sum(),
    "ypd_pooled": (z.ypd_pooled * z.deaths).sum() / z.deaths.sum(),
    "ypd_fixed": (z.ypd_fixed * z.deaths).sum() / z.deaths.sum(),
    "asyr_countryyear": z.asyr_countryyear.mean(),
    "asyr_pooled": z.asyr_pooled.mean(),
    "asyr_fixed": z.asyr_fixed.mean(),
}), include_groups=False)
cty["ypd_nat"] = g.ypd_nat
cty["asyr_nat"] = g.asyr_nat

J = {"rank_agreement": {r["comparison"]: r for r in RA.to_dict("records")}}

J["open_share_deaths_pct"] = round(float(100 * D.d_open.sum() / D.deaths.sum()), 2)
J["e_open_countryyear_lo"] = round(float(D.e_open_countryyear.min()), 3)
J["e_open_countryyear_hi"] = round(float(D.e_open_countryyear.max()), 3)
J["e_open_pooled"] = round(E_POOLED, 3)
J["e_open_fixed"] = round(E_FIXED, 3)
J["ypd_maxabs_diff_pooled"] = round(float((cty.ypd_pooled - cty.ypd_countryyear).abs().max()), 4)
J["ypd_maxabs_diff_fixed"] = round(float((cty.ypd_fixed - cty.ypd_countryyear).abs().max()), 4)
J["asyr_maxabs_diff_pooled"] = round(float((cty.asyr_pooled - cty.asyr_countryyear).abs().max()), 3)
J["asyr_maxabs_diff_fixed"] = round(float((cty.asyr_fixed - cty.asyr_countryyear).abs().max()), 3)
for tag in ["pooled", "fixed"]:
    pct = 100 * (cty["ypd_" + tag] / cty.ypd_nat - 1)
    J["off_pct_mean_" + tag] = round(float(pct.mean()), 1)
    J["off_pct_lo_" + tag] = round(float(pct.min()), 1)
    J["off_pct_hi_" + tag] = round(float(pct.max()), 1)
    a = agreement(cty.asyr_nat, cty["asyr_" + tag], "national", "GBD 2023 " + tag)
    J["rank_shift_ge2_" + tag] = a["n_shift_ge2"]
    J["rank_max_shift_" + tag] = a["max_abs_shift"]
    J["rho_" + tag] = a["spearman_rho"]
    J["tau_" + tag] = a["kendall_tau_b"]
    J["rho_vs_countryyear_" + tag] = round(
        float(stats.spearmanr(cty["asyr_" + tag], cty.asyr_countryyear).statistic), 4)
    J["rank_identical_" + tag] = bool(
        (cty["asyr_" + tag].rank(ascending=False) == cty.asyr_countryyear.rank(ascending=False)).all())

json.dump(J, open("RANKS.json", "w"), indent=1)
print(RA.to_string(index=False))
print()
print(json.dumps(J, indent=1))
