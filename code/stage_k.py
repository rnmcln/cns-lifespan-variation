"""Stage K: additional robustness quantities requested at revision.

Adds two results that the previous version described but did not quantify:

  1. Balanced panel. The primary frame excluding the two countries with
     incomplete series (Iceland, Luxembourg), so that every retained country
     contributes all seven years.
  2. Later period. The countries with complete coverage from 2012 to 2021,
     summarised over that longer window and compared with their own values in
     the primary frame.

Aggregation follows stage_d.py exactly: exposure-weighted rates, unweighted
means of life-table quantities across years within a country, death-weighted
mean age at death and share below the threshold.

Writes SENS3.json.
"""
import pandas as pd, numpy as np, json
from scipy import stats

R = pd.read_csv("results_countryyear.csv")
P = json.load(open("PANEL.json"))
NM = {"AT": "Austria", "BE": "Belgium", "BG": "Bulgaria", "HR": "Croatia", "CZ": "Czechia",
      "DK": "Denmark", "EE": "Estonia", "FI": "Finland", "FR": "France", "DE": "Germany",
      "EL": "Greece", "HU": "Hungary", "IS": "Iceland", "IE": "Ireland", "IT": "Italy",
      "LV": "Latvia", "LT": "Lithuania", "LU": "Luxembourg", "NL": "Netherlands",
      "NO": "Norway", "PL": "Poland", "PT": "Portugal", "SK": "Slovakia", "SI": "Slovenia",
      "ES": "Spain", "SE": "Sweden", "CH": "Switzerland", "UK": "United Kingdom"}
INV = {v: k for k, v in NM.items()}


def summarise(d):
    """Country-level summary, identical convention to stage_d.py."""
    return d.groupby("geo").apply(lambda s: pd.Series({
        "deaths": s.deaths.sum(), "dall": s.dall.sum(), "nyr": s.year.nunique(),
        "asmr": np.average(s.asmr, weights=s.expo),
        "asyr_nat": np.average(s.asyr_nat, weights=s.expo),
        "asyr_gbd": np.average(s.asyr_gbd, weights=s.expo),
        "yll_nat": s.yll_nat.sum(), "yll_gbd": s.yll_gbd.sum(),
        "gain_e0": s.gain_e0.mean(), "d_edag": s.d_edag.mean(),
        "edag": s.edag.mean(), "e0": s.e0.mean(), "thresh": s.thresh.median(),
        "mad": np.average(s.mad, weights=s.deaths),
        "below": np.average(s.frac_below_thresh, weights=s.deaths)}),
        include_groups=False)


T = R[R.sex == "T"]
prim = summarise(T[T.year.between(2012, 2018)])

G = {}

# ---------------------------------------------------------------- balanced panel
drop = [INV[c] for c in ["Iceland", "Luxembourg"]]
bal = prim.drop(index=drop)
G["bal_n"] = int(len(bal))
G["bal_excluded"] = ["Iceland", "Luxembourg"]
G["bal_gain_e0"] = round(float(bal.gain_e0.mean()), 3)
G["bal_d_edag"] = round(float(bal.d_edag.mean()), 4)
G["bal_d_edag_positive"] = int((bal.d_edag > 0).sum())
G["bal_gain_lo"] = round(float(bal.gain_e0.min()), 3)
G["bal_gain_hi"] = round(float(bal.gain_e0.max()), 3)
G["bal_d_edag_lo"] = round(float(bal.d_edag.min()), 4)
G["bal_d_edag_hi"] = round(float(bal.d_edag.max()), 4)
G["bal_thresh_median"] = float(bal.thresh.median())
# difference from the full 28-country frame
G["full_gain_e0"] = round(float(prim.gain_e0.mean()), 3)
G["full_d_edag"] = round(float(prim.d_edag.mean()), 4)
G["bal_gain_diff"] = round(abs(bal.gain_e0.mean() - prim.gain_e0.mean()), 4)
G["bal_d_edag_diff"] = round(abs(bal.d_edag.mean() - prim.d_edag.mean()), 5)
# ordering: rank of the 26 retained countries, full frame versus balanced frame
r_full = prim.loc[bal.index].asyr_nat.rank(ascending=False)
r_bal = bal.asyr_nat.rank(ascending=False)
G["bal_rho"] = round(float(stats.spearmanr(r_full, r_bal).statistic), 4)
G["bal_max_shift"] = int((r_full - r_bal).abs().max())
G["bal_identical_order"] = bool((r_full == r_bal).all())

# ------------------------------------------------------------------ later period
lat = [INV[c] for c in P["later_period_countries"]]
L = summarise(T[(T.geo.isin(lat)) & (T.year.between(2012, 2021))])
Lp = prim.loc[lat]
G["lat_n"] = int(len(L))
G["lat_years"] = P["later_period_years"]
G["lat_country_years"] = int(len(T[(T.geo.isin(lat)) & (T.year.between(2012, 2021))]))
G["lat_deaths"] = int(L.deaths.sum())
G["lat_gain_e0"] = round(float(L.gain_e0.mean()), 3)
G["lat_d_edag"] = round(float(L.d_edag.mean()), 4)
G["lat_d_edag_positive"] = int((L.d_edag > 0).sum())
G["lat_d_edag_lo"] = round(float(L.d_edag.min()), 4)
G["lat_d_edag_hi"] = round(float(L.d_edag.max()), 4)
G["lat_thresh_median"] = float(L.thresh.median())
G["lat_mad_lo"] = round(float(L.mad.min()), 1)
G["lat_mad_hi"] = round(float(L.mad.max()), 1)
G["lat_below_lo"] = int(round(100 * L.below.min()))
G["lat_below_hi"] = int(round(100 * L.below.max()))
# same countries, primary window, for a like-for-like comparison
G["latprim_gain_e0"] = round(float(Lp.gain_e0.mean()), 3)
G["latprim_d_edag"] = round(float(Lp.d_edag.mean()), 4)
G["lat_gain_diff"] = round(abs(L.gain_e0.mean() - Lp.gain_e0.mean()), 4)
G["lat_d_edag_diff"] = round(abs(L.d_edag.mean() - Lp.d_edag.mean()), 5)
r_lp = Lp.asyr_nat.rank(ascending=False)
r_l = L.asyr_nat.rank(ascending=False)
G["lat_rho"] = round(float(stats.spearmanr(r_lp, r_l).statistic), 4)
G["lat_max_shift"] = int((r_lp - r_l).abs().max())
# reference dependence within the later window
r_lnat = L.asyr_nat.rank(ascending=False)
r_lgbd = L.asyr_gbd.rank(ascending=False)
G["lat_ref_shift_ge2"] = int(((r_lnat - r_lgbd).abs() >= 2).sum())
G["lat_ref_max_shift"] = int((r_lnat - r_lgbd).abs().max())
# rank agreement reported without a displacement cut, as in stage_q
G["lat_ref_rho"] = round(float(stats.spearmanr(L.asyr_nat, L.asyr_gbd).statistic), 3)
G["lat_ref_tau"] = round(float(stats.kendalltau(L.asyr_nat, L.asyr_gbd).statistic), 3)
G["lat_ref_median_shift"] = float((r_lnat - r_lgbd).abs().median())
_n = len(L)
_a, _b = L.asyr_nat.values, L.asyr_gbd.values
_disc = sum(1 for i in range(_n) for j in range(i + 1, _n)
            if (_a[i] - _a[j]) * (_b[i] - _b[j]) < 0)
G["lat_ref_pairs"] = _n * (_n - 1) // 2
G["lat_ref_pairs_reversed"] = _disc
G["lat_ref_pct_pairs_reversed"] = round(100.0 * _disc / G["lat_ref_pairs"], 1)
G["lat_ypd_ratio_mean"] = round(float((L.yll_gbd / L.yll_nat).mean()), 3)
G["lat_ref_pct_mean"] = round(100 * (float((L.yll_gbd / L.yll_nat).mean()) - 1), 1)

json.dump(G, open("SENS3.json", "w"), indent=1)
for k, v in G.items():
    print(f"  {k:24s} {v}")
