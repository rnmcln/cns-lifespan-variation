"""Stage S. Size and stability of the rank reversals between reference life tables.

Stage Q counts the pairs of countries whose ordering by the age-standardised YLL rate reverses
when the national life tables are replaced by a common reference. A count alone does not say
whether the reversals involve countries with clearly different rates or countries whose rates
are almost equal. This stage reports, for malignant CNS neoplasms:

1. the size of the gap between the two countries of each reversed pair, under each reference,
   in absolute terms (per 100 000) and relative to the pair mean, compared with all pairs;
2. the proportion of Poisson replicates of the recorded deaths in which each reversal recurs,
   and the distribution of the number of reversed pairs across replicates;
3. the same rank statistics for YLL per death, as a contrast with the standardised rate.

It then repeats the rank comparison for the nine causes of stage J on the common country frame,
to show where the reference effect on country orderings is largest and smallest.

Inputs : hmd_all.csv, country_summary_all.csv, eurostat_full.csv, results_countryyear.csv, official_reference_country.csv,
         gbd/tmrlt_gbd2023.csv, raw/*_Deaths_5x1.txt, cause_age_deaths.csv,
         cause_age_weights.csv (undistributed, computed from Human Mortality Database series),
         generalisation_country.csv
Outputs: REVERSALS.json, reversal_pairs.csv, rank_by_cause.csv, asyr_change_country.csv
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

from engine import ESP, KEYS, build_lt, e_within  # noqa: F401

SEED = 20260803
NREP = 2000
THR = (5.0, 10.0)   # relative-gap thresholds, % of the pair mean

O = pd.read_csv("gbd/tmrlt_gbd2023.csv")
oa = O["Age"].values.astype(float)
oe = O["Life Expectancy"].values
e_star = lambda a: np.interp(a, oa, oe)          # noqa: E731
STARTS = np.array([0, 1] + list(range(5, 95, 5)), float)
WID = np.array([1.0, 4.0] + [5.0] * 18)
SUB_E = e_star(np.array([97.5, 102.5, 107.5, 110.0]))
OPEN = ["95-99", "100-104", "105-109", "110+"]
MAPc = {"AT": "AUT", "BE": "BEL", "BG": "BGR", "HR": "HRV", "CZ": "CZE", "DK": "DNK",
        "EE": "EST", "FI": "FIN", "FR": "FRATNP", "DE": "DEUTNP", "EL": "GRC", "HU": "HUN",
        "IS": "ISL", "IE": "IRL", "IT": "ITA", "LV": "LVA", "LT": "LTU", "LU": "LUX",
        "NL": "NLD", "NO": "NOR", "PL": "POL", "PT": "PRT", "SK": "SVK", "SI": "SVN",
        "ES": "ESP", "SE": "SWE", "CH": "CHE", "UK": "GBR_NP"}


def read_hmd(path):
    L = open(path).read().split("\n")
    h = [i for i, l in enumerate(L) if l.strip().startswith("Year")][0]
    cols = L[h].split()
    df = pd.DataFrame([l.split() for l in L[h + 1:] if l.strip()], columns=cols)
    df["Year"] = df.Year.astype(int)
    for c in cols[2:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


OW = {}
for geo, cc in MAPc.items():
    d = read_hmd("raw/%s_Deaths_5x1.txt" % cc)
    d = d[d.Age.isin(OPEN)]
    for yr, s in d.groupby("Year"):
        w = s.set_index("Age")["Total"].reindex(OPEN).fillna(0).values.astype(float)
        OW[(geo, yr)] = w / w.sum() if w.sum() > 0 else np.array([1.0, 0, 0, 0])


def pair_table(a, b, names):
    """All pairs of countries with their gaps under two rate vectors a and b."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    rows = []
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            ga, gb = a[i] - a[j], b[i] - b[j]
            rows.append(dict(
                geo_1=names[i], geo_2=names[j],
                rate_a_1=a[i], rate_a_2=a[j], rate_b_1=b[i], rate_b_2=b[j],
                gap_a_abs=abs(ga), gap_b_abs=abs(gb),
                gap_a_rel=100 * abs(ga) / ((a[i] + a[j]) / 2),
                gap_b_rel=100 * abs(gb) / ((b[i] + b[j]) / 2),
                reversed=bool(ga * gb < 0)))
    P = pd.DataFrame(rows)
    P["gap_min_rel"] = P[["gap_a_rel", "gap_b_rel"]].min(axis=1)
    return P


def agreement(a, b):
    ra = pd.Series(np.asarray(a, float)).rank(ascending=False)
    rb = pd.Series(np.asarray(b, float)).rank(ascending=False)
    P = pair_table(a, b, list(range(len(a))))
    n_rev = int(P.reversed.sum())
    rev = P[P.reversed]
    out = {
        "n": int(len(ra)), "n_pairs": int(len(P)),
        "spearman_rho": round(float(stats.spearmanr(a, b).statistic), 3),
        "kendall_tau_b": round(float(stats.kendalltau(a, b).statistic), 3),
        "max_abs_shift": int((ra - rb).abs().max()),
        "n_shift_ge2": int((ra - rb).abs().ge(2).sum()),
        "n_pairs_reversed": n_rev,
        "pct_pairs_reversed": round(100.0 * n_rev / len(P), 1),
        "all_gap_a_rel_median": round(float(P.gap_a_rel.median()), 1),
    }
    if n_rev:
        out["rev_gap_a_rel_median"] = round(float(rev.gap_a_rel.median()), 1)
        out["rev_gap_b_rel_median"] = round(float(rev.gap_b_rel.median()), 1)
        out["rev_gap_a_rel_max"] = round(float(rev.gap_a_rel.max()), 1)
        out["rev_gap_b_rel_max"] = round(float(rev.gap_b_rel.max()), 1)
        for t in THR:
            out["n_rev_both_ge%d" % t] = int((rev.gap_min_rel >= t).sum())
    else:
        out.update({k: None for k in ["rev_gap_a_rel_median", "rev_gap_b_rel_median",
                                      "rev_gap_a_rel_max", "rev_gap_b_rel_max"]})
        for t in THR:
            out["n_rev_both_ge%d" % t] = 0
    return out, P


def poisson_reversals(cells, rng):
    """cells: dict geo -> (d, w_a, w_b, expo_y) with d of shape (years, 21).
    w_* are per-cell weights such that the annual rate is (d * w).sum(axis=1);
    expo_y are the person-years used to average the annual rates.
    Returns replicate country rates under both references, shape (NREP, n_countries)."""
    geos = sorted(cells)
    A = np.empty((NREP, len(geos)))
    B = np.empty((NREP, len(geos)))
    for k, g in enumerate(geos):
        d, wa, wb, ey = cells[g]
        D = rng.poisson(np.broadcast_to(d, (NREP,) + d.shape))
        ra = (D * wa).sum(axis=2)
        rb = (D * wb).sum(axis=2)
        A[:, k] = (ra * ey).sum(axis=1) / ey.sum()
        B[:, k] = (rb * ey).sum(axis=1) / ey.sum()
    return geos, A, B


def replicate_summary(geos, A, B, P):
    """Recurrence of each observed reversal, and replicate counts of reversed pairs."""
    idx = {g: i for i, g in enumerate(geos)}
    i1 = P.geo_1.map(idx).values
    i2 = P.geo_2.map(idx).values
    sa = np.sign(P.rate_a_1.values - P.rate_a_2.values)
    sb = np.sign(P.rate_b_1.values - P.rate_b_2.values)
    da = np.sign(A[:, i1] - A[:, i2])
    db = np.sign(B[:, i1] - B[:, i2])
    same_a = (da == sa).mean(axis=0)
    same_b = (db == sb).mean(axis=0)
    both = ((da == sa) & (db == sb)).mean(axis=0)
    nrev = ((da * db) < 0).sum(axis=1)
    return same_a, same_b, both, nrev


rng = np.random.default_rng(SEED)
J = {"seed": SEED, "n_replicates": NREP, "rel_gap_thresholds_pct": list(THR)}

# ------------------------------------------------------------- 1. CNS, primary frame
H = pd.read_csv("hmd_all.csv")
E = pd.read_csv("eurostat_full.csv")
Rr = pd.read_csv("results_countryyear.csv")
g = pd.read_csv("official_reference_country.csv").set_index("geo").sort_index()

cns = E[E.which == "cns"].rename(columns={"time": "Year"})
cns = cns[cns.k.notna()]
piv = cns.pivot_table(index=["geo", "Year", "sex", "k"], values="deaths", dropna=False).reset_index()
ok = (Rr[(Rr.sex == "T") & (Rr.year.between(2012, 2018))][["geo", "sex", "year"]]
      .rename(columns={"year": "Year"}))
M = (H.merge(piv, on=["geo", "Year", "sex", "k"], how="inner")
      .merge(ok, on=["geo", "sex", "Year"]).sort_values(["geo", "Year", "ord"]))

cells = {}
chk = []
for geo, z in M.groupby("geo"):
    ds, wn, wo, ey = [], [], [], []
    for yr, s in z.groupby("Year"):
        s = s.sort_values("ord")
        d = np.nan_to_num(s.deaths.values)
        ax = s.ax.values
        expo = s.expo.values
        e_nat = e_within(build_lt(s.mx.values, ax)["ex"], ax)
        eo = np.empty(21)
        eo[:20] = e_star(STARTS + np.clip(ax[:20] / WID, 0, 1) * WID)
        eo[20] = float((OW[(geo, yr)] * SUB_E).sum())
        w = np.divide(ESP * 1e5, expo, out=np.zeros_like(expo, dtype=float), where=expo > 0)
        ds.append(d)
        wn.append(e_nat * w)
        wo.append(eo * w)
        ey.append(float(expo.sum()))
    d, wn, wo, ey = map(np.array, (ds, wn, wo, ey))
    cells[geo] = (d, wn, wo, ey)
    an = float(((d * wn).sum(axis=1) * ey).sum() / ey.sum())
    ao = float(((d * wo).sum(axis=1) * ey).sum() / ey.sum())
    chk.append((geo, an, ao))
C = pd.DataFrame(chk, columns=["geo", "asyr_nat", "asyr_off"]).set_index("geo").sort_index()
# the recomputed point estimates must equal the published ones
J["check_max_absdiff_nat"] = float(np.abs(C.asyr_nat - g.asyr_nat).max())
J["check_max_absdiff_off"] = float(np.abs(C.asyr_off - g.asyr_off).max())
assert J["check_max_absdiff_nat"] < 1e-6 and J["check_max_absdiff_off"] < 1e-6

names = list(g.index)

# ---- change in each country's age-standardised YLL rate between references
# The comparison of countries is made on this rate, so its change carries the main result; YLL per
# death is reported separately and labelled as such.
S_ = pd.read_csv("country_summary_all.csv")
e0 = S_[(S_.frame == "primary_2012_2018") & (S_.sex == "T")].set_index("geo").e0.reindex(g.index)
CH = pd.DataFrame({"asyr_nat": g.asyr_nat, "asyr_gbd": g.asyr_off, "asyr_minmort": g.asyr_bp,
                   "e0_national": e0})
CH["change_gbd_abs"] = CH.asyr_gbd - CH.asyr_nat
CH["change_gbd_pct"] = 100 * (CH.asyr_gbd / CH.asyr_nat - 1)
CH["change_minmort_pct"] = 100 * (CH.asyr_minmort / CH.asyr_nat - 1)
CH.round(4).to_csv("asyr_change_country.csv")
for tag, col in (("gbd", "change_gbd_pct"), ("minmort", "change_minmort_pct")):
    c = CH[col]
    J["asyr_change_" + tag] = {
        "median_pct": round(float(c.median()), 1), "mean_pct": round(float(c.mean()), 1),
        "lo_pct": round(float(c.min()), 1), "hi_pct": round(float(c.max()), 1),
        "lo_c": c.idxmin(), "hi_c": c.idxmax(),
        "q1_pct": round(float(c.quantile(0.25)), 1), "q3_pct": round(float(c.quantile(0.75)), 1),
        "spearman_with_e0": round(float(stats.spearmanr(c, CH.e0_national).statistic), 3)}
J["asyr_change_gbd"]["median_abs"] = round(float(CH.change_gbd_abs.median()), 1)
J["asyr_change_gbd"]["lo_abs"] = round(float(CH.change_gbd_abs.min()), 1)
J["asyr_change_gbd"]["hi_abs"] = round(float(CH.change_gbd_abs.max()), 1)
CEE_ = {"BG", "HR", "CZ", "EE", "HU", "LV", "LT", "PL", "SK", "SI"}
top = CH.change_gbd_pct.sort_values(ascending=False)
J["asyr_change_gbd"]["top5"] = list(top.index[:5])
J["asyr_change_gbd"]["top10_cee"] = int(sum(x in CEE_ for x in top.index[:10]))
J["asyr_change_gbd"]["cee_median_pct"] = round(float(top[top.index.isin(CEE_)].median()), 1)
J["asyr_change_gbd"]["other_median_pct"] = round(float(top[~top.index.isin(CEE_)].median()), 1)
res = {}
for lab, a, b in [("nat_vs_gbd", g.asyr_nat, g.asyr_off),
                  ("nat_vs_minmort", g.asyr_nat, g.asyr_bp),
                  ("ypd_nat_vs_gbd", g.ypd_nat, g.ypd_off)]:
    r, P = agreement(a.values, b.values)
    P = pair_table(a.values, b.values, names)
    res[lab] = (r, P)
    J[lab] = r

# gap distributions for the primary comparison
r, P = res["nat_vs_gbd"]
rev = P[P.reversed]
con = P[~P.reversed]
J["nat_vs_gbd"].update({
    "rev_gap_a_abs_median": round(float(rev.gap_a_abs.median()), 1),
    "rev_gap_a_abs_q1": round(float(rev.gap_a_abs.quantile(0.25)), 1),
    "rev_gap_a_abs_q3": round(float(rev.gap_a_abs.quantile(0.75)), 1),
    "rev_gap_a_abs_max": round(float(rev.gap_a_abs.max()), 1),
    "rev_gap_b_abs_median": round(float(rev.gap_b_abs.median()), 1),
    "rev_gap_b_abs_max": round(float(rev.gap_b_abs.max()), 1),
    "rev_gap_a_rel_q1": round(float(rev.gap_a_rel.quantile(0.25)), 1),
    "rev_gap_a_rel_q3": round(float(rev.gap_a_rel.quantile(0.75)), 1),
    "rev_gap_b_rel_q1": round(float(rev.gap_b_rel.quantile(0.25)), 1),
    "rev_gap_b_rel_q3": round(float(rev.gap_b_rel.quantile(0.75)), 1),
    "con_gap_a_rel_median": round(float(con.gap_a_rel.median()), 1),
    "asyr_nat_range_rel_pct": round(float(100 * (g.asyr_nat.max() / g.asyr_nat.min() - 1)), 1),
})

# Poisson recurrence of the reversals, national vs GBD 2023
geos, A, B = poisson_reversals(cells, rng)
assert geos == names
same_a, same_b, both, nrev = replicate_summary(geos, A, B, P)
P["p_order_nat"] = same_a
P["p_order_gbd"] = same_b
P["p_reversal_recurs"] = both
rev = P[P.reversed]
J["nat_vs_gbd"].update({
    "rep_n_reversed_median": int(np.median(nrev)),
    "rep_n_reversed_lo": int(np.percentile(nrev, 2.5)),
    "rep_n_reversed_hi": int(np.percentile(nrev, 97.5)),
    "rev_recur_ge95": int((rev.p_reversal_recurs >= 0.95).sum()),
    "rev_recur_ge80": int((rev.p_reversal_recurs >= 0.80).sum()),
    "rev_recur_median": round(float(rev.p_reversal_recurs.median()), 3),
    "all_nat_order_ge95": int((P.p_order_nat >= 0.95).sum()),
    "rev_nat_order_ge95": int((rev.p_order_nat >= 0.95).sum()),
})

# countries involved in reversals that recur in at least 95% of replicates
st = rev[rev.p_reversal_recurs >= 0.95]
inv = pd.concat([st.geo_1, st.geo_2]).value_counts()
J["nat_vs_gbd"]["stable_rev_countries"] = {k: int(v) for k, v in inv.sort_index().items()}
# the largest stable reversals, ordered by the smaller of the two relative gaps
top = st.sort_values("gap_min_rel", ascending=False).head(5)
J["nat_vs_gbd"]["largest_stable_reversals"] = [
    dict(higher_nat=(a1 if r1 > r2 else a2), lower_nat=(a2 if r1 > r2 else a1),
         gap_nat_rel=round(gn, 1), gap_gbd_rel=round(go, 1),
         gap_nat_abs=round(gna, 1), gap_gbd_abs=round(goa, 1), p_recur=round(p, 3))
    for a1, a2, r1, r2, gn, go, gna, goa, p in zip(
        top.geo_1, top.geo_2, top.rate_a_1, top.rate_a_2, top.gap_a_rel, top.gap_b_rel,
        top.gap_a_abs, top.gap_b_abs, top.p_reversal_recurs)]

CEE = {"BG", "HR", "CZ", "EE", "HU", "LV", "LT", "PL", "SK", "SI"}
gain = np.where(st.rate_b_1 > st.rate_b_2, st.geo_1, st.geo_2)   # higher under GBD 2023
lose = np.where(st.rate_b_1 > st.rate_b_2, st.geo_2, st.geo_1)
J["nat_vs_gbd"]["stable_rev_cee_moves_above_other"] = int(sum(
    (a in CEE) and (b not in CEE) for a, b in zip(gain, lose)))
J["nat_vs_gbd"]["stable_rev_other_moves_above_cee"] = int(sum(
    (a not in CEE) and (b in CEE) for a, b in zip(gain, lose)))
J["nat_vs_gbd"]["stable_rev_within_group"] = int(len(st)) - \
    J["nat_vs_gbd"]["stable_rev_cee_moves_above_other"] - J["nat_vs_gbd"]["stable_rev_other_moves_above_cee"]
J["nat_vs_gbd"]["stable_rev_gap_nat_rel_median"] = round(float(st.gap_a_rel.median()), 1)
J["nat_vs_gbd"]["stable_rev_gap_gbd_rel_median"] = round(float(st.gap_b_rel.median()), 1)
J["nat_vs_gbd"]["stable_rev_gap_nat_abs_median"] = round(float(st.gap_a_abs.median()), 1)
J["nat_vs_gbd"]["stable_rev_gap_gbd_abs_median"] = round(float(st.gap_b_abs.median()), 1)
J["nat_vs_gbd"]["n_pairs_nat_gap_gt_rev_max"] = int((P.gap_a_rel > rev.gap_a_rel.max()).sum())

out = P.rename(columns={"rate_a_1": "asyr_nat_1", "rate_a_2": "asyr_nat_2",
                        "rate_b_1": "asyr_gbd_1", "rate_b_2": "asyr_gbd_2",
                        "gap_a_abs": "gap_nat_abs", "gap_b_abs": "gap_gbd_abs",
                        "gap_a_rel": "gap_nat_rel_pct", "gap_b_rel": "gap_gbd_rel_pct",
                        "gap_min_rel": "gap_min_rel_pct"})
out.round(4).to_csv("reversal_pairs.csv", index=False)

# --------------------------------------------- 2. nine causes, common country frame
GC = pd.read_csv("generalisation_country.csv")
AG = pd.read_csv("cause_age_deaths.csv").merge(pd.read_csv("cause_age_weights.csv"),
                                                on=["icd", "geo", "year", "age"], validate="1:1")
common = sorted(set.intersection(*[set(z.geo) for _, z in GC.groupby("icd")]))
order = (pd.read_csv("generalisation_common.csv")[["icd", "cause", "mad"]])
rows = []
for icd, cause, mad in order.itertuples(index=False):
    z = GC[(GC.icd == icd) & GC.geo.isin(common)].set_index("geo").sort_index()
    r, Pc = agreement(z.asyr_nat.values, z.asyr_off.values)
    Pc = pair_table(z.asyr_nat.values, z.asyr_off.values, list(z.index))
    # Poisson recurrence for this cause, from its own age-specific record
    a = AG[(AG.icd == icd) & AG.geo.isin(common)]
    cc = {}
    for geo, s in a.groupby("geo"):
        yrs = sorted(s.year.unique())
        dd = np.array([s[s.year == y].deaths.values for y in yrs], float)
        ex = np.array([s[s.year == y].expo.values for y in yrs], float)
        w = np.divide(ESP * 1e5, ex, out=np.zeros_like(ex), where=ex > 0)
        wn = np.array([s[s.year == y].e_nat.values for y in yrs]) * w
        wo = np.array([s[s.year == y].e_off.values for y in yrs]) * w
        cc[geo] = (dd, wn, wo, ex.sum(axis=1))
    geos_c, Ac, Bc = poisson_reversals(cc, rng)
    assert geos_c == list(z.index)
    _, _, both_c, nrev_c = replicate_summary(geos_c, Ac, Bc, Pc)
    revc = Pc.reversed.values
    infl = 100 * (z.asyr_off / z.asyr_nat - 1)
    rows.append(dict(
        icd=icd, cause=cause, mad=round(float(mad), 1), n=r["n"],
        spearman_rho=r["spearman_rho"], kendall_tau_b=r["kendall_tau_b"],
        max_abs_shift=r["max_abs_shift"], n_shift_ge2=r["n_shift_ge2"],
        n_pairs_reversed=r["n_pairs_reversed"], pct_pairs_reversed=r["pct_pairs_reversed"],
        rev_recur_ge95=int((both_c[revc] >= 0.95).sum()),
        rev_gap_nat_rel_median=r["rev_gap_a_rel_median"],
        n_rev_both_ge10=r["n_rev_both_ge10"],
        all_gap_nat_rel_median=r["all_gap_a_rel_median"],
        infl_sd=round(float(infl.std(ddof=1)), 1),
        asyr_chg_median=round(float(infl.median()), 1),
        asyr_chg_lo=round(float(infl.min()), 1),
        asyr_chg_hi=round(float(infl.max()), 1),
        sd_log_rate_nat=round(float(np.log(z.asyr_nat).std(ddof=1)), 3),
        sd_log_ratio=round(float(np.log(z.asyr_off / z.asyr_nat).std(ddof=1)), 3),
        rep_n_reversed_median=int(np.median(nrev_c)),
        rep_n_reversed_lo=int(np.percentile(nrev_c, 2.5)),
        rep_n_reversed_hi=int(np.percentile(nrev_c, 97.5))))
RC = pd.DataFrame(rows)
RC["noise_ratio"] = (RC.sd_log_ratio / RC.sd_log_rate_nat).round(3)
RC.to_csv("rank_by_cause.csv", index=False)

J["by_cause_n_countries"] = len(common)
J["by_cause_excluded"] = sorted(set(GC.geo) - set(common))
J["by_cause_pct_rev_lo"] = float(RC.pct_pairs_reversed.min())
J["by_cause_pct_rev_hi"] = float(RC.pct_pairs_reversed.max())
J["by_cause_pct_rev_lo_cause"] = RC.loc[RC.pct_pairs_reversed.idxmin(), "cause"]
J["by_cause_pct_rev_hi_cause"] = RC.loc[RC.pct_pairs_reversed.idxmax(), "cause"]
J["by_cause_rho_lo"] = float(RC.spearman_rho.min())
J["by_cause_rho_hi"] = float(RC.spearman_rho.max())
J["by_cause_spearman_pctrev_noise_ratio"] = round(float(
    stats.spearmanr(RC.pct_pairs_reversed, RC.noise_ratio).statistic), 3)
J["by_cause_spearman_pctrev_mad"] = round(float(
    stats.spearmanr(RC.pct_pairs_reversed, RC.mad).statistic), 3)
J["by_cause_spearman_asyrchg_mad"] = round(float(
    stats.spearmanr(RC.asyr_chg_median, RC.mad).statistic), 3)
two = RC.nlargest(2, "noise_ratio")
rest = RC.drop(two.index)
J["by_cause_top2_noise_ratio_causes"] = list(two.cause)
J["by_cause_top2_are_top2_pctrev"] = bool(set(two.cause) == set(RC.nlargest(2, "pct_pairs_reversed").cause))
J["by_cause_rest_noise_ratio_lo"] = float(rest.noise_ratio.min())
J["by_cause_rest_noise_ratio_hi"] = float(rest.noise_ratio.max())
J["by_cause_rest_pct_rev_lo"] = float(rest.pct_pairs_reversed.min())
J["by_cause_rest_pct_rev_hi"] = float(rest.pct_pairs_reversed.max())
J["by_cause_rest_rho_lo"] = float(rest.spearman_rho.min())
J["by_cause_n_rev_both_ge10_total"] = int(RC.n_rev_both_ge10.sum())
cns_row = RC[RC.icd == "C70-C72"].iloc[0]
J["by_cause_cns_pct_rev"] = float(cns_row.pct_pairs_reversed)
J["by_cause_cns_rank_of_9"] = int(RC.pct_pairs_reversed.rank(ascending=False, method="min")[
    RC.icd == "C70-C72"].iloc[0])

json.dump(J, open("REVERSALS.json", "w"), indent=1)
print(RC.to_string(index=False))
print(json.dumps(J, indent=1))
