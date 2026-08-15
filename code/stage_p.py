"""Stage P: sensitivity of the results to the closure of the open age interval.

The analysis collapses ages above 95 into a single open interval and closes it with
L = l / m, so remaining expectancy at exact age 95 is 1 / m(95+). The published Human
Mortality Database series instead extends to 110 years and closes there. The two give
different remaining expectancy at 95, and because T(x) is cumulated from the oldest ages
downward, the difference propagates to every younger age. Life-table levels therefore differ
slightly between the reconstruction and the published series.

Every quantity reported in this paper is either a difference between two life tables built by
the same procedure, or a ratio in which the closure enters both numerator and denominator, so
the choice is expected to cancel. This stage measures by how much it does.

The alternative closure sets the open interval to the published remaining expectancy at exact
age 95. Two treatments of the cause-deleted table bracket the reasonable range:

  fixed   the deleted table takes the same open-interval expectancy as the observed table, so
          the open interval does not respond to cause deletion at all;
  scaled  the deleted table's open-interval expectancy is divided by R(95+), reproducing the
          proportional response implied by the L = l / m closure at the published level.

Writes CLOSURE.json and closure_countryyear.csv.
"""
import pandas as pd, numpy as np, json, os
from engine import KEYS, N, ESP, build_lt, e_within

MAPc = {"AT": "AUT", "BE": "BEL", "BG": "BGR", "HR": "HRV", "CZ": "CZE", "DK": "DNK",
        "EE": "EST", "FI": "FIN", "FR": "FRATNP", "DE": "DEUTNP", "EL": "GRC",
        "HU": "HUN", "IS": "ISL", "IE": "IRL", "IT": "ITA", "LV": "LVA", "LT": "LTU",
        "LU": "LUX", "NL": "NLD", "NO": "NOR", "PL": "POL", "PT": "PRT", "SK": "SVK",
        "SI": "SVN", "ES": "ESP", "SE": "SWE", "CH": "CHE", "UK": "GBR_NP"}


def read_hmd(p):
    L = open(p).read().split("\n")
    h = [i for i, l in enumerate(L) if l.strip().startswith("Year")][0]
    cols = L[h].split()
    df = pd.DataFrame([l.split() for l in L[h + 1:] if l.strip()], columns=cols)
    df["Year"] = df.Year.astype(int)
    for c in cols[2:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def build_lt_e95(mx, ax, e95):
    """engine.build_lt with the open interval closed at a supplied remaining expectancy."""
    n = N.copy(); k = len(mx)
    qx = np.empty(k)
    qx[:-1] = n[:-1] * mx[:-1] / (1 + (n[:-1] - ax[:-1]) * mx[:-1])
    qx = np.clip(qx, 0, 1); qx[-1] = 1.0
    lx = np.empty(k); lx[0] = 1.0
    for i in range(k - 1):
        lx[i + 1] = lx[i] * (1 - qx[i])
    dx = np.empty(k); dx[:-1] = lx[:-1] - lx[1:]; dx[-1] = lx[-1]
    Lx = np.empty(k)
    Lx[:-1] = n[:-1] * lx[1:] + ax[:-1] * dx[:-1]
    Lx[-1] = lx[-1] * e95
    Tx = np.flip(np.cumsum(np.flip(Lx))); ex = np.where(lx > 0, Tx / lx, 0.0)
    return dict(lx=lx, dx=dx, Tx=Tx, ex=ex)


def edag(lt, ax):
    return float((lt["dx"] * e_within(lt["ex"], ax)).sum() / lt["lx"][0])


H = pd.read_csv("hmd_all.csv"); H = H[H.sex == "T"]
A = pd.read_csv("agedist_primary.csv")
CY = pd.read_csv("results_countryyear.csv")
frame = CY[(CY.sex == "T") & (CY.year.between(2012, 2018))][["geo", "year"]].rename(
    columns={"year": "Year"})

rows, tot = [], dict(dcns=0.0, dall=0.0, dall_open=0.0,
                     ycns_rep=0.0, yall_rep=0.0, ycns_alt=0.0, yall_alt=0.0)
for (geo, yr), sub in H.merge(frame, on=["geo", "Year"]).groupby(["geo", "Year"]):
    pub = read_hmd(f"raw/{MAPc[geo]}_bltper_5x1.txt")
    pub = pub[pub.Year == yr]
    if pub.empty:
        raise SystemExit(f"missing published life table for {geo} {yr}")
    s = pub.set_index("Age")
    e95_pub = float(s.Tx["95-99"]) / float(s.lx["95-99"])
    e0_pub = float(s.ex["0"])

    sub = sub.sort_values("ord")
    mx, ax, expo = sub.mx.values, sub.ax.values, sub.expo.values
    dall = sub.D_all_hmd.values
    a = A[(A.geo == geo) & (A.year == yr)].set_index("k").reindex(KEYS)
    dcns = np.nan_to_num(a.d.values)
    Rx = np.where(dall > 0, 1 - np.divide(dcns, dall, out=np.zeros_like(dcns),
                                          where=dall > 0), 1.0)

    lt = build_lt(mx, ax); ltd = build_lt(mx * Rx, ax)
    lt2 = build_lt_e95(mx, ax, e95_pub)
    ltd_f = build_lt_e95(mx * Rx, ax, e95_pub)
    ltd_s = build_lt_e95(mx * Rx, ax, e95_pub / Rx[-1] if Rx[-1] > 0 else e95_pub)

    en_rep = e_within(lt["ex"], ax); en_alt = e_within(lt2["ex"], ax)
    tot["dcns"] += dcns.sum(); tot["dall"] += dall.sum(); tot["dall_open"] += dall[-1]
    tot["ycns_rep"] += (dcns * en_rep).sum(); tot["yall_rep"] += (dall * en_rep).sum()
    tot["ycns_alt"] += (dcns * en_alt).sum(); tot["yall_alt"] += (dall * en_alt).sum()

    r = np.divide(dcns * en_rep, expo, out=np.zeros_like(dcns), where=expo > 0)
    r2 = np.divide(dcns * en_alt, expo, out=np.zeros_like(dcns), where=expo > 0)

    rows.append(dict(geo=geo, year=yr, expo=expo.sum(),
                     e0_published=e0_pub, e0_reported=float(lt["ex"][0]),
                     e0_alt=float(lt2["ex"][0]),
                     e95_published=e95_pub, e95_reported=float(lt["ex"][-1]),
                     edag_reported=edag(lt, ax), edag_alt=edag(lt2, ax),
                     gain_reported=float(ltd["ex"][0] - lt["ex"][0]),
                     gain_alt_fixed=float(ltd_f["ex"][0] - lt2["ex"][0]),
                     gain_alt_scaled=float(ltd_s["ex"][0] - lt2["ex"][0]),
                     ded_reported=edag(lt, ax) - edag(ltd, ax),
                     ded_alt_fixed=edag(lt2, ax) - edag(ltd_f, ax),
                     ded_alt_scaled=edag(lt2, ax) - edag(ltd_s, ax),
                     asyr_reported=float((r * ESP).sum() * 1e5),
                     asyr_alt=float((r2 * ESP).sum() * 1e5),
                     cns_open=float(dcns[-1]), cns_total=float(dcns.sum())))

D = pd.DataFrame(rows)
D.to_csv("closure_countryyear.csv", index=False)
C = D.groupby("geo").apply(lambda x: pd.Series({
    "e0_published": x.e0_published.mean(), "e0_reported": x.e0_reported.mean(),
    "edag_reported": x.edag_reported.mean(), "edag_alt": x.edag_alt.mean(),
    "gain_reported": x.gain_reported.mean(),
    "gain_alt_fixed": x.gain_alt_fixed.mean(), "gain_alt_scaled": x.gain_alt_scaled.mean(),
    "ded_reported": x.ded_reported.mean(),
    "ded_alt_fixed": x.ded_alt_fixed.mean(), "ded_alt_scaled": x.ded_alt_scaled.mean(),
    "asyr_reported": np.average(x.asyr_reported, weights=x.expo),
    "asyr_alt": np.average(x.asyr_alt, weights=x.expo)}), include_groups=False)

sd = 100 * tot["dcns"] / tot["dall"]
G = {
 "n_country_years": int(len(D)),
 "e0_recon_minus_published_mean": round(float((D.e0_reported - D.e0_published).mean()), 4),
 "e0_recon_minus_published_max": round(float((D.e0_reported - D.e0_published).abs().max()), 4),
 "e95_published_mean": round(float(D.e95_published.mean()), 2),
 "e95_reported_mean": round(float(D.e95_reported.mean()), 2),
 "edag_reported_mean": round(float(C.edag_reported.mean()), 3),
 "edag_alt_mean": round(float(C.edag_alt.mean()), 3),
 "gain_reported": round(float(C.gain_reported.mean()), 4),
 "gain_alt_fixed": round(float(C.gain_alt_fixed.mean()), 4),
 "gain_alt_scaled": round(float(C.gain_alt_scaled.mean()), 4),
 "gain_max_country_diff": round(float(pd.concat([
     (C.gain_alt_fixed - C.gain_reported).abs(),
     (C.gain_alt_scaled - C.gain_reported).abs()]).max()), 5),
 "ded_reported": round(float(C.ded_reported.mean()), 4),
 "ded_alt_fixed": round(float(C.ded_alt_fixed.mean()), 4),
 "ded_alt_scaled": round(float(C.ded_alt_scaled.mean()), 4),
 "ded_max_country_diff": round(float(pd.concat([
     (C.ded_alt_fixed - C.ded_reported).abs(),
     (C.ded_alt_scaled - C.ded_reported).abs()]).max()), 5),
 "ded_positive_reported": int((C.ded_reported > 0).sum()),
 "ded_positive_alt_fixed": int((C.ded_alt_fixed > 0).sum()),
 "ded_positive_alt_scaled": int((C.ded_alt_scaled > 0).sum()),
 "share_deaths_pct": round(sd, 3),
 "share_yll_reported_pct": round(100 * tot["ycns_rep"] / tot["yall_rep"], 3),
 "share_yll_alt_pct": round(100 * tot["ycns_alt"] / tot["yall_alt"], 3),
 "asyr_max_country_pct_diff": round(float(
     (100 * (C.asyr_alt - C.asyr_reported).abs() / C.asyr_reported).max()), 2),
 "asyr_rank_changes": int((C.asyr_reported.rank(ascending=False)
                           != C.asyr_alt.rank(ascending=False)).sum()),
 "asyr_rank_max_shift": int((C.asyr_reported.rank(ascending=False)
                             - C.asyr_alt.rank(ascending=False)).abs().max()),
 "allcause_deaths_open_pct": round(100 * tot["dall_open"] / tot["dall"], 1),
 "cns_deaths_open_pct": round(100 * D.cns_open.sum() / D.cns_total.sum(), 2),
}
json.dump(G, open("CLOSURE.json", "w"), indent=1)
for k, v in G.items():
    print("  %-32s %s" % (k, v))
