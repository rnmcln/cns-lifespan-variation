import pandas as pd, numpy as np, json
S=pd.read_csv("country_summary_all.csv"); R=pd.read_csv("results_countryyear.csv")
U=pd.read_csv("uncertainty_country.csv").set_index("geo")
A=pd.read_csv("allcause_yll.csv"); rec=pd.read_csv("reconciliation_allcause.csv")
sh=pd.read_csv("shares_primary.csv").set_index("geo"); X=pd.read_csv("exclusions.csv")
def f(x,n=3): return round(float(x),n)
P=S[(S.frame=="primary_2012_2018")&(S.sex=="T")].set_index("geo")
rp=R[(R.sex=="T")&(R.year.between(2012,2018))]
ap=A[(A.sex=="T")&(A.year.between(2012,2018))]
J={}
J["frame"]=dict(countries=int(len(P)),country_years=int(len(rp)),deaths=int(rp.deaths.sum()),
  allcause_deaths=int(rp.dall.sum()),period="2012-2018",
  excl_unreconciled=int(X[(X.sex=='T')&X.reason.str.startswith('unreconciled')].shape[0]),
  excl_nototal=int(X[(X.sex=='T')&(X.reason=='no_all_age_total')].shape[0]))
J["shares"]=dict(deaths_pct=f(100*rp.deaths.sum()/rp.dall.sum()),
  yll_pct=f(100*ap.yll_cns_nat.sum()/ap.yll_all_nat.sum()),
  ratio=f((ap.yll_cns_nat.sum()/ap.yll_all_nat.sum())/(rp.deaths.sum()/rp.dall.sum()),2),
  ratio_min=f(sh.ratio.min(),2),ratio_max=f(sh.ratio.max(),2),
  ratio_min_c=sh.ratio.idxmin(),ratio_max_c=sh.ratio.idxmax())
J["asmr"]=dict(lo=f(P.asmr.min(),2),hi=f(P.asmr.max(),2),lo_c=P.asmr.idxmin(),hi_c=P.asmr.idxmax())
J["asyr"]=dict(nat_lo=f(P.asyr_nat.min(),0),nat_hi=f(P.asyr_nat.max(),0),
  gbd_lo=f(P.asyr_gbd.min(),0),gbd_hi=f(P.asyr_gbd.max(),0),
  nat_lo_c=P.asyr_nat.idxmin(),nat_hi_c=P.asyr_nat.idxmax())
J["ypd"]=dict(nat_lo=f(P.ypd_nat.min(),1),nat_hi=f(P.ypd_nat.max(),1),
  nat_lo_c=P.ypd_nat.idxmin(),nat_hi_c=P.ypd_nat.idxmax(),
  ratio_mean=f(P.ratio_gbd_nat.mean()),ratio_lo=f(P.ratio_gbd_nat.min()),
  ratio_hi=f(P.ratio_gbd_nat.max()),ratio_lo_c=P.ratio_gbd_nat.idxmin(),
  ratio_hi_c=P.ratio_gbd_nat.idxmax(),
  pct_mean=f(100*(P.ratio_gbd_nat.mean()-1),1),pct_lo=f(100*(P.ratio_gbd_nat.min()-1),1),
  pct_hi=f(100*(P.ratio_gbd_nat.max()-1),1))
J["gain_e0"]=dict(mean=f(P.gain_e0.mean()),lo=f(P.gain_e0.min()),hi=f(P.gain_e0.max()),
  lo_c=P.gain_e0.idxmin(),hi_c=P.gain_e0.idxmax(),ci_excl_zero=int((U.glo>0).sum()))
J["d_edag"]=dict(mean=f(P.d_edag.mean(),4),lo=f(P.d_edag.min(),4),hi=f(P.d_edag.max(),4),
  lo_c=P.d_edag.idxmin(),hi_c=P.d_edag.idxmax(),positive=int((P.d_edag>0).sum()),
  ci_excl_zero=int((U.dlo>0).sum()),edag_mean=f(P.edag.mean(),2))
J["threshold"]=dict(median=f(P.thresh.median(),0),lo=f(rp.thresh.min(),0),hi=f(rp.thresh.max(),0),
  mad_lo=f(P.mad.min(),1),mad_hi=f(P.mad.max(),1),mad_lo_c=P.mad.idxmin(),mad_hi_c=P.mad.idxmax(),
  below_lo=f(100*P.below.min(),0),below_hi=f(100*P.below.max(),0))
piv=S[S.sex=="T"].pivot_table(index="geo",columns="frame",values="asyr_nat").rank(ascending=False)
J["frames"]=dict(rho_prim_max=f(piv.primary_2012_2018.corr(piv.maximal,method="spearman")),
  rho_prim_excl=f(piv.primary_2012_2018.corr(piv.excl_2020_21,method="spearman")),
  max_shift=int((piv.primary_2012_2018-piv.maximal).abs().max()))
rn=P.asyr_nat.rank(ascending=False); rg=P.asyr_gbd.rank(ascending=False)
J["ranks"]=dict(shift_ge2=int((rg-rn).abs().ge(2).sum()),max_shift=int((rg-rn).abs().max()))
r2=rec[rec.geo.isin(P.index)]
J["reconciliation"]=dict(median_abs=f(r2.pct.abs().median()),p95=f(r2.pct.abs().quantile(.95)),
  max=f(r2.pct.abs().max(),2),max_c=r2.loc[r2.pct.abs().idxmax(),"geo"],
  max_y=int(r2.loc[r2.pct.abs().idxmax(),"Year"]))
for sx,lab in [("M","male"),("F","female")]:
    Q=S[(S.frame=="primary_2012_2018")&(S.sex==sx)].set_index("geo")
    J[lab]=dict(asmr_lo=f(Q.asmr.min(),2),asmr_hi=f(Q.asmr.max(),2),
      ypd_lo=f(Q.ypd_nat.min(),1),ypd_hi=f(Q.ypd_nat.max(),1),
      gain=f(Q.gain_e0.mean()),ded=f(Q.d_edag.mean(),4),
      ded_pos=int((Q.d_edag>0).sum()),n=int(len(Q)),
      thresh=f(Q.thresh.median(),0),mad=f(Q.mad.mean(),1))
d=pd.read_csv("decomposition.csv")
J["decomp"]=dict(var_obs=f(np.var(d.observed,ddof=1),3),
  var_ref=f(np.var(d.ref_only,ddof=1),3),var_age=f(np.var(d.age_only,ddof=1),3),
  pct_ref=f(100*np.var(d.ref_only,ddof=1)/np.var(d.observed,ddof=1),0),
  pct_age=f(100*np.var(d.age_only,ddof=1)/np.var(d.observed,ddof=1),0))
print(json.dumps(J,indent=1)[:200])
# ---- merge in panel, sensitivity and rounded/day-converted values ----
import math
PJ=json.load(open("PANEL.json")); SJ=json.load(open("SENS2.json")); FJ=json.load(open("SENSITIVITY.json"))
QJ=json.load(open("OFFICIAL.json")); GJ=json.load(open("GENERALISATION.json"))
J["panel"]=PJ; J["sens"]=SJ; J["off"]=QJ; J["gen"]=GJ
J["round"]=dict(
  deaths_pct=f"{J['shares']['deaths_pct']:.2f}", yll_pct=f"{J['shares']['yll_pct']:.2f}",
  gain_yr=f"{J['gain_e0']['mean']:.2f}", gain_days=int(round(J['gain_e0']['mean']*365.25)),
  ded_yr=f"{J['d_edag']['mean']:.3f}", ded_days=int(round(J['d_edag']['mean']*365.25)),
  gain_lo=f"{J['gain_e0']['lo']:.2f}", gain_hi=f"{J['gain_e0']['hi']:.2f}",
  ded_lo=f"{J['d_edag']['lo']:.3f}", ded_hi=f"{J['d_edag']['hi']:.3f}",
  edag=f"{J['d_edag']['edag_mean']:.1f}")
json.dump(J,open("RESULTS.json","w"),indent=1)
print("\nmerged panel + sensitivity; gain %s yr (%d days); e-dagger %s yr (%d days)"%(
  J["round"]["gain_yr"],J["round"]["gain_days"],J["round"]["ded_yr"],J["round"]["ded_days"]))

# country table for the manuscript
NAME={"AT":"Austria","BE":"Belgium","BG":"Bulgaria","HR":"Croatia","CZ":"Czechia","DK":"Denmark",
"EE":"Estonia","FI":"Finland","FR":"France","DE":"Germany","EL":"Greece","HU":"Hungary","IS":"Iceland",
"IE":"Ireland","IT":"Italy","LV":"Latvia","LT":"Lithuania","LU":"Luxembourg","NL":"Netherlands",
"NO":"Norway","PL":"Poland","PT":"Portugal","SK":"Slovakia","SI":"Slovenia","ES":"Spain","SE":"Sweden",
"CH":"Switzerland","UK":"United Kingdom"}
T=pd.DataFrame({"Country":[NAME[g] for g in P.index],"Deaths":P.deaths.astype(int).values,
 "ASMR":P.asmr.round(2).values,"ASYR_nat":P.asyr_nat.round(0).astype(int).values,
 "ASYR_gbd":P.asyr_gbd.round(0).astype(int).values,"YPD_nat":P.ypd_nat.round(1).values,
 "YPD_gbd":P.ypd_gbd.round(1).values,"Ratio_gbd":P.ratio_gbd_nat.round(3).values,
 "ShareD":sh.sd.reindex(P.index).round(3).values,"ShareY":sh.sy.reindex(P.index).round(3).values,
 "SRatio":sh.ratio.reindex(P.index).round(2).values,
 "e0":P.e0.round(1).values,"gain_e0":P.gain_e0.round(3).values,
 "edag":P.edag.round(2).values,"d_edag":P.d_edag.round(4).values,
 "d_lo":U.dlo.reindex(P.index).round(4).values,"d_hi":U.dhi.reindex(P.index).round(4).values,
 "thresh":P.thresh.astype(int).values,"mad":P.mad.round(1).values,
 "below":(100*P.below).round(0).astype(int).values}).sort_values("Country")

# ---- extra columns from the sensitivity stage, so MS_TABLES.csv is fully pipeline-generated ----
OFF=pd.read_csv("official_reference_country.csv").set_index("geo")
Dq=pd.read_csv("sensitivity_countryyear.csv")
gq=Dq.groupby("geo").apply(lambda z: pd.Series({
    "ypd_bp":(z.ypd_bp*z.deaths).sum()/z.deaths.sum(),
    "asyr_bp":np.average(z.asyr_bp,weights=z.expo),"d_es":z.d_edag_es.mean()}),include_groups=False)
gq.index=[NAME[i] for i in gq.index]
T=T.merge(gq.reset_index().rename(columns={"index":"Country"}),on="Country")
T["YPD_bp"]=T.ypd_bp.round(1); T["ASYR_bp"]=T.asyr_bp.round(0).astype(int)
T["d_edag_es"]=T.d_es.round(3)
T["nyr"]=[PJ["years_per_country"][c] for c in T.Country]
_off=OFF.copy(); _off.index=[NAME[i] for i in _off.index]
T["YPD_off"]=[round(_off.ypd_off[c],1) for c in T.Country]
T["ASYR_off"]=[int(round(_off.asyr_off[c])) for c in T.Country]
T=T.drop(columns=["ypd_bp","asyr_bp","d_es"])
T.to_csv("MS_TABLES.csv",index=False)

# one-decimal age-standardised rates for the years-of-life-lost table, taken from the
# unrounded country summaries so that ranks quoted in the text and the printed values agree
_ad=OFF.copy(); _ad.index=[NAME[i] for i in _ad.index]
AD=pd.DataFrame({"Country":T.Country,
                 "ASYR_nat":[round(_ad.asyr_nat[c],1) for c in T.Country],
                 "ASYR_off":[round(_ad.asyr_off[c],1) for c in T.Country],
                 "ASYR_bp":[round(_ad.asyr_bp[c],1) for c in T.Country]})
AD.sort_values("Country").to_csv("MS_TABLES_ASYR.csv",index=False)
print("\n",T.to_string(index=False))
