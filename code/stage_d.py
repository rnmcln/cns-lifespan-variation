import pandas as pd, numpy as np, json
from engine import *
rng=np.random.default_rng(20260803)
R=pd.read_csv("results_countryyear.csv")
FR={"primary_2012_2018":lambda d:d[(d.year>=2012)&(d.year<=2018)],
    "maximal":lambda d:d,"excl_2020_21":lambda d:d[~d.year.isin([2020,2021])]}
rows=[]
for name,f in FR.items():
    for sx in ["T","M","F"]:
        d=f(R[R.sex==sx])
        c=d.groupby("geo").apply(lambda s: pd.Series({
          "deaths":s.deaths.sum(),"dall":s.dall.sum(),"nyr":s.year.nunique(),
          "asmr":np.average(s.asmr,weights=s.expo),
          "asyr_nat":np.average(s.asyr_nat,weights=s.expo),
          "asyr_gbd":np.average(s.asyr_gbd,weights=s.expo),
          "yll_nat":s.yll_nat.sum(),"yll_gbd":s.yll_gbd.sum(),
          "gain_e0":s.gain_e0.mean(),"d_edag":s.d_edag.mean(),"edag":s.edag.mean(),
          "e0":s.e0.mean(),"thresh":s.thresh.median(),
          "mad":np.average(s.mad,weights=s.deaths),
          "below":np.average(s.frac_below_thresh,weights=s.deaths)}),include_groups=False)
        c["ypd_nat"]=c.yll_nat/c.deaths; c["ypd_gbd"]=c.yll_gbd/c.deaths
        c["ratio_gbd_nat"]=c.ypd_gbd/c.ypd_nat
        rows.append(c.assign(frame=name,sex=sx).reset_index())
S=pd.concat(rows); S.to_csv("country_summary_all.csv",index=False)
P=S[(S.frame=="primary_2012_2018")&(S.sex=="T")].set_index("geo")
print("PRIMARY 2012-2018 both sexes: %d countries, %d deaths, %d country-years"%(
  len(P),P.deaths.sum(),len(R[(R.sex=='T')&(R.year.between(2012,2018))])))
print(" ASMR %.2f-%.2f | ASYR nat %.0f-%.0f | ASYR GBD %.0f-%.0f"%(
  P.asmr.min(),P.asmr.max(),P.asyr_nat.min(),P.asyr_nat.max(),P.asyr_gbd.min(),P.asyr_gbd.max()))
print(" YLL/death nat %.1f-%.1f | GBD/nat ratio %.3f-%.3f (mean %.3f)"%(
  P.ypd_nat.min(),P.ypd_nat.max(),P.ratio_gbd_nat.min(),P.ratio_gbd_nat.max(),P.ratio_gbd_nat.mean()))
print(" gain_e0 %.3f-%.3f mean %.3f | d_edag %.4f-%.4f mean %.4f"%(
  P.gain_e0.min(),P.gain_e0.max(),P.gain_e0.mean(),P.d_edag.min(),P.d_edag.max(),P.d_edag.mean()))
piv=S[S.sex=="T"].pivot_table(index="geo",columns="frame",values="asyr_nat")
rk=piv.rank(ascending=False)
print("\nSpearman rank corr of ASYR across frames:\n",rk.corr(method='spearman').round(3).to_string())
print("max rank shift primary vs maximal:",int((rk.primary_2012_2018-rk.maximal).abs().max()))
# ---- pooled shares (primary frame) ----
d=R[(R.sex=="T")&(R.year.between(2012,2018))]
print("\nPooled share of deaths %.3f%% | share of YLL(nat) %.3f%% | ratio %.2f"%(
 100*d.deaths.sum()/d.dall.sum(), 100*d.yll_nat.sum()/ (d.dall.sum()*0+1) *0+
 100*d.yll_nat.sum()/d.yll_nat.sum()*0+0, 0))
