import pandas as pd, numpy as np, json
from engine import *
H=pd.read_csv("hmd_all.csv"); E=pd.read_csv("eurostat_full.csv")
cns=E[E.which=="cns"].rename(columns={"time":"Year"}); cns=cns[cns.k.notna()]
piv=cns.pivot_table(index=["geo","Year","sex","k"],values="deaths",dropna=False).reset_index()
R=pd.read_csv("results_countryyear.csv")
ok=R[(R.sex=="T")&(R.year.between(2012,2018))][["geo","sex","year"]].rename(columns={"year":"Year"})
M=H.merge(piv,on=["geo","Year","sex","k"],how="inner").merge(ok,on=["geo","sex","Year"]).sort_values(["geo","Year","ord"])
def stepwise_rev(mx,ax,Rx):
    prev=edagger({**build_lt(mx,ax),"ax":ax}); out=np.empty(len(mx)); m2=mx.copy()
    for i in range(len(mx)-1,-1,-1):
        m2[i]=mx[i]*Rx[i]; cur=edagger({**build_lt(m2,ax),"ax":ax}); out[i]=prev-cur; prev=cur
    return out
F=[];Rv=[];mad={96.0:{},97.5:{},100.0:{}}
for (geo,yr),s in M.groupby(["geo","Year"]):
    s=s.sort_values("ord"); d=np.nan_to_num(s.deaths.values); dall=s.D_all_hmd.values
    mx=s.mx.values; ax=s.ax.values
    Rx=1-np.divide(d,dall,out=np.zeros_like(d),where=dall>0)
    f,_=stepwise_contrib(mx,ax,Rx); F.append(f); Rv.append(stepwise_rev(mx,ax,Rx))
    for oa in mad: mad[oa].setdefault(geo,[]).append((mean_age_at_death(d,open_age=oa),d.sum()))
F=np.array(F); Rv=np.array(Rv)
b=pd.DataFrame({oa:{g:np.average([a for a,_ in x],weights=[w for _,w in x]) for g,x in v.items()} for oa,v in mad.items()})
S={"rev_sign_agreement_pct":round(100*float((np.sign(F)==np.sign(Rv)).mean()),2),
   "rev_max_abs_diff":float(np.abs(F-Rv).max()),
   "rev_total_max_diff":float(np.abs(F.sum(1)-Rv.sum(1)).max()),
   "openage_max_change_96":round(float((b[96.0]-b[97.5]).abs().max()),4),
   "openage_max_change_100":round(float((b[100.0]-b[97.5]).abs().max()),4)}
f=E[(E.which=="cns")&E.flag.notna()&(E.flag!="")]
S["suppressed_T"]=int((f.sex=="T").sum()); S["suppressed_M"]=int((f.sex=="M").sum())
S["suppressed_F"]=int((f.sex=="F").sum())
S["suppressed_countries_T"]=sorted(set(f[f.sex=="T"].geo))
S["suppressed_countries_sexspec_extra"]=sorted(set(f.geo)-set(f[f.sex=="T"].geo))
cy=R[R.year.between(2012,2018)].groupby(["geo","sex"]).size().unstack(fill_value=0)
S["sexspec_min_country_years_M"]=int(cy["M"].min()); S["sexspec_min_country_years_F"]=int(cy["F"].min())
S["sexspec_min_country_M"]=cy["M"].idxmin(); S["sexspec_min_country_F"]=cy["F"].idxmin()
json.dump(S,open("SENSITIVITY.json","w"),indent=1); print(json.dumps(S,indent=1))
