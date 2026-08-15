import pandas as pd, numpy as np, requests
from engine import *
H=pd.read_csv("hmd_all.csv"); E=pd.read_csv("eurostat_full.csv")
imp=pd.read_csv("gbd/gbd_implied_reflife.csv"); imp=imp[~imp.is_agg]
MP={"<1 year":"0","2-4 years":"1-4","95+ years":"95+"}
for a in range(5,95,5): MP[f"{a}-{a+4} years"]=f"{a}-{a+4}"
E_GBD=pd.Series({MP[k]:v for k,v in imp.groupby("age_name").e_implied.mean().items() if k in MP}).reindex(KEYS).values

cns=E[E.which=="cns"].rename(columns={"time":"Year"})
cns=cns[cns.k.notna()]
piv=cns.pivot_table(index=["geo","Year","sex","k"],values="deaths",dropna=False).reset_index()

BASE="https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/hlth_cd_aro"
def total_by(sex):
    p=[("format","JSON"),("lang","EN"),("icd10","C70-C72"),("resid","TOT_RESID"),("sex",sex),("age","TOTAL")]
    d=requests.get(BASE,params=p,timeout=180).json(); dims,sizes=d["id"],d["size"]
    idx={k:{v:kk for kk,v in d["dimension"][k]["category"]["index"].items()} for k in dims}
    out=[]
    for f,v in d["value"].items():
        n=int(f); c=[]
        for s in reversed(sizes): c.append(n%s); n//=s
        r={k:idx[k][i] for k,i in zip(dims,list(reversed(c)))}
        out.append((r["geo"],int(r["time"]),v))
    return pd.DataFrame(out,columns=["geo","Year","cns_total"]).assign(sex=sex)
TOT=pd.concat([total_by(s) for s in ["T","M","F"]])
TOT.to_csv("eurostat_cns_totals.csv",index=False)

M=H.merge(piv,on=["geo","Year","sex","k"],how="inner").sort_values(["geo","sex","Year","ord"])
M=M[(M.Year>=2011)&(M.Year<=2023)]
recs=[];contribs=[];excl=[]
starts=np.array([0,1]+list(range(5,100,5)),float)
for (geo,sx,yr),s in M.groupby(["geo","sex","Year"]):
    s=s.sort_values("ord")
    if len(s)!=21: excl.append((geo,sx,yr,"incomplete_age_rows")); continue
    tt=TOT[(TOT.geo==geo)&(TOT.Year==yr)&(TOT.sex==sx)].cns_total
    if len(tt)==0: excl.append((geo,sx,yr,"no_all_age_total")); continue
    raw=s.deaths.values; n_sup=int(np.isnan(raw).sum()); dcns=np.nan_to_num(raw)
    sf=float(tt.iloc[0]-dcns.sum())
    if abs(sf)>1e-6: excl.append((geo,sx,yr,f"unreconciled_{sf:.0f}")); continue
    dall=s.D_all_hmd.values; mx=s.mx.values; ax=s.ax.values; expo=s.expo.values
    if np.isnan(mx).any() or np.isnan(ax).any() or np.isnan(dall).any():
        excl.append((geo,sx,yr,"missing_hmd")); continue
    Rx=1-np.divide(dcns,dall,out=np.zeros_like(dcns),where=dall>0)
    if (Rx<0).any() or (Rx>1).any(): excl.append((geo,sx,yr,"invalid_Rx")); continue
    r=analyse(mx,ax,Rx,expo,dcns,E_GBD); ta,_=threshold_age(mx,ax)
    con,tc=stepwise_contrib(mx,ax,Rx)
    r.update(geo=geo,sex=sx,year=yr,thresh=ta,dall=dall.sum(),n_suppressed=n_sup,
             mad=mean_age_at_death(dcns),
             frac_below_thresh=float(dcns[starts<ta].sum()/dcns.sum()) if not np.isnan(ta) and dcns.sum()>0 else np.nan,
             stepwise_total=tc,Rx_min=float(Rx.min()))
    recs.append(r); contribs.append(pd.DataFrame({"geo":geo,"sex":sx,"year":yr,"k":KEYS,"contrib":con}))
R=pd.DataFrame(recs); C=pd.concat(contribs)
R.to_csv("results_countryyear.csv",index=False); C.to_csv("contrib_age.csv",index=False)
X=pd.DataFrame(excl,columns=["geo","sex","year","reason"]); X.to_csv("exclusions.csv",index=False)
t=R[R.sex=="T"]
print("country-years both-sex %d | M %d | F %d"%(len(t),(R.sex=='M').sum(),(R.sex=='F').sum()))
print("countries %d | CNS deaths %d"%(t.geo.nunique(),t.deaths.sum()))
print("exclusions (both-sex):\n",X[X.sex=='T'].reason.str.replace(r'_-?\d+','_n',regex=True).value_counts().to_string())
print("stepwise exactness max|resid| %.2e"%(R.d_edag-R.stepwise_total).abs().max())
print("threshold median %.0f range %.0f-%.0f"%(t.thresh.median(),t.thresh.min(),t.thresh.max()))
print("mean age CNS death %.1f-%.1f | frac below threshold %.3f-%.3f"%(t.mad.min(),t.mad.max(),t.frac_below_thresh.min(),t.frac_below_thresh.max()))
print("d_edag positive %d/%d country-years"%((t.d_edag>0).sum(),len(t)))
