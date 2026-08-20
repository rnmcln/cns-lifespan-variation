import pandas as pd, numpy as np, json, requests, time
from engine import *
BASE="https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/hlth_cd_aro"
AGES=["Y_LT1","Y1-4"]+[f"Y{a}-{a+4}" for a in range(5,95,5)]+["Y_GE95"]
AMAP=dict(zip(AGES,KEYS))
CAUSES=[("C70-C72","Malignant CNS neoplasms"),
        ("P","Perinatal conditions"),
        ("V_Y85","Transport accidents"),
        ("X60-X84_Y870","Intentional self-harm"),
        ("K70_K73_K74","Chronic liver disease"),
        ("C33_C34","Lung cancer"),
        ("I20-I25","Ischaemic heart disease"),
        ("I60-I69","Cerebrovascular disease"),
        ("F01_F03","Dementia")]
def pull(icd,ages):
    p=[("format","JSON"),("lang","EN"),("icd10",icd),("resid","TOT_RESID"),("sex","T")]+[("age",a) for a in ages]
    d=requests.get(BASE,params=p,timeout=300).json()
    dims,sizes=d["id"],d["size"]
    idx={k:{v:kk for kk,v in d["dimension"][k]["category"]["index"].items()} for k in dims}
    rows=[]
    for f,v in d["value"].items():
        n=int(f); c=[]
        for s in reversed(sizes): c.append(n%s); n//=s
        r={k:idx[k][i] for k,i in zip(dims,list(reversed(c)))}
        rows.append((r["geo"],int(r["time"]),r.get("age"),v))
    return pd.DataFrame(rows,columns=["geo","Year","age","deaths"])

# official GBD 2023 schedule
O=pd.read_csv("gbd/tmrlt_gbd2023.csv"); oa=O["Age"].values.astype(float); oe=O["Life Expectancy"].values
e_star=lambda a: np.interp(a,oa,oe)
STARTS=np.array([0,1]+list(range(5,95,5)),float); WID=np.array([1.,4.]+[5.]*18)
SUB_E=np.array([e_star(97.5),e_star(102.5),e_star(107.5),e_star(110.)])
OPEN=["95-99","100-104","105-109","110+"]
MAPc={"AT":"AUT","BE":"BEL","BG":"BGR","HR":"HRV","CZ":"CZE","DK":"DNK","EE":"EST","FI":"FIN",
"FR":"FRATNP","DE":"DEUTNP","EL":"GRC","HU":"HUN","IS":"ISL","IE":"IRL","IT":"ITA","LV":"LVA",
"LT":"LTU","LU":"LUX","NL":"NLD","NO":"NOR","PL":"POL","PT":"PRT","SK":"SVK","SI":"SVN",
"ES":"ESP","SE":"SWE","CH":"CHE","UK":"GBR_NP"}
def read_hmd(p):
    L=open(p).read().split("\n"); h=[i for i,l in enumerate(L) if l.strip().startswith("Year")][0]
    cols=L[h].split(); df=pd.DataFrame([l.split() for l in L[h+1:] if l.strip()],columns=cols)
    df["Year"]=df.Year.astype(int)
    for c in cols[2:]: df[c]=pd.to_numeric(df[c],errors="coerce")
    return df
OW={}
for geo,cc in MAPc.items():
    d=read_hmd(f"raw/{cc}_Deaths_5x1.txt"); d=d[d.Age.isin(OPEN)]
    for yr,s in d.groupby("Year"):
        w=s.set_index("Age")["Total"].reindex(OPEN).fillna(0).values.astype(float)
        OW[(geo,yr)]=w/w.sum() if w.sum()>0 else np.array([1.,0,0,0])

H=pd.read_csv("hmd_all.csv"); H=H[H.sex=="T"]
frame=pd.read_csv("results_countryyear.csv")
frame=frame[(frame.sex=="T")&(frame.year.between(2012,2018))][["geo","year"]].rename(columns={"year":"Year"})
MID=np.array([0.5,3.0]+[a+2.5 for a in range(5,95,5)]+[97.5])

out=[]; CTY=[]
for icd,name in CAUSES:
    A=pull(icd,AGES); A["k"]=A.age.map(AMAP); A=A[A.k.notna()]
    T=pull(icd,["TOTAL"]).rename(columns={"deaths":"tot"})[["geo","Year","tot"]]
    time.sleep(0.4)
    M=H.merge(A[["geo","Year","k","deaths"]],on=["geo","Year","k"],how="left").merge(frame,on=["geo","Year"])
    recs=[]
    for (geo,yr),s in M.groupby(["geo","Year"]):
        s=s.sort_values("ord")
        if len(s)!=21: continue
        tt=T[(T.geo==geo)&(T.Year==yr)].tot
        d=np.nan_to_num(s.deaths.values)
        if len(tt)==0 or abs(float(tt.iloc[0])-d.sum())>1e-6: continue   # same reconciliation rule
        if d.sum()==0: continue
        ax=s.ax.values
        lt=build_lt(s.mx.values,ax); e_nat=e_within(lt["ex"],ax)
        eo=np.empty(21); eo[:20]=e_star(STARTS+np.clip(ax[:20]/WID,0,1)*WID)
        eo[20]=float((OW[(geo,yr)]*SUB_E).sum())
        recs.append((geo,yr,d.sum(),(d*e_nat).sum(),(d*eo).sum(),float(np.average(MID,weights=d))))
    R=pd.DataFrame(recs,columns=["geo","Year","deaths","yll_nat","yll_off","mad"])
    if R.empty: print("  skipped",name); continue
    cty=R.groupby("geo").apply(lambda z:pd.Series({
        "ypd_nat":z.yll_nat.sum()/z.deaths.sum(),"ypd_off":z.yll_off.sum()/z.deaths.sum(),
        "deaths":z.deaths.sum(),
        "mad":float(np.average(z.mad,weights=z.deaths))}),include_groups=False)
    infl=100*(cty.ypd_off/cty.ypd_nat-1)
    c2=cty.reset_index(); c2.insert(0,"cause",name); c2.insert(0,"icd",icd)
    c2["infl_pct"]=infl.values
    CTY.append(c2)
    out.append(dict(icd=icd,cause=name,n_cy=len(R),n_ctry=R.geo.nunique(),deaths=int(R.deaths.sum()),
        mad=float(np.average(R.mad,weights=R.deaths)),
        ypd_nat=float(R.yll_nat.sum()/R.deaths.sum()),ypd_off=float(R.yll_off.sum()/R.deaths.sum()),
        infl_pooled=float(100*(R.yll_off.sum()/R.yll_nat.sum()-1)),
        infl_mean=float(infl.mean()),infl_lo=float(infl.min()),infl_hi=float(infl.max())))
    print("  %-26s n_cy=%3d deaths=%9d mad=%5.1f  nat=%5.1f off=%5.1f  infl=%5.1f%% (%.1f-%.1f)"%(
        name,out[-1]["n_cy"],out[-1]["deaths"],out[-1]["mad"],out[-1]["ypd_nat"],out[-1]["ypd_off"],
        out[-1]["infl_pooled"],out[-1]["infl_lo"],out[-1]["infl_hi"]))
G=pd.DataFrame(out).sort_values("mad")
G.to_csv("generalisation_causes.csv",index=False)
C=pd.concat(CTY,ignore_index=True)
C.to_csv("generalisation_country.csv",index=False)

# Common country frame: the countries available for every cause, so that the cross-cause
# comparison is not confounded by differing country coverage.
sets=[set(g.geo) for _,g in C.groupby("cause")]
common=sorted(set.intersection(*sets))
CC=C[C.geo.isin(common)]
rows=[]
for icd,name in CAUSES:
    z=CC[CC.icd==icd]
    if z.empty: continue
    rows.append(dict(icd=icd,cause=name,n_ctry=z.geo.nunique(),deaths=int(z.deaths.sum()),
        mad=float(np.average(z.mad,weights=z.deaths)),
        infl_mean=float(z.infl_pct.mean()),infl_lo=float(z.infl_pct.min()),
        infl_hi=float(z.infl_pct.max())))
GC=pd.DataFrame(rows).sort_values("mad")
GC.to_csv("generalisation_common.csv",index=False)
r=np.corrcoef(G.mad,G.infl_mean)[0,1]
rs=pd.Series(G.mad.values).corr(pd.Series(G.infl_mean.values),method="spearman")
J={"n_causes":int(len(G)),"pearson_mad_infl":round(float(r),3),"spearman_mad_infl":round(float(rs),3),
   "cns_infl":round(float(G[G.icd=="C70-C72"].infl_mean.iloc[0]),1),
   "cns_mad":round(float(G[G.icd=="C70-C72"].mad.iloc[0]),1),
   "min_infl_cause":G.iloc[G.infl_mean.values.argmin()]["cause"],
   "min_infl":round(float(G.infl_mean.min()),1),
   "max_infl_cause":G.iloc[G.infl_mean.values.argmax()]["cause"],
   "max_infl":round(float(G.infl_mean.max()),1),
   "mad_lo":round(float(G.mad.min()),1),"mad_hi":round(float(G.mad.max()),1),
   "total_deaths":int(G.deaths.sum())}
J["n_common_countries"]=len(common)
J["common_excluded"]=sorted(set(C.geo)-set(common))
J["common_spearman_mad_infl"]=round(float(pd.Series(GC.mad.values).corr(
    pd.Series(GC.infl_mean.values),method="spearman")),3)
J["common_min_infl"]=round(float(GC.infl_mean.min()),1)
J["common_max_infl"]=round(float(GC.infl_mean.max()),1)
J["common_min_infl_cause"]=GC.iloc[GC.infl_mean.values.argmin()]["cause"]
J["common_max_infl_cause"]=GC.iloc[GC.infl_mean.values.argmax()]["cause"]
J["common_mad_lo"]=round(float(GC.mad.min()),1)
J["common_mad_hi"]=round(float(GC.mad.max()),1)
J["common_cns_infl"]=round(float(GC[GC.icd=="C70-C72"].infl_mean.iloc[0]),1)
J["common_cns_mad"]=round(float(GC[GC.icd=="C70-C72"].mad.iloc[0]),1)
J["common_total_deaths"]=int(GC.deaths.sum())
J["common_order_matches_full"]=bool(list(GC.cause)==list(G.cause))
json.dump(J,open("GENERALISATION.json","w"),indent=1)
print("\n"+json.dumps(J,indent=1))
