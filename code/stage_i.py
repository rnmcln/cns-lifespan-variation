import pandas as pd, numpy as np, json
from engine import *
O=pd.read_csv("gbd/tmrlt_gbd2023.csv")
oa=O["Age"].values.astype(float); oe=O["Life Expectancy"].values
e_star=lambda a: np.interp(a,oa,oe)
G19=dict(zip([0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95],
 [88.87,88.00,84.03,79.05,74.07,69.11,64.15,59.20,54.25,49.32,44.43,39.63,34.91,30.25,25.68,
  21.29,17.10,13.24,9.99,7.62,5.92]))
print("Official GBD 2023 TMRLT vs published GBD 2019, at exact ages")
for a in [0,1,5,20,40,60,75,80,85,90,95]:
    print("  age %3d  2023 %8.3f  2019 %6.2f  diff %+6.2f"%(a,e_star(a),G19[a],e_star(a)-G19[a]))
print("  ages above 95 present only in GBD 2023: "+", ".join("%d %.2f"%(a,e) for a,e in zip(oa[oa>95],oe[oa>95])))

H=pd.read_csv("hmd_all.csv"); E=pd.read_csv("eurostat_full.csv"); Rr=pd.read_csv("results_countryyear.csv")
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
STARTS=np.array([0,1]+list(range(5,95,5)),float); WID=np.array([1.,4.]+[5.]*18)
SUB_E=np.array([e_star(97.5),e_star(102.5),e_star(107.5),e_star(110.)])

cns=E[E.which=="cns"].rename(columns={"time":"Year"}); cns=cns[cns.k.notna()]
piv=cns.pivot_table(index=["geo","Year","sex","k"],values="deaths",dropna=False).reset_index()
ok=Rr[(Rr.sex=="T")&(Rr.year.between(2012,2018))][["geo","sex","year"]].rename(columns={"year":"Year"})
M=H.merge(piv,on=["geo","Year","sex","k"],how="inner").merge(ok,on=["geo","sex","Year"]).sort_values(["geo","Year","ord"])
imp=pd.read_csv("gbd/gbd_implied_reflife.csv"); imp=imp[~imp.is_agg]
MP={"<1 year":"0","2-4 years":"1-4","95+ years":"95+"}
for a in range(5,95,5): MP[f"{a}-{a+4} years"]=f"{a}-{a+4}"
E_DER=pd.Series({MP[k]:v for k,v in imp.groupby("age_name").e_implied.mean().items() if k in MP}).reindex(KEYS).values

rows=[];diffs=[]
for (geo,yr),s in M.groupby(["geo","Year"]):
    s=s.sort_values("ord"); d=np.nan_to_num(s.deaths.values); ax=s.ax.values; expo=s.expo.values
    eo=np.empty(21)
    eo[:20]=e_star(STARTS+np.clip(ax[:20]/WID,0,1)*WID)
    eo[20]=float((OW[(geo,yr)]*SUB_E).sum())
    rows.append({"geo":geo,"year":yr,"deaths":d.sum(),
       "ypd_off":float((d*eo).sum()/d.sum()),
       "asyr_off":float((np.divide(d*eo,expo,out=np.zeros_like(d),where=expo>0)*ESP).sum()*1e5)})
    diffs.append(eo-E_DER)
D=pd.DataFrame(rows); DF=np.array(diffs)
g=D.groupby("geo").apply(lambda z:pd.Series({"ypd_off":(z.ypd_off*z.deaths).sum()/z.deaths.sum(),
   "asyr_off":z.asyr_off.mean()}),include_groups=False)
S=pd.read_csv("country_summary_all.csv"); Pp=S[(S.frame=="primary_2012_2018")&(S.sex=="T")].set_index("geo")
bp=pd.read_csv("sensitivity_countryyear.csv").groupby("geo").asyr_bp.mean()
g["ypd_nat"]=Pp.ypd_nat; g["ypd_der"]=Pp.ypd_gbd; g["asyr_nat"]=Pp.asyr_nat; g["asyr_der"]=Pp.asyr_gbd; g["asyr_bp"]=bp
Q={}
Q["off_e0"]=round(float(e_star(0)),2); Q["off_e95"]=round(float(e_star(95)),2)
Q["d2019_e0"]=round(float(e_star(0)-88.87),2); Q["d2019_e95"]=round(float(e_star(95)-5.92),2)
Q["d2019_e90"]=round(float(e_star(90)-7.62),2)
Q["max_age"]=int(oa.max()); Q["e110"]=round(float(oe[-1]),2)
Q["off_pct_mean"]=round(float(100*((g.ypd_off/g.ypd_nat).mean()-1)),1)
Q["off_pct_lo"]=round(float(100*((g.ypd_off/g.ypd_nat).min()-1)),1)
Q["off_pct_hi"]=round(float(100*((g.ypd_off/g.ypd_nat).max()-1)),1)
Q["off_pct_lo_c"]=(g.ypd_off/g.ypd_nat).idxmin(); Q["off_pct_hi_c"]=(g.ypd_off/g.ypd_nat).idxmax()
Q["ypd_off_lo"]=round(float(g.ypd_off.min()),1); Q["ypd_off_hi"]=round(float(g.ypd_off.max()),1)
Q["asyr_off_lo"]=int(round(g.asyr_off.min())); Q["asyr_off_hi"]=int(round(g.asyr_off.max()))
rn=g.asyr_nat.rank(ascending=False); ro=g.asyr_off.rank(ascending=False)
Q["rank_shift_ge2_off"]=int((ro-rn).abs().ge(2).sum()); Q["rank_shift_max_off"]=int((ro-rn).abs().max())
Q["rho_off_der"]=round(float(ro.corr(g.asyr_der.rank(ascending=False),method="spearman")),3)
Q["rho_off_bp"]=round(float(ro.corr(g.asyr_bp.rank(ascending=False),method="spearman")),3)
Q["approx_err_mean_yrs"]=round(float(np.abs(DF).mean()),2)
Q["approx_err_max_yrs"]=round(float(np.abs(DF).max()),2)
Q["approx_err_ypd_mean_pct"]=round(float(100*((g.ypd_der/g.ypd_off)-1).mean()),1)
Q["approx_err_ypd_maxabs_pct"]=round(float(100*((g.ypd_der/g.ypd_off)-1).abs().max()),1)

# ---- Additional file 2: derived vs official vs GBD 2019, per age group ----
AGES=[0,1]+list(range(5,100,5))
off_exact=np.array([e_star(a) for a in AGES])
g19=np.array([88.87,88.00,84.03,79.05,74.07,69.11,64.15,59.20,54.25,49.32,44.43,39.63,
              34.91,30.25,25.68,21.29,17.10,13.24,9.99,7.62,5.92])
mean_applied=np.array([np.mean([r for r in DF[:,i]+E_DER[i]]) for i in range(21)])
AF=pd.DataFrame({"age_group":KEYS,"exact_age":AGES,
   "gbd2023_official_at_exact_age":np.round(off_exact,3),
   "gbd2019_published_at_exact_age":g19,
   "gbd2023_minus_gbd2019":np.round(off_exact-g19,3),
   "official_applied_to_deaths_mean":np.round(mean_applied,3),
   "derived_from_published_outputs":np.round(E_DER,3),
   "derived_minus_official_mean":np.round(E_DER-mean_applied,3)})
AF.to_csv("reference_schedules.csv",index=False)
print("\nAdditional file 2 written")
print(AF.to_string(index=False))

json.dump(Q,open("OFFICIAL.json","w"),indent=1); g.to_csv("official_reference_country.csv")
print("\n"+json.dumps(Q,indent=1))
