import pandas as pd, numpy as np, json
from engine import *
H=pd.read_csv("hmd_all.csv"); E=pd.read_csv("eurostat_full.csv"); R=pd.read_csv("results_countryyear.csv")
rec=pd.read_csv("reconciliation_allcause.csv")
imp=pd.read_csv("gbd/gbd_implied_reflife.csv"); imp=imp[~imp.is_agg]
MP={"<1 year":"0","2-4 years":"1-4","95+ years":"95+"}
for a in range(5,95,5): MP[f"{a}-{a+4} years"]=f"{a}-{a+4}"
E_GBD=pd.Series({MP[k]:v for k,v in imp.groupby("age_name").e_implied.mean().items() if k in MP}).reindex(KEYS).values
cns=E[E.which=="cns"].rename(columns={"time":"Year"}); cns=cns[cns.k.notna()]
allc=E[E.which=="all"].rename(columns={"time":"Year"}); allc=allc[allc.k.notna()]
pc=cns.pivot_table(index=["geo","Year","sex","k"],values="deaths",dropna=False).reset_index()
pa=allc.pivot_table(index=["geo","Year","sex","k"],values="deaths",dropna=False).rename(columns={"deaths":"dall_es"}).reset_index()
ok=R[(R.sex=="T")&(R.year.between(2012,2018))][["geo","sex","year"]].rename(columns={"year":"Year"})
M=H.merge(pc,on=["geo","Year","sex","k"],how="inner").merge(pa,on=["geo","Year","sex","k"],how="left")\
   .merge(ok,on=["geo","sex","Year"]).sort_values(["geo","Year","ord"])
O={}
# ---------- 1. European best-practice reference (coherent life table from min m_x) ----------
mn=M.groupby("ord").mx.min().values
ax_med=M.groupby("ord").ax.median().values
lt_bp=build_lt(mn,ax_med); lt_bp["ax"]=ax_med
E_BP=e_within(lt_bp["ex"],ax_med)
O["bp_e0"]=round(float(lt_bp["ex"][0]),2)
O["bp_monotone_ex"]=bool(np.all(np.diff(lt_bp["ex"])<0))
O["bp_qx_valid"]=bool(np.all((lt_bp["qx"]>=0)&(lt_bp["qx"]<=1)))
O["bp_schedule"]={k:round(float(v),2) for k,v in zip(KEYS,E_BP)}
# ---------- 2/3. per-country-year quantities under both common refs + source sensitivity ----
rows=[];contribF=[];contribR=[]
for (geo,yr),s in M.groupby(["geo","Year"]):
    s=s.sort_values("ord"); d=np.nan_to_num(s.deaths.values)
    dh=s.D_all_hmd.values; de=s.dall_es.values
    mx=s.mx.values; ax=s.ax.values; expo=s.expo.values
    lt=build_lt(mx,ax); lt["ax"]=ax; e_nat=e_within(lt["ex"],ax)
    Rh=1-np.divide(d,dh,out=np.zeros_like(d),where=dh>0)
    r={"geo":geo,"year":yr,"deaths":d.sum(),
       "ypd_nat":float((d*e_nat).sum()/d.sum()),"ypd_gbd":float((d*E_GBD).sum()/d.sum()),
       "ypd_bp":float((d*E_BP).sum()/d.sum()),
       "asyr_bp":float((np.divide(d*E_BP,expo,out=np.zeros_like(d),where=expo>0)*ESP).sum()*1e5)}
    ltdh=build_lt(mx*Rh,ax); ltdh["ax"]=ax
    r["d_edag_hmd"]=edagger(lt)-edagger(ltdh); r["gain_hmd"]=ltdh["ex"][0]-lt["ex"][0]
    if not np.isnan(de).any():
        Re=1-np.divide(d,de,out=np.zeros_like(d),where=de>0); Re=np.clip(Re,0,1)
        ltde=build_lt(mx*Re,ax); ltde["ax"]=ax
        r["d_edag_es"]=edagger(lt)-edagger(ltde); r["gain_es"]=ltde["ex"][0]-lt["ex"][0]
    # symmetric age decomposition
    f,_=stepwise_contrib(mx,ax,Rh)
    prev=edagger(lt); rv=np.empty(21); m2=mx.copy()
    for i in range(20,-1,-1):
        m2[i]=mx[i]*Rh[i]; cur=edagger({**build_lt(m2,ax),"ax":ax}); rv[i]=prev-cur; prev=cur
    contribF.append(pd.DataFrame({"geo":geo,"year":yr,"k":KEYS,"fwd":f,"rev":rv,"sym":(f+rv)/2}))
    rows.append(r)
D=pd.DataFrame(rows); CS=pd.concat(contribF)
D.to_csv("sensitivity_countryyear.csv",index=False); CS.to_csv("contrib_symmetric.csv",index=False)
g=D.groupby("geo").apply(lambda s:pd.Series({
  "ypd_nat":np.average(s.ypd_nat,weights=s.deaths),"ypd_gbd":np.average(s.ypd_gbd,weights=s.deaths),
  "ypd_bp":np.average(s.ypd_bp,weights=s.deaths),"asyr_bp":s.asyr_bp.mean(),
  "d_hmd":s.d_edag_hmd.mean(),"d_es":s.d_edag_es.mean(),
  "g_hmd":s.gain_hmd.mean(),"g_es":s.gain_es.mean()}),include_groups=False)
S=pd.read_csv("country_summary_all.csv"); P=S[(S.frame=="primary_2012_2018")&(S.sex=="T")].set_index("geo")
g["asyr_nat"]=P.asyr_nat; g["asyr_gbd"]=P.asyr_gbd
O["bp_ratio_mean"]=round(float((g.ypd_bp/g.ypd_nat).mean()),3)
O["bp_pct_lo"]=round(float(100*((g.ypd_bp/g.ypd_nat).min()-1)),1)
O["bp_pct_hi"]=round(float(100*((g.ypd_bp/g.ypd_nat).max()-1)),1)
O["bp_pct_mean"]=round(float(100*((g.ypd_bp/g.ypd_nat).mean()-1)),1)
rn=g.asyr_nat.rank(ascending=False); rg=g.asyr_gbd.rank(ascending=False); rb=g.asyr_bp.rank(ascending=False)
O["rank_shift_ge2_gbd"]=int((rg-rn).abs().ge(2).sum()); O["rank_shift_max_gbd"]=int((rg-rn).abs().max())
O["rank_shift_ge2_bp"]=int((rb-rn).abs().ge(2).sum()); O["rank_shift_max_bp"]=int((rb-rn).abs().max())
O["rho_gbd_bp"]=round(float(rg.corr(rb,method="spearman")),3)
# source sensitivity
O["src_d_edag_hmd"]=round(float(g.d_hmd.mean()),4); O["src_d_edag_es"]=round(float(g.d_es.mean()),4)
O["src_gain_hmd"]=round(float(g.g_hmd.mean()),3); O["src_gain_es"]=round(float(g.g_es.mean()),3)
O["src_n_countries_es"]=int(g.d_es.notna().sum())
O["src_max_abs_diff_edag"]=round(float((g.d_hmd-g.d_es).abs().max()),4)
r2=rec[rec.geo.isin(g.index)&rec.Year.between(2012,2018)]
lowdisc=set(r2[r2.pct.abs()<1].geo.unique())&set(r2.groupby("geo").pct.apply(lambda s:(s.abs()<1).all())[lambda x:x].index)
O["src_lowdisc_n"]=len(lowdisc)
O["src_lowdisc_d_edag"]=round(float(g.loc[list(lowdisc),"d_hmd"].mean()),4)
O["src_lowdisc_positive"]=int((g.loc[list(lowdisc),"d_hmd"]>0).sum())
# ---------- 4. threshold perturbation size ----------
th={}
for dl in [0.01,0.001,0.0001]:
    v=[]
    for (geo,yr),s in M.groupby(["geo","Year"]):
        s=s.sort_values("ord"); a,_=threshold_age(s.mx.values,s.ax.values,delta=dl); v.append(a)
    th[str(dl)]={"median":float(np.median(v)),"min":float(np.min(v)),"max":float(np.max(v))}
O["threshold_by_delta"]=th
O["threshold_identical_across_delta"]=bool(th["0.01"]==th["0.001"]==th["0.0001"])
# ---------- 5. symmetric vs forward profiles ----------
O["sym_vs_fwd_max_diff"]=round(float((CS.sym-CS.fwd).abs().max()),6)
# ---------- 6. formal variance decomposition (exact identity + Shapley split) ----------
w=M.groupby(["geo","k"]).deaths.sum().unstack().reindex(columns=KEYS).fillna(0)
w=w.div(w.sum(axis=1),axis=0)
en=[]
for geo,s in M.groupby("geo"):
    ee=[]
    for yr,ss in s.groupby("Year"):
        ss=ss.sort_values("ord"); lt=build_lt(ss.mx.values,ss.ax.values)
        ee.append(e_within(lt["ex"],ss.ax.values))
    en.append(pd.Series(np.mean(ee,axis=0),index=KEYS,name=geo))
en=pd.DataFrame(en)
wb=w.mean(axis=0); eb=en.mean(axis=0)
y=(w.values*en.values).sum(axis=1); ybar=y.mean()
A=((w.values-wb.values)*eb.values).sum(axis=1)
Rr=(wb.values*(en.values-eb.values)).sum(axis=1)
I=((w.values-wb.values)*(en.values-eb.values)).sum(axis=1)
O["ident_max_resid"]=float(np.abs((y-ybar)-(A+Rr+I)-((y-ybar).mean()-(A+Rr+I).mean())).max())
shA=A+I/2; shR=Rr+I/2
O["var_total"]=round(float(np.var(y,ddof=1)),3)
O["var_shapley_age"]=round(float(np.var(shA,ddof=1)),3)
O["var_shapley_ref"]=round(float(np.var(shR,ddof=1)),3)
O["cov_shapley"]=round(float(2*np.cov(shA,shR,ddof=1)[0,1]),3)
O["var_sum_check"]=round(float(np.var(shA,ddof=1)+np.var(shR,ddof=1)+2*np.cov(shA,shR,ddof=1)[0,1]),3)
O["sd_shapley_age"]=round(float(np.std(shA,ddof=1)),3); O["sd_shapley_ref"]=round(float(np.std(shR,ddof=1)),3)
json.dump(O,open("SENS2.json","w"),indent=1); print(json.dumps(O,indent=1)[:3000])
