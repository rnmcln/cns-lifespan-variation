import pandas as pd, numpy as np, json
from engine import *
rng=np.random.default_rng(20260803); NREP=10000
H=pd.read_csv("hmd_all.csv"); E=pd.read_csv("eurostat_full.csv")
imp=pd.read_csv("gbd/gbd_implied_reflife.csv"); imp=imp[~imp.is_agg]
MP={"<1 year":"0","2-4 years":"1-4","95+ years":"95+"}
for a in range(5,95,5): MP[f"{a}-{a+4} years"]=f"{a}-{a+4}"
E_GBD=pd.Series({MP[k]:v for k,v in imp.groupby("age_name").e_implied.mean().items() if k in MP}).reindex(KEYS).values
cns=E[E.which=="cns"].rename(columns={"time":"Year"}); cns=cns[cns.k.notna()]
piv=cns.pivot_table(index=["geo","Year","sex","k"],values="deaths",dropna=False).reset_index()
M=H.merge(piv,on=["geo","Year","sex","k"],how="inner")
ok=pd.read_csv("results_countryyear.csv")[["geo","sex","year"]].rename(columns={"year":"Year"})
M=M.merge(ok,on=["geo","sex","Year"]).sort_values(["geo","sex","Year","ord"])
allc=[];unc=[];agedist=[]
for (geo,sx,yr),s in M.groupby(["geo","sex","Year"]):
    s=s.sort_values("ord")
    dcns=np.nan_to_num(s.deaths.values); dall=s.D_all_hmd.values
    mx=s.mx.values; ax=s.ax.values; expo=s.expo.values
    lt=build_lt(mx,ax); lt["ax"]=ax; e_nat=e_within(lt["ex"],ax)
    allc.append(dict(geo=geo,sex=sx,year=yr,yll_all_nat=float((dall*e_nat).sum()),
                     yll_cns_nat=float((dcns*e_nat).sum()),dall=dall.sum(),dcns=dcns.sum()))
    if sx=="T" and 2012<=yr<=2018:
        agedist.append(pd.DataFrame({"geo":geo,"year":yr,"k":KEYS,"d":dcns,"e_nat":e_nat}))
        Rx0=1-np.divide(dcns,dall,out=np.zeros_like(dcns),where=dall>0)
        b=analyse(mx,ax,Rx0,expo,dcns,E_GBD); tot=dcns.sum()
        # Poisson total with multinomial allocation is exactly equivalent to
        # independent Poisson draws per age group with mean equal to the observed count
        D2=rng.poisson(np.tile(dcns,(NREP,1))).astype(float)
        R2=np.clip(1-np.divide(D2,dall,out=np.zeros_like(D2),where=dall>0),0,1)
        dx,lx,ex=build_lt_2d(mx[None,:]*R2,ax)
        de=edagger(lt)-edagger_2d(dx,lx,ex,ax); g=ex[:,0]-lt["ex"][0]
        yp=(D2*e_nat).sum(axis=1)/np.maximum(D2.sum(axis=1),1)
        unc.append(dict(geo=geo,year=yr,gain_e0=b["gain_e0"],
            gain_lo=np.percentile(g,2.5),gain_hi=np.percentile(g,97.5),
            d_edag=b["d_edag"],ed_lo=np.percentile(de,2.5),ed_hi=np.percentile(de,97.5),
            ypd=b["yll_nat"]/max(tot,1),ypd_lo=np.percentile(yp,2.5),ypd_hi=np.percentile(yp,97.5)))
A=pd.DataFrame(allc); A.to_csv("allcause_yll.csv",index=False)
U=pd.DataFrame(unc); U.to_csv("uncertainty_countryyear.csv",index=False)
AD=pd.concat(agedist); AD.to_csv("agedist_primary.csv",index=False)
a=A[(A.sex=="T")&(A.year.between(2012,2018))]
sd=100*a.dcns.sum()/a.dall.sum(); sy=100*a.yll_cns_nat.sum()/a.yll_all_nat.sum()
print("PRIMARY pooled: deaths %.3f%% | all-cause YLL %.3f%% | ratio %.2f"%(sd,sy,sy/sd))
cty=a.groupby("geo").apply(lambda s:pd.Series({"sd":100*s.dcns.sum()/s.dall.sum(),
   "sy":100*s.yll_cns_nat.sum()/s.yll_all_nat.sum()}),include_groups=False)
cty["ratio"]=cty.sy/cty.sd; cty.to_csv("shares_primary.csv")
print(" country ratio range %.2f-%.2f"%(cty.ratio.min(),cty.ratio.max()))
uc=U.groupby("geo").agg(gain=("gain_e0","mean"),glo=("gain_lo","mean"),ghi=("gain_hi","mean"),
   ded=("d_edag","mean"),dlo=("ed_lo","mean"),dhi=("ed_hi","mean"),
   ypd=("ypd","mean"),ylo=("ypd_lo","mean"),yhi=("ypd_hi","mean"))
uc.to_csv("uncertainty_country.csv")
print("NREP=%d | CI excludes zero: e-dagger %d/%d, gain %d/%d"%(NREP,(uc.dlo>0).sum(),len(uc),(uc.glo>0).sum(),len(uc)))
w=AD.pivot_table(index="geo",columns="k",values="d",aggfunc="sum").reindex(columns=KEYS)
w=w.div(w.sum(axis=1),axis=0)
en=AD.pivot_table(index="geo",columns="k",values="e_nat",aggfunc="mean").reindex(columns=KEYS)
wbar=w.mean(axis=0); ebar=en.mean(axis=0)
pd.DataFrame({"geo":w.index,"observed":(w.values*en.values).sum(1),
  "ref_only":(np.tile(wbar.values,(len(w),1))*en.values).sum(1),
  "age_only":(w.values*np.tile(ebar.values,(len(w),1))).sum(1)}).to_csv("decomposition.csv",index=False)
