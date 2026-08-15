import pandas as pd, numpy as np, re, json
KEYS=["0","1-4"]+[f"{a}-{a+4}" for a in range(5,95,5)]+["95+"]
OPEN=["95-99","100-104","105-109","110+"]
MAP={"AT":"AUT","BE":"BEL","BG":"BGR","HR":"HRV","CZ":"CZE","DK":"DNK","EE":"EST","FI":"FIN",
"FR":"FRATNP","DE":"DEUTNP","EL":"GRC","HU":"HUN","IS":"ISL","IE":"IRL","IT":"ITA","LV":"LVA",
"LT":"LTU","LU":"LUX","NL":"NLD","NO":"NOR","PL":"POL","PT":"PRT","SK":"SVK","SI":"SVN",
"ES":"ESP","SE":"SWE","CH":"CHE","UK":"GBR_NP"}
def read_hmd(p):
    L=open(p).read().split("\n"); h=[i for i,l in enumerate(L) if l.strip().startswith("Year")][0]
    cols=L[h].split(); df=pd.DataFrame([l.split() for l in L[h+1:] if l.strip()],columns=cols)
    df["Year"]=df.Year.astype(int)
    for c in cols[2:]: df[c]=pd.to_numeric(df[c],errors="coerce")
    df["k"]=df.Age.where(~df.Age.isin(OPEN),"95+")
    return df
E=pd.read_csv("eurostat_full.csv")

# ---- HMD assembly (T/M/F) -------------------------------------------------
LT={"T":"bltper_5x1","M":"mltper_5x1","F":"fltper_5x1"}
EXC={"T":"Total","M":"Male","F":"Female"}
H=[]
for geo,cc in MAP.items():
    ex=read_hmd(f"raw/{cc}_Exposures_5x1.txt"); dt=read_hmd(f"raw/{cc}_Deaths_5x1.txt")
    for sx in ["T","M","F"]:
        lt=read_hmd(f"raw/{cc}_{LT[sx]}.txt")
        e=ex.groupby(["Year","k"],as_index=False)[EXC[sx]].sum().rename(columns={EXC[sx]:"expo"})
        dd=dt.groupby(["Year","k"],as_index=False)[EXC[sx]].sum().rename(columns={EXC[sx]:"D_all_hmd"})
        L=lt[~lt.Age.isin(OPEN[1:])].copy()
        L.loc[L.Age=="95-99","ex"]=L.loc[L.Age=="95-99","Tx"]/L.loc[L.Age=="95-99","lx"]
        L=L[["Year","k","mx","qx","ax","lx","dx","Lx","Tx","ex"]]
        m=e.merge(dd,on=["Year","k"]).merge(L,on=["Year","k"]).assign(geo=geo,sex=sx)
        H.append(m)
H=pd.concat(H)
H["ord"]=H.k.map({k:i for i,k in enumerate(KEYS)})
H=H.dropna(subset=["ord"]).sort_values(["geo","sex","Year","ord"])
H.to_csv("hmd_all.csv",index=False)
print("HMD rows:",len(H),"| countries:",H.geo.nunique())

# ---- Reconciliation: Eurostat all-cause vs HMD all-cause ------------------
ea=E[(E.which=="all")&(E.sex=="T")].groupby(["geo","time"],as_index=False).deaths.sum().rename(columns={"time":"Year","deaths":"D_all_es"})
ha=H[H.sex=="T"].groupby(["geo","Year"],as_index=False).D_all_hmd.sum()
rec=ea.merge(ha,on=["geo","Year"])
rec["diff"]=rec.D_all_es-rec.D_all_hmd; rec["pct"]=100*rec["diff"]/rec.D_all_hmd
rec=rec[rec.geo.isin(MAP)]
print("\n--- Eurostat vs HMD all-cause deaths, %d country-years ---"%len(rec))
print("median |%% diff| %.3f | 95th pct %.3f | max %.3f (%s %d)"%(
 rec.pct.abs().median(),rec.pct.abs().quantile(.95),rec.pct.abs().max(),
 rec.loc[rec.pct.abs().idxmax(),"geo"],rec.loc[rec.pct.abs().idxmax(),"Year"]))
print(rec.groupby("geo").pct.agg(["mean","min","max"]).round(3).to_string())
rec.to_csv("reconciliation_allcause.csv",index=False)
