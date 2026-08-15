import pandas as pd, numpy as np, json
from engine import *
H=pd.read_csv("hmd_all.csv"); E=pd.read_csv("eurostat_full.csv")
R=pd.read_csv("results_countryyear.csv"); G={}
NM={"AT":"Austria","BE":"Belgium","BG":"Bulgaria","HR":"Croatia","CZ":"Czechia","DK":"Denmark",
"EE":"Estonia","FI":"Finland","FR":"France","DE":"Germany","EL":"Greece","HU":"Hungary","IS":"Iceland",
"IE":"Ireland","IT":"Italy","LV":"Latvia","LT":"Lithuania","LU":"Luxembourg","NL":"Netherlands",
"NO":"Norway","PL":"Poland","PT":"Portugal","SK":"Slovakia","SI":"Slovenia","ES":"Spain","SE":"Sweden",
"CH":"Switzerland","UK":"United Kingdom"}
YRS=list(range(2012,2019))
t=R[(R.sex=="T")&(R.year.between(2012,2018))]
have=t.groupby("geo").year.apply(set)
G["panel_potential"]=28*7
G["panel_actual"]=int(len(t))
miss=[]
for g,s in have.items():
    for y in YRS:
        if y not in s: miss.append({"country":NM[g],"year":y})
G["panel_missing"]=miss
G["years_per_country"]={NM[g]:len(s) for g,s in have.items()}
G["balanced_countries"]=[NM[g] for g,s in have.items() if len(s)==7]
G["n_balanced"]=len(G["balanced_countries"])
# later period: countries complete 2012-2021
tt=R[(R.sex=="T")&(R.year.between(2012,2021))]
h2=tt.groupby("geo").year.apply(set)
lat=[g for g,s in h2.items() if set(range(2012,2022))<=s]
G["later_period_countries"]=[NM[g] for g in lat]; G["n_later"]=len(lat)
G["later_period_years"]="2012-2021"
json.dump(G,open("PANEL.json","w"),indent=1)
print("potential %d actual %d missing %d"%(G["panel_potential"],G["panel_actual"],len(miss)))
print("missing:",miss)
print("balanced countries (7/7): %d"%G["n_balanced"])
print("complete 2012-2021: %d ->"%G["n_later"],[NM[g] for g in lat])
