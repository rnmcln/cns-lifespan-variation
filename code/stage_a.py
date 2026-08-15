import requests, pandas as pd, numpy as np, json
H={"User-Agent":"eurostat/0.2"}
BASE="https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/hlth_cd_aro"
AGES=["Y_LT1","Y1-4"]+[f"Y{a}-{a+4}" for a in range(5,95,5)]+["Y_GE95"]
KEYS=["0","1-4"]+[f"{a}-{a+4}" for a in range(5,95,5)]+["95+"]
AMAP=dict(zip(AGES,KEYS))

def pull(sex,icd,ages):
    p=[("format","JSON"),("lang","EN"),("icd10",icd),("resid","TOT_RESID"),("sex",sex)]+[("age",a) for a in ages]
    d=requests.get(BASE,params=p,headers=H,timeout=300).json()
    dims,sizes=d["id"],d["size"]
    idx={k:{v:kk for kk,v in d["dimension"][k]["category"]["index"].items()} for k in dims}
    def dec(flat):
        f=int(flat); c=[]
        for s in reversed(sizes): c.append(f%s); f//=s
        return {k:idx[k][i] for k,i in zip(dims,list(reversed(c)))}
    rows=[dec(f)|{"deaths":v,"flag":d.get("status",{}).get(f,"")} for f,v in d["value"].items()]
    rows+=[dec(f)|{"deaths":np.nan,"flag":v} for f,v in d.get("status",{}).items() if f not in d["value"]]
    x=pd.DataFrame(rows); x["time"]=x.time.astype(int)
    x["k"]=x.age.map(AMAP)
    return x

out={}
for sx in ["T","M","F"]:
    out[("cns",sx)]=pull(sx,"C70-C72",AGES)
    out[("all",sx)]=pull(sx,"TOTAL",AGES)
    print(sx,"cns",len(out[("cns",sx)]),"all",len(out[("all",sx)]))
frames=[]
for (w,sx),v in out.items():
    v=v.assign(which=w,sex=sx); frames.append(v[["geo","time","k","age","deaths","flag","which","sex"]])
E=pd.concat(frames)
E.to_csv("eurostat_full.csv",index=False)
print("saved", E.shape)
