import requests, re, os, time, sys
env=dict(l.strip().split("=",1) for l in open(".hmdenv") if "=" in l)
s=requests.Session(); s.headers.update({"User-Agent":"Mozilla/5.0 (research data retrieval)"})
L="https://www.mortality.org/Account/Login"
tok=re.search(r'name="__RequestVerificationToken"[^>]*value="([^"]+)"',s.get(L,timeout=60).text).group(1)
s.post(L,data={"Email":env["HMD_USER"],"Password":env["HMD_PASS"],"ReturnUrl":"",
               "__RequestVerificationToken":tok},timeout=60)

MAP={"AT":"AUT","BE":"BEL","BG":"BGR","HR":"HRV","CZ":"CZE","DK":"DNK","EE":"EST","FI":"FIN",
"FR":"FRATNP","DE":"DEUTNP","EL":"GRC","HU":"HUN","IS":"ISL","IE":"IRL","IT":"ITA","LV":"LVA",
"LT":"LTU","LU":"LUX","NL":"NLD","NO":"NOR","PL":"POL","PT":"PRT","SK":"SVK","SI":"SVN",
"ES":"ESP","SE":"SWE","CH":"CHE","UK":"GBR_NP"}
FILES=["Deaths_5x1","Exposures_5x1","bltper_5x1","fltper_5x1","mltper_5x1"]
os.makedirs("raw",exist_ok=True)
manifest=[]
import hashlib
for geo,cc in MAP.items():
    for f in FILES:
        url=f"https://www.mortality.org/File/GetDocument/hmd.v6/{cc}/STATS/{f}.txt"
        r=s.get(url,timeout=180)
        p=f"raw/{cc}_{f}.txt"
        if r.status_code==200 and r.content[:1]!=b"<":
            open(p,"wb").write(r.content)
            hdr=r.text.split("\n")[0].strip()
            lm=re.search(r"Last modified:\s*([^;]+)",hdr)
            manifest.append((geo,cc,f,r.status_code,len(r.content),
                             lm.group(1).strip() if lm else "",
                             hashlib.sha256(r.content).hexdigest()[:16]))
        else:
            manifest.append((geo,cc,f,r.status_code,len(r.content),"FAIL",""))
        time.sleep(0.3)
import csv
with open("hmd_manifest.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["geo","hmd_code","file","status","bytes","last_modified","sha256_16"])
    w.writerows(manifest)
bad=[m for m in manifest if m[3]!=200 or m[5]=="FAIL"]
print("files:",len(manifest),"| failures:",len(bad))
for b in bad: print("  FAIL",b[:4])
