import numpy as np, pandas as pd

KEYS = ["0","1-4"]+[f"{a}-{a+4}" for a in range(5,95,5)]+["95+"]
N    = np.array([1.,4.]+[5.]*18+[np.nan])
MID  = np.array([0.5,3.0]+[a+2.5 for a in range(5,95,5)]+[np.nan])   # 95+ handled separately
ESP  = np.array([1000,4000,5500,5500,5500,6000,6000,6500,7000,7000,7000,7000,6500,6000,
                 5500,5000,4000,2500,1500,800,200],float); ESP = ESP/ESP.sum()

def build_lt(mx, ax):
    """Abridged period life table from age-specific rates and observed a_x.
    Open interval (last): L = l/m, e = 1/m."""
    n = N.copy(); k = len(mx)
    qx = np.empty(k); qx[:-1] = n[:-1]*mx[:-1]/(1+(n[:-1]-ax[:-1])*mx[:-1])
    qx = np.clip(qx, 0, 1); qx[-1] = 1.0
    lx = np.empty(k); lx[0] = 1.0
    for i in range(k-1): lx[i+1] = lx[i]*(1-qx[i])
    dx = np.empty(k); dx[:-1] = lx[:-1]-lx[1:]; dx[-1] = lx[-1]
    Lx = np.empty(k)
    Lx[:-1] = n[:-1]*lx[1:] + ax[:-1]*dx[:-1]
    Lx[-1] = lx[-1]/mx[-1] if mx[-1] > 0 else 0.0
    Tx = np.flip(np.cumsum(np.flip(Lx))); ex = np.where(lx > 0, Tx/lx, 0.0)
    return dict(qx=qx, lx=lx, dx=dx, Lx=Lx, Tx=Tx, ex=ex)

def e_within(ex, ax):
    """Reference remaining expectancy applied to deaths inside each interval:
    linear interpolation between e_x and e_{x+n} at position a_x/n.
    Open interval uses e at exact age 95."""
    en = np.append(ex[1:], ex[-1]); f = np.clip(ax/np.where(np.isnan(N),1,N), 0, 1)
    out = ex + (en-ex)*f; out[-1] = ex[-1]; return out

def edagger(lt):
    """Life disparity: mean remaining life expectancy at death, computed
    entirely from the life table supplied (its own d_x, e_x, a_x)."""
    return float((lt["dx"]*e_within(lt["ex"], lt["ax"])).sum()/lt["lx"][0])

def analyse(mx, ax, Rx, expo, dcns, e_gbd):
    """Rx = 1 - (cause deaths / all-cause deaths), both from the same source."""
    lt  = build_lt(mx, ax); lt["ax"] = ax
    ltd = build_lt(mx*Rx, ax); ltd["ax"] = ax
    e_nat = e_within(lt["ex"], ax)
    rate  = np.divide(dcns, expo, out=np.zeros_like(dcns,dtype=float), where=expo>0)
    D     = dcns.sum()
    return dict(
        deaths=D, expo=expo.sum(),
        asmr=float((rate*ESP).sum()*1e5),
        yll_nat=float((dcns*e_nat).sum()), yll_gbd=float((dcns*e_gbd).sum()),
        asyr_nat=float((np.divide(dcns*e_nat,expo,out=np.zeros_like(dcns,dtype=float),
                                  where=expo>0)*ESP).sum()*1e5),
        asyr_gbd=float((np.divide(dcns*e_gbd,expo,out=np.zeros_like(dcns,dtype=float),
                                  where=expo>0)*ESP).sum()*1e5),
        e0=float(lt["ex"][0]), e0_del=float(ltd["ex"][0]),
        gain_e0=float(ltd["ex"][0]-lt["ex"][0]),
        edag=edagger(lt), edag_del=edagger(ltd),
        d_edag=float(edagger(lt)-edagger(ltd)))

def threshold_age(mx, ax, delta=0.01):
    """Numerical threshold age: the youngest age group at which a small
    proportional REDUCTION in mortality raises e-dagger. Returns the lower
    bound of that age group, or nan if no sign change occurs."""
    base = edagger({**build_lt(mx,ax), "ax":ax})
    sign = np.empty(len(mx))
    for i in range(len(mx)):
        m2 = mx.copy(); m2[i] = m2[i]*(1-delta)
        sign[i] = edagger({**build_lt(m2,ax), "ax":ax}) - base
    pos = np.where(sign > 0)[0]
    if len(pos) == 0: return np.nan, sign
    starts = np.array([0,1]+list(range(5,100,5)),float)
    return float(starts[pos[0]]), sign

def stepwise_contrib(mx, ax, Rx):
    """Exact stepwise decomposition of the finite cause-deletion effect on
    e-dagger: delete the cause cumulatively from the youngest age upward and
    take successive differences. Contributions sum exactly to the total."""
    base = edagger({**build_lt(mx,ax), "ax":ax}); prev = base
    out = np.empty(len(mx)); m2 = mx.copy()
    for i in range(len(mx)):
        m2[i] = mx[i]*Rx[i]
        cur = edagger({**build_lt(m2,ax), "ax":ax})
        out[i] = prev-cur; prev = cur
    return out, base-prev

def mean_age_at_death(d, open_age=97.5):
    mid = MID.copy(); mid[-1] = open_age
    return float(np.average(mid, weights=d)) if d.sum() > 0 else np.nan


# --- vectorised over replicates (rows) for resampling; identical algebra to build_lt ---
def build_lt_2d(MX, ax):
    n=N.copy(); k=MX.shape[1]
    qx=np.empty_like(MX)
    qx[:,:-1]=n[:-1]*MX[:,:-1]/(1+(n[:-1]-ax[:-1])*MX[:,:-1])
    qx=np.clip(qx,0,1); qx[:,-1]=1.0
    lx=np.empty_like(MX); lx[:,0]=1.0
    lx[:,1:]=np.cumprod(1-qx[:,:-1],axis=1)
    dx=np.empty_like(MX); dx[:,:-1]=lx[:,:-1]-lx[:,1:]; dx[:,-1]=lx[:,-1]
    Lx=np.empty_like(MX)
    Lx[:,:-1]=n[:-1]*lx[:,1:]+ax[:-1]*dx[:,:-1]
    Lx[:,-1]=np.where(MX[:,-1]>0, lx[:,-1]/np.where(MX[:,-1]>0,MX[:,-1],1), 0.0)
    Tx=np.flip(np.cumsum(np.flip(Lx,axis=1),axis=1),axis=1)
    ex=np.where(lx>0,Tx/lx,0.0)
    return dx,lx,ex

def e_within_2d(EX, ax):
    en=np.concatenate([EX[:,1:],EX[:,-1:]],axis=1)
    f=np.clip(ax/np.where(np.isnan(N),1,N),0,1)
    out=EX+(en-EX)*f; out[:,-1]=EX[:,-1]; return out

def edagger_2d(dx,lx,ex,ax):
    return (dx*e_within_2d(ex,ax)).sum(axis=1)/lx[:,0]
