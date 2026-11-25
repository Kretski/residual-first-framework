#!/usr/bin/env python3
"""
Shammass et al. (2012) R0 implementation, specification v1.0.3.

Safety gates:
  * default stage is audit (hashes, counts, WPD conversion reproduction only)
  * curve-only computes §3/§4.3/§6 diagnostics without reading experimental f
  * synthetic gate requires --authorize-gate
  * real residual analysis requires --authorize-real AND a gate verdict file
No stage silently runs the synthetic gate or opens experiment-minus-GPE residuals.

Dependencies: numpy, scipy. Python >=3.10.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.stats import f as f_dist, beta

SPEC = "SHAMMASS_R0_PREANALYSIS_SPEC_v1.0.3"
FILES = {
    "exp_large": ("exp_large.csv", 16, "203637ca80adee49d01bc01a2c9456db5fb5bcbf9e8efac768e8e2537adbb39b"),
    "exp_small": ("exp_small.csv", 9, "a39239136a7011a794a8e72e2f50f0176fef4e04eb5362076d6ed512e190aa02"),
    "gpe_large": ("gpe_large.csv", 38, "ddb5e53199d8048aba8ccd08c6905b44151c63e573d11ff7514657f703d2e4d9"),
    "gpe_small": ("gpe_small.csv", 28, "c19213b54d9c5433265140f8cad1517102733bc5031a581b3bc7e7c26b500066"),
}
PAIR = {"large": ("exp_large", "gpe_large"), "small": ("exp_small", "gpe_small")}
AK, AF, RATIO = 0.0065830, 0.0032438, 2.0294
ALPHA, COND_MAX = 0.01, 1e8
SEED = 20261125
DELTAS = (1, 2, 4, 8)
MS = (2, 3, 5, 8)
SHAPES = (1, 2, 3)
SIGNS = (-1, 1)
N_NULL, N_SIGNAL = 10_000, 500

class R0Error(RuntimeError): pass

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def json_safe(value):
    """Convert NumPy values and non-finite floats to strict-JSON values."""
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        x = float(value)
        return x if math.isfinite(x) else None
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value

def json_text(value, **kwargs):
    """Serialize strict RFC-compatible JSON (NaN/Infinity are forbidden)."""
    return json.dumps(json_safe(value), allow_nan=False, **kwargs)

def load_csv(path: Path, columns: set[str] | None = None) -> dict[str, np.ndarray]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f)
        required = {"k_um_inv", "omega_over_2pi_kHz", "px", "py"}
        if not required.issubset(rd.fieldnames or []):
            raise R0Error(f"{path}: required columns absent")
        rows = list(rd)
    out = {}
    use = required if columns is None else set(columns)
    if not use.issubset(required): raise R0Error(f"Unsupported requested columns: {use-required}")
    for c in use:
        try: out[c] = np.asarray([float(r[c]) for r in rows], dtype=float)
        except Exception as e: raise R0Error(f"{path}: invalid numeric column {c}: {e}")
        if not np.all(np.isfinite(out[c])): raise R0Error(f"{path}: non-finite {c}")
    out["_n"] = np.asarray([len(rows)])
    return out

@dataclass
class Calibration:
    A: np.ndarray
    c: np.ndarray
    pix: np.ndarray
    data: np.ndarray

    def pixel_to_data(self, xy: np.ndarray) -> np.ndarray:
        return np.asarray(xy) @ self.A.T + self.c

    def data_to_pixel(self, uv: np.ndarray) -> np.ndarray:
        return (np.asarray(uv) - self.c) @ np.linalg.inv(self.A).T

def calibration_from_wpd(path: Path, jitter: np.ndarray | None = None) -> Calibration:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if len(obj.get("axesColl", [])) != 1: raise R0Error("Expected exactly one axesColl")
    ax = obj["axesColl"][0]
    if ax.get("type") != "XYAxes" or ax.get("isLogX") or ax.get("isLogY") or ax.get("noRotation") is not False:
        raise R0Error("Expected linear XYAxes with noRotation=false")
    cps = ax.get("calibrationPoints", [])
    if len(cps) != 4: raise R0Error("Expected calibration points X1,X2,Y1,Y2")
    pix = np.array([[float(q["px"]), float(q["py"])] for q in cps], dtype=float)
    if jitter is not None:
        if np.asarray(jitter).shape != (4, 2): raise R0Error("Calibration jitter must be 4x2")
        pix = pix + jitter
    data = np.array([[float(cps[0]["dx"]), 0.0],
                     [float(cps[1]["dx"]), 0.0],
                     [0.0, float(cps[2]["dy"])],
                     [0.0, float(cps[3]["dy"])]], dtype=float)
    xmin, xmax = data[0,0], data[1,0]
    ymin, ymax = data[2,1], data[3,1]
    x1,x2,x3,x4 = pix[:,0]; y1,y2,y3,y4 = pix[:,1]
    dat_mat = np.array([[xmin-xmax, 0.0], [0.0, ymin-ymax]])
    pix_mat = np.array([[x1-x2, x3-x4], [y1-y2, y3-y4]])
    if abs(np.linalg.det(pix_mat)) < 1e-15: raise R0Error("Singular calibration")
    A = dat_mat @ np.linalg.inv(pix_mat)
    c = np.array([xmin-A[0,0]*x1-A[0,1]*y1,
                  ymin-A[1,0]*x3-A[1,1]*y3])
    return Calibration(A, c, pix, data)

def locate(root: Path, filename: str) -> Path:
    candidates = list(root.rglob(filename))
    if len(candidates) != 1: raise R0Error(f"Expected one {filename}, found {len(candidates)}")
    return candidates[0]

def load_inputs(root: Path, verify_hashes: bool = True, curve_only: bool = False):
    exports = root / "realdata/shammass2012/unified_v1/exports"
    if not exports.exists():
        # permit passing realdata/shammass2012 or unified_v1 directly
        candidates = [root / "unified_v1/exports", root / "exports"]
        exports = next((x for x in candidates if x.exists()), exports)
    if not exports.exists(): raise R0Error(f"Cannot locate unified_v1/exports below {root}")
    ds = {}
    for key, (fn, n, digest) in FILES.items():
        path = exports / fn
        if not path.is_file(): raise R0Error(f"Missing input: {path}")
        actual = sha256(path)
        if verify_hashes and actual != digest:
            raise R0Error(f"SHA-256 mismatch for {fn}: expected {digest}, got {actual}")
        cols = {"k_um_inv"} if curve_only and key.startswith("exp_") else None
        d = load_csv(path, cols)
        if len(d["k_um_inv"]) != n: raise R0Error(f"Point count mismatch for {fn}: {len(d['k_um_inv'])} != {n}")
        ds[key] = d
    wpd = locate(root, "wpd.json")
    return ds, wpd

def prepare_curve(curve: dict[str,np.ndarray]) -> dict[str,Any]:
    k, f = curve["k_um_inv"].copy(), curve["omega_over_2pi_kHz"].copy()
    order = np.argsort(k, kind="mergesort"); k, f = k[order], f[order]
    kk, ff, merged = [], [], 0
    i = 0
    while i < len(k):
        j = i + 1
        while j < len(k) and abs(k[j]-k[i]) < 1e-9: j += 1
        kk.append(float(np.mean(k[i:j]))); ff.append(float(np.mean(f[i:j])))
        if j-i > 1: merged += j-i-1
        i = j
    k, f = np.asarray(kk), np.asarray(ff)
    noninc = np.flatnonzero(np.diff(f) <= 0)
    return {"k":k, "f":f, "merged":merged, "nonincreasing":noninc}

def bracket(curve, x):
    k = curve["k"]; idx = np.searchsorted(k, x, side="right") - 1
    idx = np.clip(idx, 0, len(k)-2)
    dx = k[idx+1]-k[idx]
    t = (x-k[idx])/dx
    return idx, t, dx

def interpolate(curve, x, method="linear"):
    k,f=curve["k"],curve["f"]
    if np.any(x < k[0]) or np.any(x > k[-1]): raise R0Error("Extrapolation requested")
    if method == "linear": return np.interp(x,k,f)
    if method == "pchip": return PchipInterpolator(k,f,extrapolate=False)(x)
    raise ValueError(method)

def curve_quantities(curve, x):
    idx,t,dx=bracket(curve,x)
    slope=(curve["f"][idx+1]-curve["f"][idx])/dx
    gap=dx
    return slope,t,gap,idx

def loo_rms(curve):
    k,f=curve["k"],curve["f"]
    if len(k)<3: return float("nan"),0
    errs=[]
    for j in range(1,len(k)-1):
        kp=np.delete(k,j); fp=np.delete(f,j)
        pred=np.interp(k[j],kp,fp)
        errs.append(f[j]-pred)
    rms=float(np.sqrt(np.mean(np.square(errs)))) if errs else float("nan")
    return rms,len(errs)

def weights(curve, x):
    slope,t,gap,idx=curve_quantities(curve,x)
    v=(1.0+(RATIO*slope)**2)*(1.0+(1.0-t)**2+t**2)
    return 1.0/v, slope,t,gap,idx

def design_diagnostics(curve, x):
    F=interpolate(curve,x)
    cols=[np.ones_like(F),F,F*x,F*x**2,F*x**3]
    w,*_=weights(curve,x)
    corr=np.eye(len(cols))
    for i in range(len(cols)):
        for j in range(len(cols)):
            a,b=cols[i],cols[j]
            ma=np.sum(w*a)/np.sum(w); mb=np.sum(w*b)/np.sum(w)
            da=a-ma; db=b-mb
            den=np.sqrt(np.sum(w*da*da)*np.sum(w*db*db))
            corr[i,j]=np.sum(w*da*db)/den if den>0 else float("nan")
    cond={}
    for q in SHAPES:
        X=np.column_stack((cols[0],cols[1],cols[q+1]))
        Xs=X*np.sqrt(w)[:,None]
        norms=np.linalg.norm(Xs,axis=0)
        if np.any(norms==0): cond[str(q)]=float("inf")
        else: cond[str(q)]=float(np.linalg.cond(Xs/norms))
    return corr,cond

def fit_models(y, F, x, w):
    # y is residual. Base columns [1,F]; extensions add F*k^q.
    B=np.column_stack((np.ones_like(F),F))
    sw=np.sqrt(w)
    def solve(X):
        Z=X*sw[:,None]; norms=np.linalg.norm(Z,axis=0)
        if np.any(norms==0): return None, float("inf"), None
        Zs=Z/norms
        cond=float(np.linalg.cond(Zs))
        coef_scaled, *_=np.linalg.lstsq(Zs,y*sw,rcond=None)
        coef=coef_scaled/norms
        resid=(y-X@coef)*sw
        rss=float(resid@resid)
        return coef,cond,rss
    cb,condb,rssb=solve(B)
    N=len(y); df2=N-3
    out={}
    if cb is None or df2<=0: raise R0Error("Base model singular or insufficient degrees of freedom")
    for q in SHAPES:
        X=np.column_stack((B,F*x**q))
        ce,cond,rss=solve(X)
        if ce is None or cond>COND_MAX:
            out[q]={"condition":cond,"ill_conditioned":True,"p":1.0,"gamma":float("nan"),"fstat":float("nan")}
            continue
        improvement=max(0.0,rssb-rss)
        fstat=(improvement/1.0)/(rss/df2) if rss>0 else (float("inf") if improvement>0 else 0.0)
        pv=float(f_dist.sf(fstat,1,df2))
        out[q]={"condition":cond,"ill_conditioned":False,"p":pv,"gamma":float(ce[2]),"fstat":float(fstat)}
    if any(v.get("ill_conditioned") for v in out.values()):
        return out, cb, rssb, True
    return out, cb, rssb, False

def holm_decision(fits):
    ps=np.array([fits[q]["p"] for q in SHAPES],float)
    order=sorted(range(3),key=lambda i:(ps[i],SHAPES[i]))
    rejected=set(); still=True
    for rank,i in enumerate(order):
        if still and ps[i] <= ALPHA/(3-rank): rejected.add(SHAPES[i])
        else: still=False
    if not rejected: return {"outcome":"NO_DETECTION","q_star":None,"sign":0,"rejected":[]}
    q=min(rejected,key=lambda z:(fits[z]["p"],z))
    gamma=fits[q]["gamma"]
    sign=0 if gamma==0 else (1 if gamma>0 else -1)
    return {"outcome":"DETECTED","q_star":q,"sign":sign,"rejected":sorted(rejected)}

def clopper_pearson(x,n):
    if n<=0: return [None,None]
    lo=0.0 if x==0 else float(beta.ppf(.025,x,n-x+1))
    hi=1.0 if x==n else float(beta.ppf(.975,x+1,n-x))
    return [lo,hi]

def analyze_one(curve, exp, use_exp_frequency=True, method="linear", weighted=True):
    x=exp["k_um_inv"]
    F=interpolate(curve,x,method)
    w,slope,t,gap,idx=weights(curve,x)
    if not np.all(gap <= .5): raise R0Error("An experimental point has bracketing segment >0.5 um^-1")
    if weighted is False: w=np.ones_like(x)
    if not use_exp_frequency: raise R0Error("Synthetic truth must be supplied by caller")
    residual=exp["omega_over_2pi_kHz"]-F
    fits,base,rss,ill=fit_models(residual,F,x,w)
    decision=holm_decision(fits) if not ill else {"outcome":"NOT_TESTABLE","q_star":None,"sign":0,"rejected":[]}
    sigma=float(np.sqrt(rss/max(1,len(x)-2)))
    return {"N":len(x),"F":F,"residual":residual,"weights":w,"slope":slope,"t":t,"gap":gap,
            "fits":fits,"decision":decision,"ill_conditioned":ill,"sigma_kHz":sigma,
            "sigma_px":sigma/AF,"base_rss":rss}

def audit(root, outdir):
    ds,wpd=load_inputs(root)
    cal=calibration_from_wpd(wpd)
    checks={}
    for key in FILES:
        d=ds[key]
        rec=cal.pixel_to_data(np.column_stack((d["px"],d["py"])))
        dk=np.abs(rec[:,0]-d["k_um_inv"]); df=np.abs(rec[:,1]-d["omega_over_2pi_kHz"])
        # Apply the same WPD reproduction tolerance to all four CSV datasets.
        if np.any(dk > 1e-9+1e-9*np.abs(d["k_um_inv"])) or np.any(df > 1e-9+1e-9*np.abs(d["omega_over_2pi_kHz"])):
            raise R0Error(f"WPD pixel-to-value reproduction failed for {key}")
        checks[key]={"n":len(dk),"max_abs_dk":float(dk.max()),"max_abs_df":float(df.max())}
    result={"spec":SPEC,"stage":"audit","inputs":{k:{"sha256":sha256((root/"realdata/shammass2012/unified_v1/exports"/FILES[k][0])) if (root/"realdata/shammass2012/unified_v1/exports"/FILES[k][0]).exists() else FILES[k][2],"n":FILES[k][1]} for k in FILES},
            "wpd_sha256":sha256(wpd),"calibration_A":cal.A.tolist(),"calibration_c":cal.c.tolist(),"reproduction":checks,
            "note":"No interpolation, residual, physical fit, or synthetic gate executed."}
    outdir.mkdir(parents=True,exist_ok=True); (outdir/"r0_audit.json").write_text(json_text(result,indent=2),encoding="utf-8")
    print(json_text(result,indent=2))

def curve_only(root,outdir):
    ds,wpd=load_inputs(root,curve_only=True); results={}
    for series,(ek,gk) in PAIR.items():
        c=prepare_curve(ds[gk]); x=ds[ek]["k_um_inv"]
        # No experimental frequencies are read in this stage.
        if np.any(x<c["k"][0]) or np.any(x>c["k"][-1]): raise R0Error(f"{series}: experimental k outside curve domain")
        slope,t,gap,idx=curve_quantities(c,x)
        if np.any(gap > .5):
            raise R0Error(f"{series}: an experimental point has bracketing segment >0.5 um^-1")
        rms,nloo=loo_rms(c); corr,conds=design_diagnostics(c,x)
        results[series]={"curve_points_after_merge":len(c["k"]),"merged_duplicate_count":c["merged"],
            "nonincreasing_segments_k_indices":c["nonincreasing"].tolist(),"nonincreasing_count":len(c["nonincreasing"]),
            "experimental_points_in_domain":len(x),"max_bracket_gap":float(gap.max()),
            "slope_khz_per_um_inv":slope.tolist(),"interpolation_t":t.tolist(),
            "loo_rms_khz":rms,"loo_rms_px":rms/AF,"loo_points":nloo,
            "weighted_design_correlations":corr.tolist(),"condition_numbers":conds}
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/"r0_curve_only.json").write_text(json_text({"spec":SPEC,"stage":"curve-only","series":results,
       "note":"Experimental frequencies were not accessed; no experiment-minus-GPE residual computed."},indent=2),encoding="utf-8")
    print(json_text(results,indent=2))

def perturb_and_recalibrate(cal, delta, rng):
    # draw order is fixed: calibration (8), experimental points, curve nodes
    jit=rng.normal(0.0,1.0,size=(4,2))
    # calibration_from_wpd is reconstructed by caller from saved WPD json
    return jit

def simulate_replicate(series_data, cal0, wpd_path, delta, q, m, sign, rng, null=False):
    exp,curve=series_data
    c0=prepare_curve(curve)
    xobs=exp["k_um_inv"]
    F0=interpolate(c0,xobs)
    T=F0*xobs**q
    A=0.0 if null else sign*(m*delta*AF)/float(np.max(np.abs(T)))
    ftrue=F0+A*T
    # Pixel coordinates are generated using archived calibration.
    exp_pix=cal0.data_to_pixel(np.column_stack((xobs,ftrue)))
    curve_pix=cal0.data_to_pixel(np.column_stack((curve["k_um_inv"],curve["omega_over_2pi_kHz"])))
    # Required random draw order: calibration jitter, exp x/y, curve x/y.
    cal_jit=rng.normal(0.0,1.0,size=(4,2))
    exp_jit=rng.normal(0.0,delta,size=exp_pix.shape)
    curve_jit=rng.normal(0.0,delta,size=curve_pix.shape)
    c_jittered=calibration_from_wpd(wpd_path,cal_jit)
    exp_val=c_jittered.pixel_to_data(exp_pix+exp_jit)
    curve_val=c_jittered.pixel_to_data(curve_pix+curve_jit)
    exp_sim={"k_um_inv":exp_val[:,0],"omega_over_2pi_kHz":exp_val[:,1]}
    cur_sim={"k_um_inv":curve_val[:,0],"omega_over_2pi_kHz":curve_val[:,1]}
    cp=prepare_curve(cur_sim); xx=exp_sim["k_um_inv"]
    if np.any(xx<cp["k"][0]) or np.any(xx>cp["k"][-1]): return {"outcome":"NOT_TESTABLE","q_star":None,"sign":0,"ill":True}
    _,_,gaps,_=curve_quantities(cp,xx)
    if np.any(gaps>.5): return {"outcome":"NOT_TESTABLE","q_star":None,"sign":0,"ill":True}
    F=interpolate(cp,xx); w,*_=weights(cp,xx)
    residual=exp_sim["omega_over_2pi_kHz"]-F
    fits,_,_,ill=fit_models(residual,F,xx,w)
    if ill: return {"outcome":"NOT_TESTABLE","q_star":None,"sign":0,"ill":True}
    d=holm_decision(fits)
    d["ill"]=False
    return d

def gate(root,outdir):
    ds,wpd=load_inputs(root); cal=calibration_from_wpd(wpd)
    series_inputs={s:(ds[e],ds[g]) for s,(e,g) in PAIR.items()}
    # Check exported GPE coordinates against WPD conversion before any draws.
    for _,g in PAIR.values():
        d=ds[g]; rec=cal.pixel_to_data(np.column_stack((d["px"],d["py"])))
        if np.any(np.abs(rec[:,0]-d["k_um_inv"])>1e-9+1e-9*np.abs(d["k_um_inv"])) or np.any(np.abs(rec[:,1]-d["omega_over_2pi_kHz"])>1e-9+1e-9*np.abs(d["omega_over_2pi_kHz"])):
            raise R0Error("§8.2 WPD conversion check failed; gate aborted.")
    root_seed=np.random.SeedSequence(SEED)
    children=root_seed.spawn(2*4*(1+3*4*2))
    child_i=0; allout={}
    for series in ("large","small"):
        exp,curve=series_inputs[series]
        c=prepare_curve(curve); x=exp["k_um_inv"]
        corr,conds=design_diagnostics(c,x)
        series_out={"design_correlations":corr.tolist(),"condition_numbers":conds,"deltas":{}}
        for delta in DELTAS:
            delta_out={"null":None,"signals":{}}
            # one child RNG per case; sequential replicates
            rng=np.random.default_rng(children[child_i]); child_i+=1
            null_counts={"detected":0,"not_testable":0}
            for _ in range(N_NULL):
                z=simulate_replicate((exp,curve),cal,wpd,delta,2,2,1,rng,null=True)
                null_counts["detected"] += z["outcome"]=="DETECTED"
                null_counts["not_testable"] += z["outcome"]=="NOT_TESTABLE"
            null_counts["n"]=N_NULL; null_counts["cp95"]=clopper_pearson(null_counts["detected"],N_NULL)
            delta_out["null"]=null_counts
            for q in SHAPES:
                for m in MS:
                    for sign in SIGNS:
                        rng=np.random.default_rng(children[child_i]); child_i+=1
                        counts={"detected":0,"not_testable":0,"q_star":{},"sign":{}}
                        for _ in range(N_SIGNAL):
                            z=simulate_replicate((exp,curve),cal,wpd,delta,q,m,sign,rng)
                            counts["detected"] += z["outcome"]=="DETECTED"
                            counts["not_testable"] += z["outcome"]=="NOT_TESTABLE"
                            if z["outcome"]=="DETECTED":
                                qs=str(z["q_star"]); ss=str(z["sign"])
                                counts["q_star"][qs]=counts["q_star"].get(qs,0)+1
                                counts["sign"][ss]=counts["sign"].get(ss,0)+1
                        counts["n"]=N_SIGNAL; counts["cp95"]=clopper_pearson(counts["detected"],N_SIGNAL)
                        delta_out["signals"][f"q{q}_m{m}_s{sign}"]=counts
            series_out["deltas"][str(delta)]=delta_out
        allout[series]=series_out
    # Derive G1/G2/G3 exactly per spec.
    verdicts={}
    for series,res in allout.items():
        g1=True
        for d in (1,2,4):
            n=res["deltas"][str(d)]["null"]; lo,hi=n["cp95"]
            g1 &= (lo<=.010 and hi<=.015 and n["not_testable"]==0)
        mstar=None; g2=False
        for m in MS:
            ok=True
            for d in (1,2,4):
                for s in SIGNS:
                    z=res["deltas"][str(d)]["signals"][f"q2_m{m}_s{s}"]
                    lo,hi=z["cp95"]
                    if lo is None or lo<.80: ok=False; continue
                    det=z["detected"]
                    q2=z["q_star"].get("2",0); correct=z["sign"].get(str(s),0)
                    if det==0 or clopper_pearson(q2,det)[0]<.80 or clopper_pearson(correct,det)[0]<.98: ok=False
            if ok: mstar=m; g2=True; break
        g3=True; g3_not_estimable=[]
        if mstar is None: g3=False
        else:
            for q in (1,3):
                for d in (1,2,4):
                    for s in SIGNS:
                        z=res["deltas"][str(d)]["signals"][f"q{q}_m{mstar}_s{s}"]
                        det=z["detected"]
                        if det==0:
                            g3=False; g3_not_estimable.append(f"q{q},d{d},s{s}")
                        else:
                            wrong=z["q_star"].get("2",0)
                            if clopper_pearson(wrong,det)[1]>.20: g3=False
        if not g1 or not g2: verdict="NOT_TESTABLE"
        elif g3: verdict="TESTABLE"
        else: verdict="SIGN_ONLY"
        verdicts[series]={"G1":bool(g1),"G2":bool(g2),"G3":bool(g3),"G3_not_estimable_cells":g3_not_estimable,"m_star":mstar,"verdict":verdict}
    outdir.mkdir(parents=True,exist_ok=True)
    payload={"spec":SPEC,"stage":"synthetic-gate","seed":SEED,"replicates_per_series_delta":{"null":N_NULL,"signal":N_SIGNAL},
             "verdicts":verdicts,"results":allout,"warning":"Gate results only; real residual not opened."}
    (outdir/"r0_gate_results.json").write_text(json_text(payload,indent=2),encoding="utf-8")
    print(json_text(verdicts,indent=2))

def real_analysis(root,outdir,gate_file):
    gate_data=json.loads(gate_file.read_text(encoding="utf-8"))
    ds,wpd=load_inputs(root)
    out={}
    for series,(ek,gk) in PAIR.items():
        gv=gate_data["verdicts"][series]["verdict"]
        if gv=="NOT_TESTABLE":
            out[series]={"gate_verdict":gv,"real_analysis":"NOT_RUN"}
            continue
        primary=analyze_one(prepare_curve(ds[gk]),ds[ek],True,"linear",True)
        pchip=analyze_one(prepare_curve(ds[gk]),ds[ek],True,"pchip",True)
        ols=analyze_one(prepare_curve(ds[gk]),ds[ek],True,"linear",False)
        def compact(z):
            return {"decision":z["decision"],"fits":z["fits"],"sigma_kHz":z["sigma_kHz"],"sigma_px":z["sigma_px"],"N":z["N"]}
        labels=[]
        if (primary["decision"] != pchip["decision"]): labels.append("interpolation-sensitive")
        if (primary["decision"] != ols["decision"]): labels.append("weight-sensitive")
        if primary["sigma_px"]>8: labels.append("outside_gate_range")
        out[series]={"gate_verdict":gv,"primary":compact(primary),"pchip_sensitivity":compact(pchip),
                     "ols_sensitivity":compact(ols),"labels":labels}
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/"r0_real_analysis.json").write_text(json_text({"spec":SPEC,"stage":"real-analysis","series":out},indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("root",type=Path,help="R0_v0.5.2 folder containing realdata/shammass2012, or shammass2012 folder")
    ap.add_argument("--stage",choices=("audit","curve-only","gate","real"),default="audit")
    ap.add_argument("--out",type=Path,default=Path("R0_RESULTS"))
    ap.add_argument("--authorize-gate",action="store_true",help="Explicitly authorize the single synthetic gate run")
    ap.add_argument("--authorize-real",action="store_true",help="Explicitly authorize opening experiment-minus-GPE residual")
    ap.add_argument("--gate-file",type=Path,help="Gate result JSON required for real analysis")
    args=ap.parse_args()
    try:
        if args.stage=="audit": audit(args.root.resolve(),args.out)
        elif args.stage=="curve-only": curve_only(args.root.resolve(),args.out)
        elif args.stage=="gate":
            if not args.authorize_gate: raise R0Error("Blocked: synthetic gate requires --authorize-gate and separate explicit permission.")
            gate(args.root.resolve(),args.out)
        elif args.stage=="real":
            if not args.authorize_real: raise R0Error("Blocked: real residual requires --authorize-real and separate explicit permission.")
            if not args.gate_file or not args.gate_file.is_file(): raise R0Error("Real analysis requires a valid gate result file.")
            real_analysis(args.root.resolve(),args.out,args.gate_file)
        return 0
    except Exception as e:
        print(f"R0 ERROR: {e}",file=sys.stderr); return 2

if __name__=="__main__": raise SystemExit(main())
