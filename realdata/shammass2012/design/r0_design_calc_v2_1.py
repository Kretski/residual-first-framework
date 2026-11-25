"""Design calculation v2.1 (post-gate, design only). v2.1: exact 95% Clopper-Pearson
intervals on every simulated rate; SEM-lognormal sensitivity as an ensemble (new profile per
replicate) next to the single fixed profile, which is labelled conditional.
 Reads ONLY GPE curves and the k column
of the experimental CSVs; no experimental frequency is read, no residual is formed.
Fixes v1: full partial-F with candidate-dependent denominator + Holm (R0 rule);
separates R1 (single pre-declared F*k^2 test, alpha=0.01); adds misspecification
sensitivity (correlated errors, unequal SEM) for the raw-number scenario."""
import csv, json, sys
import numpy as np
from scipy import stats, optimize
ROOT, OUT = sys.argv[1], sys.argv[2]
AK, AF = 0.0065830, 0.0032438
SEED = 20261125

def curve(n):
    r = list(csv.DictReader(open(f"{ROOT}/gpe_{n}.csv")))
    k = np.array([float(x["k_um_inv"]) for x in r]); f = np.array([float(x["omega_over_2pi_kHz"]) for x in r])
    o = np.argsort(k); return k[o], f[o]
def exp_k(n): return np.array([float(x["k_um_inv"]) for x in csv.DictReader(open(f"{ROOT}/exp_{n}.csv"))])
def gate_w(kc, fc, kx):
    j = np.clip(np.searchsorted(kc, kx) - 1, 0, len(kc) - 2)
    s = (fc[j+1]-fc[j])/(kc[j+1]-kc[j]); t = (kx-kc[j])/(kc[j+1]-kc[j]); c = (AK/AF*s)**2
    return 1/((1+c)*(1+(1-t)**2+t**2))

class Design:
    def __init__(s, F, kx, w):
        s.n = len(kx); s.sw = np.sqrt(w)
        B = np.column_stack([np.ones(s.n), F]) * s.sw[:, None]
        Q, _ = np.linalg.qr(B); s.P = np.eye(s.n) - Q @ Q.T
        s.T = {q: F * kx**q for q in (1, 2, 3)}
        s.E = {}
        for q in (1, 2, 3):
            u = s.P @ (s.T[q] * s.sw); s.E[q] = u / np.linalg.norm(u)
    def amp(s, q, m):  # max|A T_q| = m (noise unit = 1 on the frequency scale of a weight-1 point)
        return m / np.max(np.abs(s.T[q]))
    def lam(s, q_true, q_test, m):
        y = s.P @ (s.amp(q_true, m) * s.T[q_true] * s.sw); return float((s.E[q_test] @ y)**2)

def fstats(D, Y):
    """Y: (reps, n) whitened data. Returns partial F for q=1,2,3 (candidate-dependent RSS_q)."""
    Yr = Y @ D.P.T; rss0 = np.sum(Yr**2, axis=1); out = []
    for q in (1, 2, 3):
        pr = Yr @ D.E[q]; rssq = rss0 - pr**2
        out.append(pr**2 / (rssq / (D.n - 3)))
    return np.column_stack(out)

def r0_rule(D, Y, alpha=0.01):
    """Holm over p1..p3; q* = smallest raw p among rejected (ties->smaller q). 0 = no detection."""
    p = stats.f.sf(fstats(D, Y), 1, D.n - 3)
    o = np.argsort(p, axis=1, kind="stable"); ps = np.take_along_axis(p, o, 1)
    rej_sorted = np.zeros_like(ps, bool); ok = np.ones(len(p), bool)
    for i, a in enumerate((alpha/3, alpha/2, alpha)):
        ok &= ps[:, i] <= a; rej_sorted[:, i] = ok
    det = rej_sorted[:, 0]  # Holm: any rejection iff first step rejects
    qstar = np.where(det, o[:, 0] + 1, 0)  # smallest raw p among rejected = first in order
    return qstar

def r1_rule(D, Y, alpha=0.01):
    F = fstats(D, Y)[:, 1]; return stats.f.sf(F, 1, D.n - 3) <= alpha

def cp(x, n):
    lo = 0.0 if x == 0 else float(stats.beta.ppf(0.025, x, n - x + 1))
    hi = 1.0 if x == n else float(stats.beta.ppf(0.975, x + 1, n - x))
    return [lo, hi]

def rate(x, n): return {"rate": x / n, "x": int(x), "n": int(n), "cp95": cp(int(x), int(n))}

def noise(rng, reps, n, rho=0.0, sem_shape=None, sem_lognorm_ensemble=None):
    z = rng.standard_normal((reps, n))
    if rho > 0:  # AR(1) along the k-ordered points
        for i in range(1, n): z[:, i] = rho * z[:, i-1] + np.sqrt(1 - rho**2) * z[:, i]
    if sem_shape is not None: z = z * sem_shape
    if sem_lognorm_ensemble is not None: z = z * np.exp(sem_lognorm_ensemble * rng.standard_normal((reps, n)))
    return z

def confusion(D, rng, m_list, reps=20000):
    res = {}
    for m in m_list:
        row = {}
        for qt in (1, 2, 3):
            sig = D.amp(qt, m) * D.T[qt] * D.sw
            for sg in (-1, 1):
                qs = r0_rule(D, sg * sig + noise(rng, reps, D.n))
                c = np.bincount(qs, minlength=4)
                row[f"q{qt}_s{sg:+d}"] = {k: rate(c[i], reps) for i, k in enumerate(("none", "q1", "q2", "q3"))}
        res[f"m{m:g}"] = row
    return res

def r1_power_curve(D, rng, m_list, reps=20000, **nk):
    out = {}
    for m in m_list:
        e = {}
        for qt in (1, 2, 3):
            sig = D.amp(qt, m) * D.T[qt] * D.sw
            hits = np.concatenate([r1_rule(D, sg*sig + noise(rng, reps, D.n, **nk)) for sg in (-1, 1)])
            e[f"q{qt}"] = rate(hits.sum(), hits.size)
        out[f"m{m:g}"] = e
    return out

def m80_analytic(D, alpha, which):
    d2 = D.n - 3; c = stats.f.ppf(1 - alpha, 1, d2)
    return float(optimize.brentq(lambda m: stats.ncf.sf(c, 1, d2, D.lam(2, 2, m)) - 0.8, 0.1, 1e4))

rng = np.random.default_rng(SEED); out = {"seed": SEED}
for name in ("large", "small"):
    kc, fc = curve(name); kx = np.sort(exp_k(name)); F = np.interp(kx, kc, fc)
    sc = {"gate_like": Design(F, kx, gate_w(kc, fc, kx)), "raw_N_obs": Design(F, kx, np.ones(len(kx)))}
    kd = np.linspace(kx.min(), kx.max(), 32); sc["raw_N32_idealized"] = Design(np.interp(kd, kc, fc), kd, np.ones(32))
    r = {}
    for key, D in sc.items():
        r[key] = {"N": D.n,
                  "R0_m80_first_step_analytic": m80_analytic(D, 0.01/3, key),
                  "R1_m80_analytic_alpha001": m80_analytic(D, 0.01, key),
                  "R0_confusion_full_partialF_Holm": confusion(D, rng, (8, 16, 32, 64, 128)),
                  "R1_power_sim_independent": r1_power_curve(D, rng, (0, 8, 16, 32))}
    # misspecification (analysis assumes independent, equal SEM) for raw scenarios
    for key in ("raw_N_obs", "raw_N32_idealized"):
        D = sc[key]; ms = {}
        for rho in (0.3, 0.6):
            ms[f"AR1_rho{rho}"] = r1_power_curve(D, rng, (0, 16), rho=rho)
        kk = kx if key == "raw_N_obs" else kd
        ms["SEM_rising_1to3"] = r1_power_curve(D, rng, (0, 16), sem_shape=1 + 2*(kk-kk.min())/(kk.max()-kk.min()))
        ms["SEM_lognorm0.5_single_fixed_profile_CONDITIONAL"] = r1_power_curve(D, rng, (0, 16), sem_shape=np.exp(0.5*np.random.default_rng(SEED+1).standard_normal(D.n)))
        ms["SEM_lognorm0.5_ensemble_new_profile_per_replicate"] = r1_power_curve(D, rng, (0, 16), sem_lognorm_ensemble=0.5)
        r[key]["R1_misspecification"] = ms
    out[name] = r
json.dump(out, open(OUT, "w"), indent=1, default=float)
