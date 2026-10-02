#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
residual_check.py — Module 1 plumbing check (no f^p templates, no Λ estimate).

For one event and its model labels: load the GWOSC strain, take the PE analysis
segment (duration from the PE configuration, merger 2 s before the end), use the
released PSDs, generate the reference waveform (maximum-likelihood sample, the
label's own settings), project it onto each detector (antenna response and
arrival time), refine a common time shift (default ±100 ms) and a common complex amplitude
(phase and distance), and form the residual r = d − h.

Check: the reduced χ² of the whitened residual in the analysis band should match
that of off-source segments (pure noise), while the data before subtraction must
lie above it. The PSD is scaled by the window power (as in bilby). Note: with a
least-squares complex amplitude the χ² drop equals SNR_mf² identically, so it is not
used as a check; SNR_mf is compared with the catalog SNR instead.

  python residual_check.py --event GW150914
  python residual_check.py --event all --out residual_check.csv
"""
import argparse
import csv
import os
import sys

import h5py
import numpy as np
from scipy.signal.windows import tukey

import lal
import lalsimulation as ls

import f3_derivatives as f3

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

POST_TRIGGER = 2.0          # s after merger (bilby default)
ROLL_OFF = 0.4              # s, Tukey roll-off (bilby default)
N_OFF = 8                   # off-source segments
OFF_GAP = 64.0              # s between consecutive off-source segments


def read_strain(path):
    with h5py.File(path, "r") as fh:
        x = fh["strain/Strain"][()]
        dt = fh["strain/Strain"].attrs["Xspacing"]
        t0 = fh["meta/GPSstart"][()]
    return np.asarray(x, float), float(dt), float(t0)


def segment(x, dt, t0, start, duration):
    i0 = int(round((start - t0) / dt))
    n = int(round(duration / dt))
    if i0 < 0 or i0 + n > len(x):
        raise ValueError("segment outside the file")
    seg = x[i0:i0 + n]
    if not np.all(np.isfinite(seg)):
        raise ValueError("NaN in segment")
    return seg


def polarizations(label, sample, st):
    """h+, hx on the rfft grid of the analysis segment (merger at t = 0)."""
    approx = ("IMRPhenomXPHM" if "XPHM" in label else
              "SEOBNRv4PHM" if "v4PHM" in label else "SEOBNRv5PHM")
    p = {k: sample[k] for k in f3.STEPS}
    m1, m2, s, iota = f3.physical(p, st)
    dist = p["luminosity_distance"] * 1e6 * lal.PC_SI
    srate, dur = st["srate"], st["duration"]
    N = int(round(srate * dur))
    f = np.fft.rfftfreq(N, 1 / srate)
    unknown = set()
    d = f3.laldict(st["wf_args"], unknown)
    if approx == "IMRPhenomXPHM":
        hp, hc = ls.SimInspiralChooseFDWaveform(
            m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *s, dist, iota, p["phase"], 0, 0, 0,
            1.0 / dur, st["f_wf"], srate / 2, st["f_ref"], d, ls.IMRPhenomXPHM)
        out = []
        for h in (hp, hc):
            arr = np.zeros(len(f), complex)
            n = min(h.data.length, len(f))
            arr[:n] = h.data.data[:n]
            out.append(arr)
        return f, out[0], out[1]
    def gen_td(dt):
        if approx == "SEOBNRv4PHM":
            return ls.SimInspiralChooseTDWaveform(
                m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *s, dist, iota, p["phase"], 0, 0, 0,
                dt, st["f_wf"], st["f_ref"], d, ls.SEOBNRv4PHM)
        from pyseobnr.generate_waveform import GenerateWaveform
        gp = {"mass1": m1, "mass2": m2, "spin1x": s[0], "spin1y": s[1], "spin1z": s[2],
              "spin2x": s[3], "spin2y": s[4], "spin2z": s[5], "deltaT": dt,
              "f22_start": st["f_wf"], "f_ref": st["f_ref"], "phi_ref": p["phase"],
              "distance": p["luminosity_distance"], "inclination": iota,
              "approximant": "SEOBNRv5PHM"}
        gp.update(st["wf_args"])
        return GenerateWaveform(gp).generate_td_polarizations()

    # The EOB model may refuse the PE sampling rate (ringdown above Nyquist): generate
    # at 2× or 4× the rate; only frequencies up to f_high are used, so the result on the
    # PE frequency grid is unchanged.
    for factor in (1, 2, 4):
        dt = 1.0 / (srate * factor)
        try:
            hp, hc = gen_td(dt)
            break
        except RuntimeError as e:
            if factor == 4 or "domain" not in str(e).lower():
                raise
    Nk = N * factor
    fk = np.fft.rfftfreq(Nk, dt)
    out = []
    for h in (hp, hc):
        x = np.asarray(h.data.data, float)
        t = float(h.epoch) + dt * np.arange(len(x))
        # Only the part inside the analysis segment matters (merger at t = 0, segment
        # [−(dur − POST_TRIGGER), +POST_TRIGGER]); crop and taper the start like the data.
        keep = t >= -(dur - POST_TRIGGER)
        x, t = x[keep], t[keep]
        nt = min(len(x), int(round(ROLL_OFF / dt)))
        if nt > 1:
            x[:nt] *= 0.5 * (1 - np.cos(np.pi * np.arange(nt) / nt))
        if len(x) > Nk:
            raise ValueError(f"TD waveform ({len(x)}) longer than the segment ({Nk}) after cropping")
        X = np.fft.rfft(x, n=Nk) * dt * np.exp(-2j * np.pi * fk * t[0])
        out.append(X[:len(f)])
    return f, out[0], out[1]


def project(f, hp, hc, ifo, sample, seg_start):
    det = ls.DetectorPrefixToLALDetector(ifo)
    tc = sample["geocent_time"]
    gmst = lal.GreenwichMeanSiderealTime(tc)
    Fp, Fc = lal.ComputeDetAMResponse(det.response, sample["ra"], sample["dec"], sample["psi"], gmst)
    dt_det = lal.TimeDelayFromEarthCenter(det.location, sample["ra"], sample["dec"], tc)
    return (Fp * hp + Fc * hc) * np.exp(-2j * np.pi * f * (tc + dt_det - seg_start))


def check(event, args, verbose=True):
    f5row = next(r for r in csv.DictReader(open(args.f5, encoding="utf-8")) if r["commonName"] == event)
    manifest = {(r["event"], r["ifo"]): r for r in
                csv.DictReader(open(os.path.join(args.strain_dir, "strain_manifest.csv"), encoding="utf-8"))}
    pe_file = f5row["file"].split(" ")[0]
    labels = [l for l in (f5row["xphm_label"], f5row["eob_label"]) if l and f5row.get(
        "xphm_ok" if l == f5row["xphm_label"] else "eob_ok") == "True"]
    ifos = [i for i in f5row["ifos"].split(",") if i]
    rows = []
    with h5py.File(os.path.join(args.pe_dir, pe_file), "r") as fh:
        for label in labels:
            try:
                sample, cfg, _ = f3.read_label(fh, label)
                st = f3.settings(cfg)
                psds = {i: np.asarray(fh[label]["psds"][i][()]) for i in ifos if i in fh[label]["psds"]}
                dur, srate = st["duration"], st["srate"]
                seg_start = sample["geocent_time"] + POST_TRIGGER - dur
                f, hp, hc = polarizations(label, sample, st)
                band = (f >= st["f_an"]) & (f <= st["f_high"])
                df = f[1] - f[0]
                win = tukey(int(round(dur * srate)), alpha=2 * ROLL_OFF / dur)
                w2 = np.mean(win ** 2)
                data, model, Seff, offs = {}, {}, {}, {}
                for ifo in ifos:
                    x, dt, t0 = read_strain(os.path.join(args.strain_dir, manifest[(event, ifo)]["file"]))
                    if abs(1 / dt - srate) > 1e-6:
                        from scipy.signal import resample_poly
                        x = resample_poly(x, int(srate), int(round(1 / dt)))
                        dt = 1.0 / srate
                    data[ifo] = np.fft.rfft(segment(x, dt, t0, seg_start, dur) * win) * dt
                    model[ifo] = project(f, hp, hc, ifo, sample, seg_start)
                    Seff[ifo] = np.interp(f, psds[ifo][:, 0], psds[ifo][:, 1], left=np.inf, right=np.inf) * w2
                    offs[ifo] = []
                    for k in range(1, 4 * N_OFF):
                        for side in (-1, +1):
                            if len(offs[ifo]) >= N_OFF:
                                break
                            s0 = seg_start + side * k * (dur + OFF_GAP)
                            try:
                                offs[ifo].append(np.fft.rfft(segment(x, dt, t0, s0, dur) * win) * dt)
                            except ValueError:
                                pass
                    if len(offs[ifo]) < 4:
                        raise ValueError(f"only {len(offs[ifo])} valid off-source segments in {ifo}")

                def ip(a, b, ifo):
                    return 4 * df * np.sum((np.conj(a) * b / Seff[ifo])[band])
                hh = sum(ip(model[i], model[i], i).real for i in ifos)
                best = (0.0, 0j)
                for tshift in np.arange(-args.tmax, args.tmax + 1e-9, 5e-5):
                    ph = np.exp(-2j * np.pi * f * tshift)
                    z = sum(ip(model[i] * ph, data[i], i) for i in ifos)
                    if abs(z) > abs(best[1]):
                        best = (tshift, z)
                tshift, z = best
                a = z / hh
                ph = np.exp(-2j * np.pi * f * tshift)
                snr_mf = abs(z) / np.sqrt(hh)

                def chi2(vec, ifo):
                    return float((4 * df * np.abs(vec[band]) ** 2 / Seff[ifo][band]).sum())
                zres = []
                for i in ifos:
                    nb = 2 * band.sum()
                    off = [chi2(o, i) / nb for o in offs[i]]
                    zres.append(((chi2(data[i] - a * model[i] * ph, i) / nb) - np.mean(off)) / np.std(off, ddof=1))
                edge = abs(abs(tshift) - args.tmax) < 1e-4
                row = {"event": event, "label": label, "snr_opt": np.sqrt(hh), "snr_mf": snr_mf,
                       "dt_ms": tshift * 1e3, "abs_a": abs(a), "arg_a": float(np.angle(a)),
                       "snr_ratio": snr_mf / args.catalog_snr.get(event, float("nan")),
                       "res_z_max": max(abs(v) for v in zres), "n_off": min(len(offs[i]) for i in ifos),
                       "edge": edge, "error": ""}
            except Exception as e:
                row = {"event": event, "label": label, "error": f"{type(e).__name__}: {str(e)[:120]}"}
            rows.append(row)
            if verbose:
                if row.get("error"):
                    print(f"  {event:18s} {label:30s} ERROR {row['error']}")
                else:
                    print(f"  {event:18s} {label:30s} SNR_mf {row['snr_mf']:6.2f} (opt {row['snr_opt']:6.2f})  "
                          f"Δt {row['dt_ms']:+7.2f} ms{' EDGE' if row['edge'] else '     '}  |a| {row['abs_a']:.3f}  "
                          f"arg {row['arg_a']:+.2f}  SNR_mf/catalog {row['snr_ratio']:.2f}  max|res z| {row['res_z_max']:.2f}"
                          f"  (off-source {row['n_off']})",
                          flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", default="GW150914", help="event name, or 'all' for every F5-eligible event")
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--strain-dir", default=os.path.expanduser("~/gwdata/strain"))
    ap.add_argument("--f5", default="f5/f5_results.csv")
    ap.add_argument("--tmax", type=float, default=0.1,
                    help="half-width of the time-refinement window [s] (default 100 ms)")
    ap.add_argument("--events-list", default="event_selection/event_list_v1.csv")
    ap.add_argument("--out", default="residual_check.csv")
    args = ap.parse_args()
    args.catalog_snr = {r["commonName"]: float(r["network_matched_filter_snr"] or "nan")
                        for r in csv.DictReader(open(args.events_list, encoding="utf-8"))
                        if r["status"] == "SELECTED"}
    if args.event == "all":
        events = [r["commonName"] for r in csv.DictReader(open(args.f5, encoding="utf-8"))
                  if r["status"] in ("BOTH_OK", "SINGLE_MODEL")]
    else:
        events = [args.event]
    rows = []
    for ev in events:
        rows += check(ev, args)
    keys = ["event", "label", "snr_opt", "snr_mf", "snr_ratio", "dt_ms", "abs_a", "arg_a",
            "res_z_max", "n_off", "edge", "error"]
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows([{k: r.get(k, "") for k in keys} for r in rows])
    ok = [r for r in rows if not r.get("error")]
    print(f"\n{len(rows)} label checks, {len(rows) - len(ok)} errors; "
          f"at window edge: {sum(r['edge'] for r in ok)}; "
          f"SNR_mf/catalog outside 0.8–1.2: {sum(not (0.8 <= r['snr_ratio'] <= 1.2) for r in ok)}; "
          f"written {args.out}")


if __name__ == "__main__":
    main()
