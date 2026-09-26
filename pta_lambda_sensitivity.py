"""
Order-of-magnitude sensitivity of PTAs to the Lambda_GW dispersion of Paper 3.

  ΔΨ(f) = -4π³ Λ_GW K(z) f³ / c³ ,   [Λ_GW] = m³/s, [K] = s  (K ≈ D/c assumed)
  equivalent to ω² = c²k² + α k⁴ with α = Λ_GW·c  →  |v_ph/c - 1| ≈ 2π² |Λ_GW| f² / c³

Run: python pta_lambda_sensitivity.py
"""
import numpy as np

c = 299792458.0
LAMBDA = 7.24e-3                  # |Λ_GW| upper edge of the published range [m^3/s]
GPC, MPC, KPC = 3.0857e25, 3.0857e22, 3.0857e19


def dpsi(f, D):                    # |ΔΨ| in rad for a source at distance D [m]
    return 4 * np.pi**3 * LAMBDA * (D / c) * f**3 / c**3


def dv(f):                         # |v_ph/c - 1|
    return 2 * np.pi**2 * LAMBDA * f**2 / c**3


if __name__ == "__main__":
    print("PTA band, source at 1 Gpc, pulsar at 1 kpc")
    print(f"{'f [nHz]':>8} {'|ΔΨ| [rad]':>12} {'|v/c-1|':>10} {'HD effect':>11}")
    for f in (1e-9, 3e-9, 1e-8, 3e-8):
        print(f"{f*1e9:8.0f} {dpsi(f, GPC):12.1e} {dv(f):10.1e} {dv(f)*2*np.pi*f*KPC/c:11.1e}")
    print("\nOrder-of-magnitude reference points (not detectability claims):")
    for f, D, lab in ((1e-2, GPC, "LISA band, 10 mHz, 1 Gpc"),
                      (100, 40 * MPC, "100 Hz, 40 Mpc"),
                      (1000, GPC, "1 kHz, 1 Gpc")):
        print(f"  {lab:<28s} |ΔΨ| = {dpsi(f, D):.1e} rad")
