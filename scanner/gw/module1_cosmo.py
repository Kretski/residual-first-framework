#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module1_cosmo.py — cosmology for Module 1 (MODULE1_DESIGN §2).

Planck15 as used by LVK for the dispersion tests: H0 = 67.74 km/s/Mpc,
Omega_m = 0.3075, flat, no radiation (Omega_Lambda = 1 - Omega_m).

  I4(z)  = (c/H0) * int_0^z (1+z')^2 / E(z') dz'          [m]
  d_L(z) = (1+z) (c/H0) * int_0^z 1 / E(z') dz'             [m]
  z_from_dl(d_L) inverts d_L(z) (bisection).

Template prefactor (group velocity, §2):  K(z) = -4 pi^3 I4(z) / c^3,
T_p(f) = i * K * f^p * h(f), so that r = Lambda * T_3 with Lambda in m^2.

Self-test (python module1_cosmo.py): reproduces the F4 value
I4(z = 0.085) = 1.23632e25 m (f4_coefficient.py, Planck15).
"""
import numpy as np

C_SI = 299792458.0
MPC_SI = 3.085677581491367e22
H0 = 67.74 * 1e3 / MPC_SI          # 1/s
OM = 0.3075
OL = 1.0 - OM


def E(z):
    return np.sqrt(OM * (1 + z) ** 3 + OL)


def _integral(fun, z, n=4001):
    zz = np.linspace(0.0, z, n)
    return np.trapezoid(fun(zz), zz) if hasattr(np, "trapezoid") else np.trapz(fun(zz), zz)


def I4(z):
    return (C_SI / H0) * _integral(lambda x: (1 + x) ** 2 / E(x), z)


def d_L(z):
    return (1 + z) * (C_SI / H0) * _integral(lambda x: 1 / E(x), z)


def z_from_dl(dl_m, zmax=10.0):
    lo, hi = 0.0, zmax
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if d_L(mid) < dl_m:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def K(z):
    """Prefactor of T_p: T_p = i * K(z) * f^p * h."""
    return -4 * np.pi ** 3 * I4(z) / C_SI ** 3


if __name__ == "__main__":
    ref = 1.23632e25
    got = I4(0.085)
    print(f"I4(0.085) = {got:.5e} m   (F4 / f4_coefficient.py, Planck15: {ref:.5e}; "
          f"ratio {got / ref:.5f})")
    z = z_from_dl(400.0 * MPC_SI)
    print(f"z(d_L = 400 Mpc) = {z:.5f}   (F4 used z = 0.08500 for D_L = 400 Mpc)")
    print(f"K(0.085) = {K(0.085):.4e} s^3/m^2 ... units: m / (m/s)^3 = s^3/m^2")
    lam = 7.4e-12
    for f in (20.0, 100.0, 300.0):
        print(f"  dPsi at Lambda = {lam:.1e} m^2, f = {f:5.0f} Hz: {lam * K(0.085) * f ** 3:+.3e} rad")
