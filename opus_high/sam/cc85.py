"""Chevalier & Clegg (1985) steady adiabatic wind with uniform mass and
energy injection inside a sphere of radius R_inj.

All inputs/outputs in cgs.
"""
import numpy as np
from scipy.optimize import brentq

from .units import GAMMA


def _mach_inside(x, g=GAMMA):
    """Mach number for x = r/R < 1 (subsonic branch)."""
    def f(M):
        a = ((3 * g + 1 / M**2) / (1 + 3 * g)) ** (-(3 * g + 1) / (5 * g + 1))
        b = ((g - 1 + 2 / M**2) / (1 + g)) ** ((g + 1) / (2 * (5 * g + 1)))
        return np.log(a * b) - np.log(x)
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    return brentq(f, 1e-12, 1.0, xtol=1e-14)


def _mach_outside(x, g=GAMMA):
    """Mach number for x = r/R > 1 (supersonic branch)."""
    def f(M):
        return (2 / (g - 1)) * np.log(M) + (g + 1) / (2 * (g - 1)) * np.log(
            (g - 1 + 2 / M**2) / (1 + g)) - 2 * np.log(x)
    if x <= 1:
        return 1.0
    return brentq(f, 1.0, 1e8, xtol=1e-14)


def mach(x, g=GAMMA):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    out = np.empty_like(x)
    for i, xi in enumerate(x):
        out[i] = _mach_inside(xi, g) if xi < 1 else _mach_outside(xi, g)
    return out


def wind_profile(r, Mdot, L, R_inj, g=GAMMA):
    """Return rho, v, P (cgs) of the CC85 solution at radii r [cm]."""
    r = np.atleast_1d(np.asarray(r, dtype=float))
    r = np.maximum(r, 1e-8 * R_inj)
    vinf = np.sqrt(2 * L / Mdot)
    x = r / R_inj
    M = mach(x, g)
    v = vinf * M / np.sqrt(M**2 + 2 / (g - 1))
    q = Mdot / (4.0 / 3.0 * np.pi * R_inj**3)
    rho = np.where(x < 1, q * r / (3 * v), Mdot / (4 * np.pi * r**2 * v))
    c = v / M
    P = rho * c**2 / g
    return rho, v, P
