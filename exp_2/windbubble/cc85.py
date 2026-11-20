"""Chevalier & Clegg (1985) steady adiabatic wind from a uniform source region.

Mass and energy are deposited uniformly (per unit volume) inside a sphere of
radius R at total rates Mdot and Edot.  The steady solution is transonic, with
M = 1 exactly at r = R.  Bernoulli: u^2/2 + c^2/(gamma-1) = Edot/Mdot.
"""
import numpy as np
from scipy.optimize import brentq

from .units import GAMMA


def _lhs_inner(M, g=GAMMA):
    return (((3 * g + 1.0 / M**2) / (1 + 3 * g)) ** (-(3 * g + 1) / (5 * g + 1))
            * ((g - 1 + 2.0 / M**2) / (1 + g)) ** ((g + 1) / (2 * (5 * g + 1))))


def _lhs_outer(M, g=GAMMA):
    return (M ** (2 / (g - 1))
            * ((g - 1 + 2.0 / M**2) / (1 + g)) ** ((g + 1) / (2 * (g - 1))))


def mach(r, R, g=GAMMA):
    """Mach number of the CC85 solution at radius r (array OK)."""
    r = np.atleast_1d(np.asarray(r, dtype=float))
    out = np.empty_like(r)
    for i, ri in enumerate(r):
        x = ri / R
        if x <= 0:
            out[i] = 0.0
        elif x < 1.0:
            f = lambda M: _lhs_inner(M, g) - x
            out[i] = brentq(f, 1e-12, 1.0, xtol=1e-15, rtol=1e-13)
        elif x == 1.0:
            out[i] = 1.0
        else:
            f = lambda M: np.log(_lhs_outer(M, g)) - 2 * np.log(x)
            out[i] = brentq(f, 1.0, 1e8, xtol=1e-14, rtol=1e-13)
    return out


def wind_state(r, Mdot, Edot, R, g=GAMMA):
    """Return (rho, u, P) of the CC85 wind at radius r [cgs]."""
    r = np.atleast_1d(np.asarray(r, dtype=float))
    M = mach(r, R, g)
    vinf2 = 2.0 * Edot / Mdot
    # u^2/2 + u^2/(M^2 (g-1)) = vinf^2/2
    u = np.sqrt(vinf2 / (1.0 + 2.0 / ((g - 1) * np.maximum(M, 1e-300) ** 2)))
    u = np.where(M > 0, u, 0.0)
    mdot_enc = Mdot * np.minimum(r / R, 1.0) ** 3
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = np.where(r > 0, mdot_enc / (4 * np.pi * r**2 * u), 0.0)
        c2 = np.where(M > 0, (u / M) ** 2, 0.0)
    # centre: limit r->0 of rho = Mdot/(4 pi R^3) * r/u, u ~ r -> finite
    if np.any(r <= 0):
        rs = 1e-6 * R
        rho0, _, P0 = wind_state(rs, Mdot, Edot, R, g)
        rho = np.where(r > 0, rho, rho0[0])
        c2 = np.where(r > 0, c2, g * P0[0] / rho0[0])
    P = rho * c2 / g
    return rho, u, P
