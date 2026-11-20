"""Validation tests of the hydro code: Sod shock tube and CC85 wind."""
import numpy as np
from scipy.optimize import brentq

from hydro1d import Hydro1D
from hydro1d.solver import uniform_injection_source
from windbubble.cc85 import wind_state


def sod_exact(x, t, x0=0.5, g=1.4, WL=(1, 0, 1), WR=(0.125, 0, 0.1)):
    rL, uL, PL = WL
    rR, uR, PR = WR
    cL, cR = np.sqrt(g * PL / rL), np.sqrt(g * PR / rR)

    def f(P, r, Pk, c):
        if P > Pk:
            A, B = 2 / ((g + 1) * r), (g - 1) / (g + 1) * Pk
            return (P - Pk) * np.sqrt(A / (P + B))
        return 2 * c / (g - 1) * ((P / Pk) ** ((g - 1) / (2 * g)) - 1)
    Ps = brentq(lambda P: f(P, rL, PL, cL) + f(P, rR, PR, cR) + uR - uL, 1e-8, 10)
    us = 0.5 * (uL + uR) + 0.5 * (f(Ps, rR, PR, cR) - f(Ps, rL, PL, cL))
    rsL = rL * (Ps / PL) ** (1 / g)
    rsR = rR * ((Ps / PR + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * Ps / PR + 1))
    S = uR + cR * np.sqrt((g + 1) / (2 * g) * Ps / PR + (g - 1) / (2 * g))
    csL = cL * (Ps / PL) ** ((g - 1) / (2 * g))
    xi = (x - x0) / t
    rho = np.where(xi < uL - cL, rL, 0.0)
    u = np.zeros_like(x); P = np.zeros_like(x)
    fan = (xi >= uL - cL) & (xi < us - csL)
    rho = np.where(fan, rL * (2 / (g + 1) + (g - 1) / ((g + 1) * cL) * (uL - xi)) ** (2 / (g - 1)), rho)
    u = np.where(fan, 2 / (g + 1) * (cL + (g - 1) / 2 * uL + xi), u)
    P = np.where(fan, PL * (2 / (g + 1) + (g - 1) / ((g + 1) * cL) * (uL - xi)) ** (2 * g / (g - 1)), P)
    m = xi < uL - cL
    u[m], P[m] = uL, PL
    m = (xi >= us - csL) & (xi < us); rho[m], u[m], P[m] = rsL, us, Ps
    m = (xi >= us) & (xi < S); rho[m], u[m], P[m] = rsR, us, Ps
    m = xi >= S; rho[m], u[m], P[m] = rR, uR, PR
    return rho, u, P


def run_sod(n, order):
    h = Hydro1D(0.0, 1.0, n, order=order, geometry="planar", gamma=1.4,
                inner_bc="outflow", outer_bc="outflow", cfl=0.4)
    rho = np.where(h.r < 0.5, 1.0, 0.125)
    P = np.where(h.r < 0.5, 1.0, 0.1)
    h.set_primitive(rho, 0.0, P)
    h.evolve(0.2)
    rho, u, P = h.primitive()
    re, ue, Pe = sod_exact(h.r, 0.2)
    return np.mean(np.abs(rho - re))


def test_sod_convergence():
    for order in (1, 2):
        e = [run_sod(n, order) for n in (100, 200, 400)]
        print("Sod L1(rho) order", order, e)
        assert e[2] < e[1] < e[0]
    assert run_sod(400, 2) < run_sod(400, 1)


def test_cc85_injection():
    """Spherical steady wind from resolved uniform injection vs CC85."""
    R, mdot, edot = 1.0, 1.0, 0.5   # v_inf^2 = 2 edot/mdot = 1
    h = Hydro1D(0.0, 10.0, 400, order=2, inner_bc="reflect", outer_bc="outflow",
                source=uniform_injection_source(R, mdot, edot),
                rho_floor=1e-14, P_floor=1e-16)
    h.set_primitive(1e-6, 0.0, 1e-8)
    h.evolve(40.0)
    rho, u, P = h.primitive()
    re, ue, Pe = wind_state(h.r, mdot, edot, R)
    # exclude the sonic point r = R where the exact solution has a kink
    m = (h.r > 0.2) & (h.r < 8) & (np.abs(h.r - R) > 0.2)
    err_rho = np.max(np.abs(rho[m] / re[m] - 1))
    err_u = np.max(np.abs(u[m] / ue[m] - 1))
    print("CC85 max rel err rho=%.3e u=%.3e" % (err_rho, err_u))
    assert err_rho < 0.01 and err_u < 0.01


if __name__ == "__main__":
    test_sod_convergence()
    test_cc85_injection()
