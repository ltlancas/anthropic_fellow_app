"""Semi-analytic model of an adiabatic wind-blown bubble.

Regions (inside -> out):
  1. free wind   r < R1      : Chevalier & Clegg (1985) steady solution
  2. inner shock r = R1      : Rankine-Hugoniot, P_post = P_b
  3. hot bubble  R1 < r < Rc : uniform pressure P_b, Weaver et al. (1977)
                                velocity field, density from entropy history
  4. shell       Rc < r < R2 : self-similar adiabatic shocked ambient gas
  5. ambient     r > R2

All internal computations in cgs.  See assumptions.md.
"""
import numpy as np
from scipy.optimize import brentq

from . import cc85
from .shell import shell_solution
from .units import GAMMA, AU, MSUN, YR, nH_to_rho, pressure_from_T, temperature

ALPHA = 0.6           # R2 ∝ t^(3/5)
ALPHA_1 = 0.4         # R1 ∝ t^(2/5) (ram-pressure balance with P_b ∝ t^-4/5)


class WindBubble:
    def __init__(self, Mdot_msun_yr, L_erg_s, nH, R_inj_cm=100 * AU,
                 T_amb=100.0, energy_correction=True):
        self.Mdot = Mdot_msun_yr * MSUN / YR
        self.L = L_erg_s
        self.rho0 = nH_to_rho(nH)
        self.P0 = pressure_from_T(self.rho0, T_amb)
        self.R_inj = R_inj_cm
        self.vinf = np.sqrt(2 * self.L / self.Mdot)
        self.energy_correction = energy_correction
        self.ss = shell_solution(ALPHA, GAMMA)
        a, g, s = ALPHA, GAMMA, self.ss
        self.xi2_ss = (1 / (4 * np.pi * a**2 * (
            s["I_shell"] + s["lam_c"]**3 * s["Pi_c"] / (3 * (g - 1))))) ** 0.2

    # ------------------------------------------------------------------
    def _post_shock(self, R1, t):
        """RH jump of the CC85 wind at a shock at R1 moving at 0.4 R1/t."""
        g = GAMMA
        rho1, v1, P1 = cc85.wind_profile(R1, self.Mdot, self.L, self.R_inj)
        rho1, v1, P1 = rho1[0], v1[0], P1[0]
        Vs = ALPHA_1 * R1 / t
        u1 = v1 - Vs
        M2 = u1**2 / (g * P1 / rho1)
        comp = (g + 1) * M2 / ((g - 1) * M2 + 2)
        P2 = P1 * (2 * g * M2 - (g - 1)) / (g + 1)
        rho2 = rho1 * comp
        v2 = Vs + u1 / comp
        return rho2, v2, P2

    def _inner_shock(self, Pb, Rc, t):
        f = lambda R1: self._post_shock(R1, t)[2] - Pb
        lo = self.R_inj * 1.0001
        if f(lo) < 0:
            return lo
        if f(Rc * 0.999) > 0:
            return Rc * 0.999
        return brentq(f, lo, Rc * 0.999, rtol=1e-10)

    def _bubble_density(self, r, R1, Rc, rho2):
        """Density in the hot bubble from adiabatic Lagrangian history."""
        y = np.clip((r / Rc) ** 3, 0, 1)
        y0 = (R1 / Rc) ** 3
        lo = np.full_like(y, np.log(y0 ** (5 / 3)))
        hi = np.zeros_like(y)
        f = lambda ls: np.exp(1.32 * ls) - y0 * np.exp(0.72 * ls) - (1 - y)
        for _ in range(80):                     # vectorised bisection
            mid = 0.5 * (lo + hi)
            pos = f(mid) > 0
            hi = np.where(pos, mid, hi)
            lo = np.where(pos, lo, mid)
        s = np.exp(0.5 * (lo + hi))
        return rho2 * s ** (-8.0 / 25.0)

    # ------------------------------------------------------------------
    def _structure(self, t, xi2):
        g, s = GAMMA, self.ss
        R2 = xi2 * (self.L / self.rho0) ** 0.2 * t ** 0.6
        R2dot = ALPHA * R2 / t
        Pb = self.rho0 * R2dot**2 * s["Pi_c"]
        Rc = s["lam_c"] * R2
        R1 = self._inner_shock(Pb, Rc, t)
        rho2, v2, _ = self._post_shock(R1, t)
        return dict(t=t, xi2=xi2, R2=R2, R2dot=R2dot, Pb=Pb, Rc=Rc, R1=R1,
                    rho2=rho2, v2=v2)

    def _energy(self, st, n=2000):
        g, s = GAMMA, self.ss
        E_shell = 4 * np.pi * st["R2"]**3 * self.rho0 * st["R2dot"]**2 * s["I_shell"]
        # hot bubble
        rb = np.linspace(st["R1"], st["Rc"], n)
        rho_b = self._bubble_density(rb, st["R1"], st["Rc"], st["rho2"])
        v_b = self._bubble_velocity(rb, st)
        E_b = np.trapezoid(4 * np.pi * rb**2 * (st["Pb"] / (g - 1) + 0.5 * rho_b * v_b**2), rb)
        # free wind
        rw = np.concatenate([np.linspace(1e-6 * self.R_inj, self.R_inj, 200),
                             np.geomspace(self.R_inj, st["R1"], 400)[1:]])
        rho_w, v_w, P_w = cc85.wind_profile(rw, self.Mdot, self.L, self.R_inj)
        E_w = np.trapezoid(4 * np.pi * rw**2 * (0.5 * rho_w * v_w**2 + P_w / (g - 1)), rw)
        return E_shell + E_b + E_w

    def _bubble_velocity(self, r, st):
        t, Rc = st["t"], st["Rc"]
        return (r / t) * (4.0 / 25.0 + (11.0 / 25.0) * (Rc / r) ** 3)

    def state(self, t):
        """Global bubble parameters at time t [s]."""
        if not self.energy_correction:
            return self._structure(t, self.xi2_ss)
        f = lambda xi2: self._energy(self._structure(t, xi2)) / (self.L * t) - 1
        xi2 = brentq(f, 0.7 * self.xi2_ss, 1.3 * self.xi2_ss, rtol=1e-8)
        return self._structure(t, xi2)

    # ------------------------------------------------------------------
    def profile(self, r, t, st=None):
        """rho [g/cm^3], v [cm/s], P [dyn/cm^2], T [K] at radii r [cm]."""
        st = st or self.state(t)
        r = np.asarray(r, dtype=float)
        rho = np.full_like(r, self.rho0)
        v = np.zeros_like(r)
        P = np.full_like(r, self.P0)
        s = self.ss

        m = r < st["R1"]
        if m.any():
            rho[m], v[m], P[m] = cc85.wind_profile(r[m], self.Mdot, self.L, self.R_inj)
        m = (r >= st["R1"]) & (r < st["Rc"])
        if m.any():
            rho[m] = self._bubble_density(r[m], st["R1"], st["Rc"], st["rho2"])
            v[m] = self._bubble_velocity(r[m], st)
            P[m] = st["Pb"]
        m = (r >= st["Rc"]) & (r < st["R2"])
        if m.any():
            lam = r[m] / st["R2"]
            rho[m] = self.rho0 * np.interp(lam, s["lam"], s["G"])
            v[m] = st["R2dot"] * np.interp(lam, s["lam"], s["V"])
            P[m] = self.rho0 * st["R2dot"]**2 * np.interp(lam, s["lam"], s["Pi"])
        return rho, v, P, temperature(P, rho)
