"""Semi-analytic model for an adiabatic (non-radiative, non-conducting) wind-blown bubble.

Structure (inside -> out), see sam_outputs/assumptions.md for derivations:
  1. r < R1   : steady free wind with distributed mass/energy injection inside R_inj
                (Chevalier & Clegg 1985, "CC85").
  2. R1       : wind termination (reverse) shock, Rankine-Hugoniot jump.
  3. R1<r<Rc  : shocked wind ("hot bubble"). Quasi-steady subsonic layer behind R1
                matched to a uniform-pressure Lagrangian (entropy-history) solution.
  4. Rc       : contact discontinuity.
  5. Rc<r<R2  : shocked ambient medium, exact self-similar solution for an
                energy-driven (R2 ~ t^{3/5}) adiabatic blast wave with a piston
                (Weaver et al. 1977 adiabatic phase; Castor, McCray & Weaver 1975).
  6. r > R2   : undisturbed ambient medium.
"""
import numpy as np
from functools import lru_cache
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

from . import units as u

G = u.GAMMA


# ----------------------------------------------------------------------------
# 1. Self-similar shocked-ambient shell
# ----------------------------------------------------------------------------
def _shell_rhs(xi, y, delta=-2.0 / 3.0, gam=G):
    Uv, Gv, Pv = y
    W = Uv - xi
    c2 = gam * Pv / Gv
    dU = (2.0 * c2 * Uv / xi + 2.0 * delta * Pv / Gv - delta * Uv * W) / (W * W - c2)
    dG = -Gv * (dU + 2.0 * Uv / xi) / W
    dP = Pv * (-gam * (dU + 2.0 * Uv / xi) - 2.0 * delta) / W
    return [dU, dG, dP]


@lru_cache(maxsize=None)
def shell_solution(gam=G, n=4000):
    """Integrate the similarity ODEs from the shock (xi=1) inward to the contact.

    Returns dict with xi grid (increasing), U, Gd, Pi, xi_c, Pi_c, A (radius constant).
    """
    y0 = [2.0 / (gam + 1.0), (gam + 1.0) / (gam - 1.0), 2.0 / (gam + 1.0)]

    def ev(xi, y):
        return (y[0] - xi) / xi - 1e-9
    ev.terminal = True

    sol = solve_ivp(_shell_rhs, [1.0, 0.3], y0, events=ev, rtol=1e-11, atol=1e-13,
                    method="DOP853", dense_output=True, max_step=1e-3)
    xi_c = sol.t_events[0][0] if len(sol.t_events[0]) else sol.t[-1]
    xi = np.linspace(xi_c, 1.0, n)
    Y = sol.sol(xi)
    Uv, Gv, Pv = Y
    # Shell energy integral  I_s = int xi^2 (G U^2/2 + Pi/(gam-1)) dxi
    f = lambda x: (lambda y: x * x * (0.5 * y[1] * y[0] ** 2 + y[2] / (gam - 1.0)))(sol.sol(x))
    I_s = quad(f, xi_c, 1.0, limit=200, epsabs=1e-12)[0]
    Pi_c = float(sol.sol(xi_c)[2])
    I_b = xi_c ** 3 * Pi_c / (3.0 * (gam - 1.0))
    J = I_s + I_b
    # L t = 4 pi rho0 R^3 Rdot^2 J,  R = A (L/rho0)^{1/5} t^{3/5}
    A = (25.0 / (9.0 * 4.0 * np.pi * J)) ** 0.2
    return dict(xi=xi, U=Uv, G=Gv, Pi=Pv, xi_c=float(xi_c), Pi_c=Pi_c, A=float(A),
                I_s=I_s, I_b=I_b, f_bubble=I_b / J, sol=sol)


# ----------------------------------------------------------------------------
# 2. CC85 free wind
# ----------------------------------------------------------------------------
def _cc85_mach(x, gam=G):
    """Mach number at x = r/R_inj (subsonic inside, supersonic outside)."""
    if x <= 0:
        return 0.0
    if x < 1.0:
        a1 = -(3 * gam + 1) / (5 * gam + 1)
        a2 = (gam + 1) / (2 * (5 * gam + 1))
        f = lambda M: ((3 * gam + 1 / M ** 2) / (1 + 3 * gam)) ** a1 * \
            ((gam - 1 + 2 / M ** 2) / (1 + gam)) ** a2 - x
        return brentq(f, 1e-12, 1.0, xtol=1e-14, rtol=1e-13)
    if x == 1.0:
        return 1.0
    a = (gam + 1) / (2 * (gam - 1))
    f = lambda M: np.log(M) - a * np.log(((gam - 1) * M * M + 2) / (gam + 1)) + 2 * np.log(x)
    return brentq(f, 1.0, 1e8, xtol=1e-14, rtol=1e-13)


def cc85_state(r, mdot, L, R_inj, gam=G):
    """Free-wind (rho, v, P) at radius r (scalar)."""
    vinf = np.sqrt(2 * L / mdot)
    M = _cc85_mach(r / R_inj, gam)
    c2 = vinf ** 2 / (M * M + 2.0 / (gam - 1.0))
    v = M * np.sqrt(c2)
    if r < R_inj:
        q = mdot / (4.0 / 3.0 * np.pi * R_inj ** 3)
        if r == 0:
            # rho finite at centre: rho*v = q r/3, v ~ M c, M ~ r/R_inj * const
            r = 1e-8 * R_inj
            M = _cc85_mach(r / R_inj, gam)
            c2 = vinf ** 2 / (M * M + 2.0 / (gam - 1.0))
            v = M * np.sqrt(c2)
        rho = q * r / 3.0 / v
    else:
        rho = mdot / (4 * np.pi * r * r * v)
    return rho, v, rho * c2 / gam


# ----------------------------------------------------------------------------
# 3. Reverse shock and quasi-steady layer
# ----------------------------------------------------------------------------
def rh_jump(rho1, v1, P1, gam=G):
    """Stationary normal shock; returns post-shock (rho, v, P)."""
    c1 = np.sqrt(gam * P1 / rho1)
    M1 = v1 / c1
    r = (gam + 1) * M1 ** 2 / ((gam - 1) * M1 ** 2 + 2)
    p = (2 * gam * M1 ** 2 - (gam - 1)) / (gam + 1)
    return rho1 * r, v1 / r, P1 * p


def _subsonic_mach(ratio, gam=G):
    """Solve M * (((g-1)M^2+2)/(g+1))^{-(g+1)/(2(g-1))} = ratio on subsonic branch."""
    a = (gam + 1) / (2 * (gam - 1))
    f = lambda M: np.log(M) - a * np.log(((gam - 1) * M * M + 2) / (gam + 1)) - np.log(ratio)
    return brentq(f, 1e-14, 1.0, xtol=1e-16, rtol=1e-13)


# ----------------------------------------------------------------------------
# Main model
# ----------------------------------------------------------------------------
class WindBubble:
    """Adiabatic wind bubble.

    Parameters (cgs unless noted): mdot [g/s], L [erg/s], nH [cm^-3],
    R_inj [cm], T_amb [K] (only used to fill the ambient pressure outside R2;
    the dynamics assume a strong outer shock).
    """

    def __init__(self, mdot, L, nH, R_inj=100 * u.AU, T_amb=10.0, mu=u.MU_H, gam=G):
        self.mdot, self.L, self.nH = mdot, L, nH
        self.R_inj, self.T_amb, self.mu, self.gam = R_inj, T_amb, mu, gam
        self.rho0 = u.nH_to_rho(nH, mu)
        self.P0 = u.pressure_from_T(self.rho0, T_amb, mu)
        self.vinf = np.sqrt(2 * L / mdot)
        self.ss = shell_solution(gam)
        # characteristic radius/time at which swept-up mass ~ wind mass
        self.R_ch = np.sqrt(mdot / (self.rho0 * self.vinf))
        self.t_ch = self.R_ch / self.vinf

    @classmethod
    def from_astro(cls, mdot_msun_yr, L, nH, R_inj_AU=100.0, **kw):
        return cls(u.msun_per_yr_to_cgs(mdot_msun_yr), L, nH, R_inj_AU * u.AU, **kw)

    # --- global quantities -------------------------------------------------
    def R2(self, t):
        return self.ss["A"] * (self.L / self.rho0) ** 0.2 * t ** 0.6

    def R2dot(self, t):
        return 0.6 * self.R2(t) / t

    def Rc(self, t):
        return self.ss["xi_c"] * self.R2(t)

    def P_bubble(self, t):
        return self.rho0 * self.R2dot(t) ** 2 * self.ss["Pi_c"]

    def _layer(self, R1):
        """Post-shock state and stagnation pressure for reverse shock at R1."""
        rho1, v1, P1 = cc85_state(R1, self.mdot, self.L, self.R_inj, self.gam)
        rs, vs, Ps = rh_jump(rho1, v1, P1, self.gam)
        Ms = vs / np.sqrt(self.gam * Ps / rs)
        Pst = Ps * (1 + 0.5 * (self.gam - 1) * Ms ** 2) ** (self.gam / (self.gam - 1))
        return rs, vs, Ps, Ms, Pst

    def R1(self, t):
        Pb = self.P_bubble(t)
        f = lambda lr: np.log(self._layer(np.exp(lr))[4] / Pb)
        lo, hi = np.log(1.0001 * self.R_inj), np.log(1e3 * self.R2(t))
        if f(lo) < 0:
            raise ValueError("bubble pressure exceeds wind pressure at R_inj")
        return float(np.exp(brentq(f, lo, hi, xtol=1e-12)))

    def summary(self, t):
        R1 = self.R1(t)
        return dict(t=t, R1=R1, Rc=self.Rc(t), R2=self.R2(t), P_b=self.P_bubble(t),
                    R2dot=self.R2dot(t), t_over_tch=t / self.t_ch,
                    valid=bool(R1 < self.Rc(t)))

    # --- profiles ----------------------------------------------------------
    def profile(self, r, t):
        """Return dict of rho, v, P, T arrays at radii r (cm) and time t (s)."""
        r = np.atleast_1d(np.asarray(r, dtype=float))
        gam = self.gam
        R1, Rc, R2 = self.R1(t), self.Rc(t), self.R2(t)
        # Outside validity (R1 >= Rc, t <~ t_ch): no hot bubble; the free wind is
        # taken to extend to the contact discontinuity (documented fallback).
        R1 = min(R1, Rc)
        Pb, R2d = self.P_bubble(t), self.R2dot(t)
        rho = np.empty_like(r); v = np.empty_like(r); P = np.empty_like(r)

        # (1) free wind
        m = r < R1
        for i in np.where(m)[0]:
            rho[i], v[i], P[i] = cc85_state(r[i], self.mdot, self.L, self.R_inj, gam)

        # (3) hot bubble
        if R1 < Rc:
            rs, vs, Ps, Ms, Pst = self._layer(R1)
            K = Ps / rs ** gam
            a = (gam + 1) / (2 * (gam - 1))
            C0 = Ms * (((gam - 1) * Ms ** 2 + 2) / (gam + 1)) ** (-a)  # times R1^2
            h0 = 0.5 * vs ** 2 + gam / (gam - 1) * Ps / rs              # Bernoulli const
            rho_st = (Pst / K) ** (1 / gam)
            V1 = 4 * np.pi / 3 * R1 ** 3
            dV = 4 * np.pi / 3 * (Rc ** 3 - R1 ** 3)
            dV1dt = 4 * np.pi * R1 ** 2 * 0.4 * R1 / t
            uo_1 = (dV1dt + (33.0 / 25.0) * dV / t)  # outer 4 pi r^2 u at x -> 1
            mb = (r >= R1) & (r < Rc)
            for i in np.where(mb)[0]:
                ri = r[i]
                # inner quasi-steady layer
                M = _subsonic_mach(C0 * (R1 / ri) ** 2, gam)
                c2 = h0 / (0.5 * M * M + 1 / (gam - 1))
                rin = (c2 / (gam * K)) ** (1 / (gam - 1))   # isentropic: c^2 = gam K rho^(gam-1)
                pin = K * rin ** gam
                vin = M * np.sqrt(c2)
                # outer Lagrangian uniform-pressure solution
                x = max(1.0 - (4 * np.pi / 3 * (ri ** 3 - R1 ** 3)) / dV, 1e-12) ** (25.0 / 33.0)
                x33 = x ** (33.0 / 25.0)
                uo = (dV1dt + (dV / t) * (1.8 * (1 - x33) + (33.0 / 25.0) * x33))
                rho[i] = rin * x ** (-8.0 / 25.0)
                P[i] = pin
                v[i] = vin + (uo - uo_1) / (4 * np.pi * ri ** 2)
        # (5) shell
        ms = (r >= Rc) & (r < R2) & (r >= R1)
        if ms.any():
            xi = r[ms] / R2
            Y = self.ss["sol"].sol(np.clip(xi, self.ss["xi_c"], 1.0))
            v[ms] = R2d * Y[0]
            rho[ms] = self.rho0 * Y[1]
            P[ms] = np.maximum(self.rho0 * R2d ** 2 * Y[2], self.P0)
        # (6) ambient
        ma = r >= R2
        rho[ma] = self.rho0; v[ma] = 0.0; P[ma] = self.P0
        T = u.temperature(rho, P, self.mu)
        return dict(r=r, rho=rho, v=v, P=P, T=T)
