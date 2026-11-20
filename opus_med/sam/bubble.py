"""Adiabatic wind-blown bubble semi-analytic model.

Structure (inside -> out):
  (a) free wind, r < R1     : Chevalier & Clegg (1985) steady wind with uniform
                              mass/energy injection inside R_inj.
  (b) shocked wind, R1<r<Rc : uniform pressure P_b(t) (Weaver et al. 1977),
                              velocity u = 11/25 Rc^3/(t r^2) + 4/25 r/t,
                              density from adiabatic Lagrangian history.
  (c) shocked ISM shell, Rc<r<R2 : self-similar adiabatic solution (gamma=5/3)
                              of a strong shock driven by constant luminosity.
  (d) ambient medium, r > R2.
All quantities in CGS internally.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

from .units import m_p, k_B, M_sun, yr, kyr, pc, AU, km, MU_H, MU, GAMMA


# ---------------------------------------------------------------------------
# (c) self-similar shell
# ---------------------------------------------------------------------------
def _shell_rhs(lam, y, gam):
    """d/dlam of (G, V, Pi) for the self-similar shell, R2 ~ t^(3/5)."""
    G, V, Pi = y
    w = V - lam
    # Linear system:
    #  w G' + G V'               = -2 G V / lam
    #  w V' + Pi'/G              = 2/3 V
    #  w Pi'/Pi - gam w G'/G     = 4/3
    A = np.array([[w, G, 0.0],
                  [0.0, w, 1.0 / G],
                  [-gam * w / G, 0.0, w / Pi]])
    b = np.array([-2.0 * G * V / lam, 2.0 / 3.0 * V, 4.0 / 3.0])
    return np.linalg.solve(A, b)


_SIM_CACHE = {}


def shell_similarity(gam=GAMMA, n=4000):
    """Integrate the shell similarity ODEs from the shock (lam=1) inward to the
    contact discontinuity (V = lam).  Returns dict with arrays and constants."""
    if gam in _SIM_CACHE:
        return _SIM_CACHE[gam]
    y0 = [(gam + 1) / (gam - 1), 2 / (gam + 1), 2 / (gam + 1)]

    def contact(lam, y, gam):
        return (y[1] - lam) / lam - 1e-10
    contact.terminal = True

    sol = solve_ivp(_shell_rhs, [1.0, 0.3], y0, args=(gam,), events=contact,
                    rtol=1e-11, atol=1e-13, dense_output=True,
                    max_step=1e-3)
    lam_c = sol.t_events[0][0] if len(sol.t_events[0]) else sol.t[-1]
    lam = np.linspace(1.0, lam_c, n)
    G, V, Pi = sol.sol(lam)
    # energies per (4 pi R2^3 rho0 U^2)
    ek = np.trapz(0.5 * G * V**2 * lam**2, lam) * -1.0
    et = np.trapz(Pi / (gam - 1) * lam**2, lam) * -1.0
    Pi_c = Pi[-1]
    # bubble thermal energy per (4 pi R2^3 rho0 U^2): Pi_c/(gam-1) * lam_c^3/3
    eb = Pi_c / (gam - 1) * lam_c**3 / 3.0
    K = 4 * np.pi * (ek + et + eb)
    # L t = R2^3 rho0 U^2 K, U = 3/5 R2/t  ->  R2 = xi0 (L t^3/rho0)^(1/5)
    xi0 = (25.0 / (9.0 * K))**0.2
    out = dict(sol=sol, lam=lam[::-1], G=G[::-1], V=V[::-1], Pi=Pi[::-1],
               lam_c=lam_c, Pi_c=Pi_c, xi0=xi0,
               frac_kin_shell=4*np.pi*ek/K, frac_th_shell=4*np.pi*et/K,
               frac_bubble=4*np.pi*eb/K)
    _SIM_CACHE[gam] = out
    return out


# ---------------------------------------------------------------------------
# (a) Chevalier & Clegg free wind
# ---------------------------------------------------------------------------
def _cc85_mach(x, gam=GAMMA):
    """Mach number as function of x = r/R_inj (vectorised)."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    M = np.empty_like(x)

    def f_in(M, xx):
        a = ((3 * gam + 1 / M**2) / (1 + 3 * gam))**(-(3 * gam + 1) / (5 * gam + 1))
        b = ((gam - 1 + 2 / M**2) / (1 + gam))**((gam + 1) / (2 * (5 * gam + 1)))
        return np.log(a * b) - np.log(xx)

    def f_out(M, xx):
        lhs = (M**(2 / (gam - 1)) *
               ((gam - 1 + 2 / M**2) / (1 + gam))**((gam + 1) / (2 * (gam - 1))))
        return np.log(lhs) - 2 * np.log(xx)

    for i, xx in enumerate(x):
        if xx <= 0:
            M[i] = 0.0
        elif xx < 1:
            M[i] = brentq(f_in, 1e-12, 1.0, args=(xx,), xtol=1e-14)
        elif xx == 1:
            M[i] = 1.0
        else:
            M[i] = brentq(f_out, 1.0, 1e6, args=(xx,), xtol=1e-14)
    return M


def cc85_wind(r, Mdot, L, R_inj, gam=GAMMA):
    """Return rho, v, P (cgs) of the CC85 steady wind at radii r [cm]."""
    r = np.atleast_1d(np.asarray(r, dtype=float))
    vinf = np.sqrt(2 * L / Mdot)
    M = _cc85_mach(r / R_inj, gam)
    v = vinf / np.sqrt(1 + 2 / ((gam - 1) * np.maximum(M, 1e-300)**2))
    v = np.where(M > 0, v, 0.0)
    inside = r < R_inj
    rhov = np.where(inside, Mdot * r / (4 * np.pi * R_inj**3),
                    Mdot / (4 * np.pi * np.maximum(r, 1e-300)**2))
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = rhov / v
        c2 = (v / M)**2
    # centre: v->0, M->0; use limit c^2 -> (gam-1) vinf^2/2
    c2 = np.where(M > 0, c2, (gam - 1) * 0.5 * vinf**2)
    # rho at r=0: limit of rhov/v with v ~ r -> finite
    if np.any(M == 0):
        eps = 1e-6 * R_inj
        rho0, _, _ = cc85_wind(np.array([eps]), Mdot, L, R_inj, gam)
        rho = np.where(M > 0, rho, rho0[0])
    P = rho * c2 / gam
    return rho, v, P


# ---------------------------------------------------------------------------
# Full model
# ---------------------------------------------------------------------------
class WindBubble:
    """Adiabatic wind-blown bubble.

    Parameters (CGS unless noted)
    ----------
    Mdot_msun_yr : wind mass-loss rate [Msun/yr]
    L : wind mechanical luminosity [erg/s]
    nH : ambient hydrogen number density [cm^-3]
    R_inj : injection radius [cm]
    T_amb : ambient temperature [K]
    """

    def __init__(self, Mdot_msun_yr, L, nH, R_inj=100 * AU, T_amb=100.0,
                 gam=GAMMA, mu_H=MU_H, mu=MU):
        self.Mdot = Mdot_msun_yr * M_sun / yr
        self.L = L
        self.nH = nH
        self.R_inj = R_inj
        self.T_amb = T_amb
        self.gam = gam
        self.mu_H = mu_H
        self.mu = mu
        self.rho0 = mu_H * m_p * nH
        self.P_amb = self.rho0 / (mu * m_p) * k_B * T_amb
        self.vinf = np.sqrt(2 * L / self.Mdot)
        self.sim = shell_similarity(gam)

    # ---- characteristic radii / pressure --------------------------------
    def R2(self, t):
        return self.sim["xi0"] * (self.L * t**3 / self.rho0)**0.2

    def Rc(self, t):
        return self.sim["lam_c"] * self.R2(t)

    def U(self, t):
        return 0.6 * self.R2(t) / t

    def P_bubble(self, t):
        """Hot (shocked wind) gas pressure [dyn cm^-2]."""
        return self.rho0 * self.U(t)**2 * self.sim["Pi_c"]

    def _jump(self, r1, t):
        """Post-shock state of the reverse shock located at r1 at time t."""
        g = self.gam
        rho1, v1, P1 = cc85_wind(np.array([r1]), self.Mdot, self.L, self.R_inj, g)
        rho1, v1, P1 = rho1[0], v1[0], P1[0]
        vs = 0.4 * r1 / t  # R1 ~ t^(2/5)
        u1 = v1 - vs
        M1sq = rho1 * u1**2 / (g * P1)
        P2 = P1 * (2 * g * M1sq - (g - 1)) / (g + 1)
        rho2 = rho1 * (g + 1) * M1sq / ((g - 1) * M1sq + 2)
        v2 = vs + u1 * rho1 / rho2
        return rho2, v2, P2

    def R1(self, t):
        """Reverse shock radius: post-shock pressure equals P_bubble."""
        Pb = self.P_bubble(t)
        f = lambda lr: np.log(self._jump(np.exp(lr), t)[2] / Pb)
        lo = np.log(1.0001 * self.R_inj)
        hi = np.log(0.999 * self.Rc(t))
        if f(lo) < 0:
            # wind not strong enough to be supersonic shocked: shock at R_inj
            return np.exp(lo)
        return np.exp(brentq(f, lo, hi, xtol=1e-12))

    # ---- profiles ---------------------------------------------------------
    def profiles(self, r, t):
        """rho [g/cc], v [cm/s], P [dyn/cm^2], T [K] at radii r [cm], time t [s]."""
        r = np.atleast_1d(np.asarray(r, dtype=float))
        g = self.gam
        R1, Rc, R2 = self.R1(t), self.Rc(t), self.R2(t)
        Pb = self.P_bubble(t)
        rho = np.full_like(r, self.rho0)
        v = np.zeros_like(r)
        P = np.full_like(r, self.P_amb)

        # (a) free wind
        a = r < R1
        if np.any(a):
            rho[a], v[a], P[a] = cc85_wind(r[a], self.Mdot, self.L, self.R_inj, g)

        # (b) shocked wind
        b = (r >= R1) & (r < Rc)
        if np.any(b):
            rb = r[b]
            v[b] = 11.0 / 25.0 * Rc**3 / (t * rb**2) + 4.0 / 25.0 * rb / t
            P[b] = Pb
            rho[b] = self._shocked_wind_density(rb, t)

        # (c) shell
        c = (r >= Rc) & (r < R2)
        if np.any(c):
            lam = r[c] / R2
            G, V, Pi = self.sim["sol"].sol(lam)
            U = self.U(t)
            rho[c] = self.rho0 * G
            v[c] = U * V
            P[c] = self.rho0 * U**2 * Pi + 0.0
        T = P * self.mu * m_p / (rho * k_B)
        return rho, v, P, T

    def _shocked_wind_density(self, rb, t):
        """Density in region (b) from the Lagrangian history of each element.

        With x = r/Rc, pathlines obey dx/dlnt = 11/25 (1-x^3)/x^2, so
        (1 - x^3) t^(33/25) = const.  An element at (r, t) was shocked at t'
        where x1(t') = R1(t')/Rc(t').  Since P_b ~ t^(-4/5) and the flow is
        adiabatic: rho(r,t) = rho_s(t') (t/t')^(-12/25)."""
        Rc = self.Rc(t)
        x = np.minimum(rb / Rc, 1 - 1e-12)
        A = (1 - x**3) * t**(33.0 / 25.0)
        rho = np.empty_like(rb)
        # x1(t') ~ t'^(-1/5) approximately; compute exactly via R1
        for i, Ai in enumerate(A):
            def f(ltp):
                tp = np.exp(ltp)
                x1 = self.R1(tp) / self.Rc(tp)
                return np.log(1 - x1**3) + 33.0 / 25.0 * ltp - np.log(Ai)
            lt_hi = np.log(t)
            if f(lt_hi) <= 0:
                tp = t
            else:
                lt_lo = lt_hi - 1.0
                while f(lt_lo) > 0:
                    lt_lo -= 1.0
                tp = np.exp(brentq(f, lt_lo, lt_hi, xtol=1e-10))
            r1p = self.R1(tp)
            rho_s, _, _ = self._jump(r1p, tp)
            rho[i] = rho_s * (t / tp)**(-12.0 / 25.0)
        return rho

    # ---- convenience ------------------------------------------------------
    def profiles_physical_units(self, r_pc, t_kyr):
        """Return r [pc], rho [g/cc], v [km/s], P/k [K/cc], T [K]."""
        rho, v, P, T = self.profiles(np.asarray(r_pc) * pc, t_kyr * kyr)
        return np.asarray(r_pc), rho, v / km, P / k_B, T

    def summary(self, t_kyr):
        t = t_kyr * kyr
        return dict(R1_pc=self.R1(t) / pc, Rc_pc=self.Rc(t) / pc,
                    R2_pc=self.R2(t) / pc, Pb_over_k=self.P_bubble(t) / k_B,
                    vinf_kms=self.vinf / km)
