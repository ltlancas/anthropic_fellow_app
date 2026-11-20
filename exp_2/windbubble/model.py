"""Semi-analytic model of an adiabatic wind-blown bubble in a uniform medium.

Structure (inside -> out):

  I   r < R1      free wind: Chevalier & Clegg (1985) steady solution for
                  uniform mass+energy injection inside r_inj.
  II  R1 < r < Rc shocked wind: isobaric (P = P_b(t)); each fluid element
                  keeps the entropy it received at the reverse shock.
  III Rc < r < R2 shocked ambient medium: quasi-self-similar adiabatic shell
                  (see shell.py) parametrised by the instantaneous
                  deceleration parameter q = R2 R2'' / R2'^2.
  IV  r > R2      undisturbed ambient medium (rho_0, T_0, u = 0).

Time evolution: ODEs for y = (R2, V2, R1, E_II, M_II)

  dR2/dt   = V2
  dV2/dt   = q V2^2 / R2,     q from the closure  P_b = Pi_c(q) rho_0 V2^2
  dR1/dt   = u_w - M1 c_w,    M1 from Rankine-Hugoniot with P_post = P_b
  dE_II/dt = F_E(R1) - P_b 4 pi Rc^2 dRc/dt,   dRc/dt = xi_c(q) V2
  dM_II/dt = 4 pi R1^2 rho_w (u_w - dR1/dt)

and at each instant q is found from the energy content of region II:
  E_II = P_b V_II/(g-1) + K_II,  V_II = 4pi/3 (Rc^3 - R1^3),  Rc = xi_c(q) R2.

See assumptions.md for the derivation.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import units as U
from .cc85 import mach as cc85_mach
from .shell import ShellTable

_SHELL_TABLE = None


def shell_table():
    global _SHELL_TABLE
    if _SHELL_TABLE is None:
        _SHELL_TABLE = ShellTable()
    return _SHELL_TABLE


class CC85Wind:
    """Fast tabulated evaluation of the CC85 wind."""

    def __init__(self, Mdot, Edot, R, g=U.GAMMA):
        self.Mdot, self.Edot, self.R, self.g = Mdot, Edot, R, g
        self.vinf = np.sqrt(2 * Edot / Mdot)
        self._lx = np.linspace(np.log(1e-4), np.log(1e9), 6000)
        self._lM = np.log(cc85_mach(np.exp(self._lx), 1.0, g))

    def mach(self, r):
        x = np.maximum(np.asarray(r, dtype=float) / self.R, 1e-4)
        return np.exp(np.interp(np.log(x), self._lx, self._lM))

    def state(self, r):
        """(rho, u, P) at radius r (cgs)."""
        r = np.asarray(r, dtype=float)
        g = self.g
        M = self.mach(r)
        u = self.vinf / np.sqrt(1.0 + 2.0 / ((g - 1) * M**2))
        rr = np.maximum(r, 1e-4 * self.R)
        mdot_enc = self.Mdot * np.minimum(rr / self.R, 1.0) ** 3
        rho = mdot_enc / (4 * np.pi * rr**2 * u)
        P = rho * (u / M) ** 2 / g
        u = np.where(r > 1e-4 * self.R, u, u * r / (1e-4 * self.R))
        return rho, u, P


class WindBubble:
    """Adiabatic wind bubble in a uniform medium.

    Parameters (astrophysical units)
    ----------
    mdot : wind mass-loss rate [Msun/yr]
    Lw   : wind mechanical luminosity [erg/s]
    nH   : ambient hydrogen number density [cm^-3]
    r_inj_au : radius of the uniform injection region [AU]
    T0   : ambient temperature [K] (only used for region IV; P_0 is neglected
           in the dynamics)
    """

    def __init__(self, mdot, Lw, nH, r_inj_au=100.0, T0=100.0,
                 mu_H=U.MU_H, mu=U.MU, g=U.GAMMA):
        self.Mdot = U.mdot_cgs(mdot)
        self.L = float(Lw)
        self.nH = nH
        self.rho0 = U.rho_from_nH(nH, mu_H)
        self.r_inj = r_inj_au * U.AU
        self.T0 = T0
        self.mu = mu
        self.g = g
        self.P0 = U.pressure_from_T(T0, self.rho0, mu)
        self.wind = CC85Wind(self.Mdot, self.L, self.r_inj, g)
        self.vw = self.wind.vinf
        self.tab = shell_table()
        # characteristic time at which swept mass ~ wind mass
        self.t_ch = np.sqrt(self.Mdot / (self.rho0 * self.vw**3))
        self.sol = None

    # ------------------------------------------------------------------ #
    # reverse shock jump
    def _reverse_shock(self, R1, Pb):
        """Given R1 and post-shock pressure, return dict of shock quantities."""
        g = self.g
        rw, uw, Pw = (float(a) for a in self.wind.state(R1))
        cw = np.sqrt(g * Pw / rw)
        Pb = max(Pb, Pw)        # no shock (M1 = 1) if Pb <= P_wind
        M1sq = max(((g + 1) * Pb / Pw + (g - 1)) / (2 * g), 1.0)
        M1 = np.sqrt(M1sq)
        w1 = M1 * cw                       # upstream speed in shock frame
        R1dot = uw - w1
        comp = (g + 1) * M1sq / ((g - 1) * M1sq + 2)
        rho2 = rw * comp
        u2 = R1dot + w1 / comp
        mflux = 4 * np.pi * R1**2 * rw * w1
        Ew = rw * (0.5 * uw**2) + Pw / (g - 1)
        Eflux = 4 * np.pi * R1**2 * (Ew * w1 + Pw * uw)
        K = Pb / rho2**g
        return dict(R1dot=R1dot, rho2=rho2, u2=u2, mflux=mflux, Eflux=Eflux,
                    K=K, rw=rw, uw=uw, Pw=Pw, Ppost=Pb)

    # ------------------------------------------------------------------ #
    # ------------------------------------------------------------------ #
    @staticmethod
    def _u_coeffs(R1, Rc, u1, uc):
        """Coefficients of r^2 u = a + b r^3 through (R1,u1) and (Rc,uc)."""
        b = (Rc**2 * uc - R1**2 * u1) / (Rc**3 - R1**3)
        a = R1**2 * u1 - b * R1**3
        return a, b

    def _bernoulli_P(self, r, R1, Rc, rs, Rcdot, nint=128):
        """Pressure and velocity in region II.

        Velocity: r^2 u = a + b r^3 through (R1, u2) and (Rc, dRc/dt).
        Pressure: quasi-steady momentum balance in the frame of the local
        pattern speed v_p(r) (linear between dR1/dt at R1 and dRc/dt at Rc),
            dP/dr = -rho (u - v_p) du/dr,  rho = rho2 (P/P_post)^(1/g).
        With Z = (P/P_post)^((g-1)/g) this integrates exactly to
            Z(r) = 1 - (g-1)/g (rho2/P_post) int_R1^r (u - v_p) u' dr'.
        For v_p = const this is Bernoulli behind the reverse shock.
        """
        g = self.g
        r = np.asarray(r, dtype=float)
        if Rc <= R1 * (1 + 1e-12):
            return np.full_like(r, rs["Ppost"]), np.full_like(r, rs["u2"])
        a, b = self._u_coeffs(R1, Rc, rs["u2"], Rcdot)
        rg = np.linspace(R1, Rc, nint + 1)
        u = (a + b * rg**3) / rg**2
        du = b - 2 * a / rg**3
        vp = rs["R1dot"] + (rg - R1) / (Rc - R1) * (Rcdot - rs["R1dot"])
        f = (u - vp) * du
        I = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(rg))])
        Z = 1.0 - (g - 1) / g * rs["rho2"] / rs["Ppost"] * I
        Z = np.maximum(Z, 1e-6)
        Pg = rs["Ppost"] * Z ** (g / (g - 1))
        P = np.interp(r, rg, Pg)
        uu = (a + b * r**3) / r**2
        return P, uu

    def _region2(self, R1, Rc, Rcdot, Pc, M2, n=24):
        """Solve for the reverse-shock state given the contact pressure Pc.

        Returns (rs, Eth, KE) where rs contains the shock jump and Ppost."""
        g = self.g
        Ppost = Pc
        for _ in range(50):
            rs = self._reverse_shock(R1, Ppost)
            if Rc <= R1:
                break
            Pcon, _ = self._bernoulli_P(np.array([Rc]), R1, Rc, rs, Rcdot)
            fac = Pcon[0] / rs["Ppost"]
            Pnew = Pc / fac
            if abs(Pnew / Ppost - 1) < 1e-11:
                Ppost = Pnew
                break
            Ppost = Pnew
        rs = self._reverse_shock(R1, Ppost)
        if Rc <= R1:
            return rs, 0.0, 0.0
        x, wq = np.polynomial.legendre.leggauss(n)
        r = 0.5 * (Rc - R1) * x + 0.5 * (Rc + R1)
        P, u = self._bernoulli_P(r, R1, Rc, rs, Rcdot)
        dV = 0.5 * (Rc - R1) * wq * 4 * np.pi * r**2
        Vol = dV.sum()
        Eth = np.sum(P * dV) / (g - 1)
        KE = 0.5 * (M2 / Vol) * np.sum(u**2 * dV)
        return rs, Eth, KE

    def _closure(self, y):
        """Solve for q (and derived quantities) from the state vector."""
        R2, V2, R1, E2, M2 = y
        g, tab = self.g, self.tab

        def parts(q):
            xic, pic = tab.xic(q), tab.pic(q)
            Pc = pic * self.rho0 * V2**2
            Rc = xic * R2
            Vol = 4 * np.pi / 3 * max(Rc**3 - R1**3, 0.0)
            rs, Eth, KE = self._region2(R1, Rc, xic * V2, Pc, M2)
            return Pc, Rc, Vol, rs, Eth, KE

        def f(q):
            Pc, Rc, Vol, rs, Eth, KE = parts(q)
            return Eth + KE - E2

        qlo, qhi = tab.q[0], tab.q[-1]
        flo, fhi = f(qlo), f(qhi)
        if flo * fhi > 0:
            q = qlo if abs(flo) < abs(fhi) else qhi
            self._clipped = True
        else:
            q = brentq(f, qlo, qhi, xtol=1e-12, rtol=1e-12)
        Pc, Rc, Vol, rs, Eth, KE = parts(q)
        return dict(q=q, Pb=Pc, Rc=Rc, Vol=Vol, rs=rs, KE=KE, Eth=Eth,
                    xic=tab.xic(q), R2=R2, V2=V2, R1=R1, E2=E2, M2=M2)

    def _rhs(self, t, y):
        R2, V2, R1, E2, M2 = y
        c = self._closure(y)
        rs = c["rs"]
        Rcdot = c["xic"] * V2
        dE2 = rs["Eflux"] - c["Pb"] * 4 * np.pi * c["Rc"]**2 * Rcdot
        return [V2, c["q"] * V2**2 / R2, rs["R1dot"], dE2, rs["mflux"]]

    # ------------------------------------------------------------------ #
    def initial_state(self, t0):
        """Early wind-dominated state at time t0 (transient is forgotten)."""
        tab = self.tab
        q = 0.0
        V2 = 0.98 * self.vw
        R2 = V2 * t0
        xic = tab.xic(q)
        Rc = xic * R2
        R1 = 0.97 * Rc
        Pc = tab.pic(q) * self.rho0 * V2**2
        Vol = 4 * np.pi / 3 * (Rc**3 - R1**3)
        rs = self._reverse_shock(R1, Pc)
        M2 = self.Mdot * t0 - self.wind_mass(R1)
        M2 = max(M2, rs["rho2"] * Vol)
        rs, Eth, KE = self._region2(R1, Rc, xic * V2, Pc, M2)
        E2 = Eth + KE
        return np.array([R2, V2, R1, E2, M2])

    def wind_mass(self, R):
        r = np.linspace(0, R, 4001)
        rho, _, _ = self.wind.state(r)
        return np.trapz(4 * np.pi * r**2 * rho, r)

    def wind_energy(self, R):
        r = np.linspace(0, R, 4001)
        rho, u, P = self.wind.state(r)
        return np.trapz(4 * np.pi * r**2 * (0.5 * rho * u**2 + P / (self.g - 1)), r)

    def solve(self, t_end, t0=None, rtol=1e-8):
        """Integrate the bubble evolution from t0 to t_end [s]."""
        if t0 is None:
            # start when R2 ~ 100 r_inj (hypersonic wind at R1)
            t0 = 100 * self.r_inj / self.vw
        self._clipped = False
        y0 = self.initial_state(t0)
        self.t0 = t0
        self.sol = solve_ivp(self._rhs, (t0, t_end), y0, method="LSODA",
                             rtol=rtol, atol=[1e-6 * y0[0], 1e-6 * y0[1],
                                              1e-6 * y0[2], 1e-12 * self.L * t_end,
                                              1e-12 * self.Mdot * t_end],
                             dense_output=True)
        if not self.sol.success:
            raise RuntimeError(self.sol.message)
        self.t_end = t_end
        self._build_history()
        return self

    def _build_history(self, n=4000):
        """Record entropy and mass of shocked wind as a function of shock time."""
        ts = np.geomspace(self.t0, self.t_end, n)
        Y = self.sol.sol(ts)
        K = np.empty(n)
        for i in range(n):
            c = self._closure(Y[:, i])
            K[i] = c["rs"]["K"]
        self._hist_t = ts
        self._hist_K = K
        self._hist_M = Y[4]

    # ------------------------------------------------------------------ #
    def state(self, t):
        """Global quantities at time t [s]."""
        y = self.sol.sol(t)
        c = self._closure(y)
        g = self.g
        sh = self.tab.sols[int(np.argmin(abs(self.tab.q - c["q"])))]
        E_shell = self.rho0 * c["R2"]**3 * c["V2"]**2 * (sh.kin + sh.th)
        E_wind = self.wind_energy(c["R1"])
        c.update(t=t, E_shell=E_shell, E_wind=E_wind,
                 E_tot=E_shell + E_wind + c["E2"], E_inj=self.L * t,
                 M_wind=self.wind_mass(c["R1"]),
                 T_hot=U.temperature(c["Pb"], c["M2"] / c["Vol"], self.mu))
        return c

    def profile(self, r, t):
        """Fluid variables at radii r [cm] and time t [s].

        Returns dict with rho [g/cm^3], P [dyn/cm^2], T [K], u [cm/s].
        """
        r = np.asarray(r, dtype=float)
        c = self.state(t)
        g = self.g
        R1, Rc, R2, Pb, V2, q = c["R1"], c["Rc"], c["R2"], c["Pb"], c["V2"], c["q"]
        rho = np.empty_like(r)
        P = np.empty_like(r)
        u = np.empty_like(r)

        # I: free wind
        m = r < R1
        rho[m], u[m], P[m] = self.wind.state(r[m])

        # II: shocked wind -- Lagrangian reconstruction from entropy history
        m = (r >= R1) & (r < Rc)
        if np.any(m):
            rs = c["rs"]
            Rcdot = c["xic"] * V2
            sel = self._hist_t <= t
            Mh = np.concatenate([[0.0], self._hist_M[sel], [c["M2"]]])
            Kh = np.concatenate([[self._hist_K[sel][0]], self._hist_K[sel], [rs["K"]]])
            dm = np.diff(Mh)
            Kmid = 0.5 * (Kh[1:] + Kh[:-1])
            # earliest-shocked gas sits next to the contact; iterate because
            # the local pressure P(r) (Bernoulli) depends on position
            P_loc = np.full_like(Kmid, Pb)
            for _ in range(4):
                rho_p = (P_loc / Kmid) ** (1 / g)
                dV = dm / rho_p
                fV = c["Vol"] / max(dV.sum(), 1e-300)
                Vcum = np.concatenate([[0.0], np.cumsum(dV * fV)])
                r_edge = np.cbrt(np.maximum(Rc**3 - 3 * Vcum / (4 * np.pi), R1**3))
                r_mid = 0.5 * (r_edge[1:] + r_edge[:-1])
                P_loc, _ = self._bernoulli_P(r_mid, R1, Rc, rs, Rcdot)
            # rescale to the geometric volume (consistency of the closure)
            rho_p = rho_p / fV
            o = np.argsort(r_mid)
            rho[m] = np.interp(r[m], r_mid[o], rho_p[o])
            P[m], u[m] = self._bernoulli_P(r[m], R1, Rc, rs, Rcdot)
            c["vol_rescale"] = fV

        # III: shell
        m = (r >= Rc) & (r < R2)
        if np.any(m):
            G, Uu, Pi = self.tab.profile(q, r[m] / R2)
            rho[m] = self.rho0 * G
            u[m] = V2 * Uu
            P[m] = self.rho0 * V2**2 * Pi

        # IV: ambient
        m = r >= R2
        rho[m], u[m], P[m] = self.rho0, 0.0, self.P0

        T = U.temperature(P, rho, self.mu)
        return dict(rho=rho, P=P, T=T, u=u, state=c)

    def output_grid(self, t, n_uniform=4000, n_log=600, rmax_fac=1.5):
        """Radial grid [cm] from ~0 to rmax_fac*R2 with extra points near
        the injection region and on both sides of each discontinuity."""
        c = self.state(t)
        R1, Rc, R2 = c["R1"], c["Rc"], c["R2"]
        rmax = rmax_fac * R2
        r = [np.linspace(0, rmax, n_uniform),
             np.geomspace(1e-2 * self.r_inj, R1, n_log),
             np.linspace(R1, Rc, 400), np.linspace(Rc, R2, 400)]
        eps = 1e-6
        for Rd in (R1, Rc, R2):
            r.append([Rd * (1 - eps), Rd * (1 + eps)])
        r = np.unique(np.concatenate([np.ravel(a) for a in r]))
        return r[(r >= 0) & (r <= rmax)]


def weaver_adiabatic(Lw, nH, t, mu_H=U.MU_H):
    """Classic self-similar adiabatic Weaver et al. (1977) values (R1 -> 0).

    Returns R2 [cm], P_b [dyn cm^-2] using the self-similar shell constant.
    """
    tab = shell_table()
    q = -2.0 / 3.0
    xic, pic = tab.xic(q), tab.pic(q)
    a = 0.6
    A = (5 / (22 * np.pi * xic**3 * pic * a**2)) ** 0.2
    rho0 = U.rho_from_nH(nH, mu_H)
    R2 = A * (Lw * t**3 / rho0) ** 0.2
    V2 = a * R2 / t
    return R2, pic * rho0 * V2**2, A
