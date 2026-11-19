"""
Semi-analytic model of a spherically symmetric wind-blown bubble.

Structure (inside -> out)
-------------------------
  0  <  r < R1 : free-streaming, cold, supersonic wind
                  rho = Mdot / (4 pi r^2 v_w),  v = v_w
  R1           : reverse (wind termination) shock, strong, adiabatic
  R1 <  r < Rc : shocked wind ("hot bubble"); isobaric, P = P_b(t);
                  each fluid element evolves adiabatically after being shocked
  Rc           : contact discontinuity
  Rc <  r < R2 : swept-up ambient gas, thin isothermal (radiative) shell
  R2           : forward shock
  r  >  R2     : undisturbed uniform ambient medium

Global dynamics (thin-shell + isobaric-bubble approximation, as in
Castor, McCray & Weaver 1975; Weaver et al. 1977), *without* the Weaver
assumptions R1 << R2 and M_wind << M_swept:

  d(M_s V)/dt = 4 pi Rc^2 (P_b - P_0),          M_s = (4 pi / 3) rho_0 Rc^3
  dE_b/dt     = L_w (1 - dR1/dt / v_w) - 4 pi Rc^2 P_b V
  E_b         = P_b V_b / (gamma - 1),          V_b = (4 pi / 3)(Rc^3 - R1^3)
  P_b         = 2/(gamma+1) rho_w(R1) (v_w - dR1/dt)^2     (strong-shock jump)

The last relation is solved for dR1/dt, closing the system for
y = (Rc, V, R1, E_b).  The first term of dE_b/dt is the wind energy flux
through the *moving* reverse shock.

Interior profiles of the shocked wind are reconstructed from the solution
history in Lagrangian fashion: gas shocked at time t' has post-shock density
rho_s(t') = (gamma+1)/(gamma-1) rho_w(R1(t')) and pressure P_b(t'); it then
evolves adiabatically so that rho(t) = rho_s(t') [P_b(t)/P_b(t')]^(1/gamma).
Stacking those elements outward from R1 gives r(t') and therefore rho(r),
T(r); v(r) follows analytically from the continuity equation.

All internal quantities are CGS; see ``wind_bubble.units``.
"""
import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq

from . import units as u
from . import weaver

# Simple fit to the collisional-ionisation-equilibrium cooling function of
# solar-metallicity gas: Draine (2011, eq. 34.4) for 1e5 K < T < 2e7 K, joined
# continuously to a free-free (bremsstrahlung, Lambda ~ T^1/2) branch above.
# Used ONLY for a-posteriori diagnostics -- the dynamics are adiabatic.
def cooling_function(T):
    """Lambda(T) [erg cm^3 s^-1] such that the loss rate is n_H^2 Lambda."""
    T6 = np.clip(np.asarray(T, dtype=float) / 1e6, 0.1, None)
    line = 1.1e-22 * np.minimum(T6, 20.0) ** -0.7
    ff = 1.1e-22 * 20.0**-0.7 * (np.maximum(T6, 20.0) / 20.0) ** 0.5
    return np.where(T6 <= 20.0, line, ff)


@dataclass
class BubbleState:
    """Global (radius-independent) quantities at one time, CGS."""
    t: float
    R1: float           # reverse-shock radius
    Rc: float           # contact discontinuity (inner edge of shell)
    R2: float           # forward shock (outer edge of shell)
    V: float            # shell velocity dRc/dt
    dR1dt: float        # reverse-shock velocity
    P_b: float          # hot-bubble pressure
    E_b: float          # hot-bubble thermal energy
    M_b: float          # mass of shocked wind
    M_s: float          # mass of swept-up shell
    P_shell_out: float  # pressure just behind the forward shock
    extra: dict = field(default_factory=dict)


class WindBubble:
    """Semi-analytic, spherically symmetric wind-blown bubble.

    Parameters
    ----------
    L_w : float
        Mechanical luminosity of the wind [erg/s].
    Mdot_w : float
        Mass-loss rate [g/s].
    n_H : float
        Ambient hydrogen-nucleus number density [cm^-3].
    mu_H : float
        Mean mass per H nucleus in units of m_p (rho = mu_H m_p n_H).
    T_ambient : float
        Ambient temperature [K] (sets ambient pressure P_0).
    ambient_state : str
        "atomic", "molecular" or "ionized" -- sets mean particle mass of the
        ambient gas for computing P_0 and T.
    T_shell : float
        Temperature of the radiatively cooled, swept-up shell [K]; the shell
        gas is taken to be (photo)ionised.
    T_wind : float
        Temperature of the free-streaming wind [K] (only used to report a
        pressure/temperature in region r < R1; dynamically the wind is cold).
    gamma : float
        Adiabatic index of the hot gas.
    """

    def __init__(self, L_w, Mdot_w, n_H, mu_H=1.4, T_ambient=100.0,
                 ambient_state="atomic", T_shell=1.0e4, T_wind=1.0e4,
                 gamma=5.0 / 3.0):
        self.L_w = float(L_w)
        self.Mdot_w = float(Mdot_w)
        self.n_H = float(n_H)
        self.mu_H = float(mu_H)
        self.gamma = float(gamma)
        self.T_ambient = float(T_ambient)
        self.T_shell = float(T_shell)
        self.T_wind = float(T_wind)

        self.rho0 = self.mu_H * u.m_p * self.n_H
        self.mu_amb = u.mean_mass_per_particle(self.mu_H, state=ambient_state)
        self.mu_ion = u.mean_mass_per_particle(self.mu_H, state="ionized")
        self.P0 = self.rho0 * u.k_B * self.T_ambient / (self.mu_amb * u.m_p)
        self.cs2_shell = u.k_B * self.T_shell / (self.mu_ion * u.m_p)

        # terminal wind velocity (all injected energy is kinetic)
        self.v_w = np.sqrt(2.0 * self.L_w / self.Mdot_w)

        # characteristic scales: swept-up mass = wind mass at t_sw
        #   (4 pi/3) rho0 (v_w t)^3 = Mdot t
        self.t_sw = np.sqrt(3.0 * self.Mdot_w / (4.0 * np.pi * self.rho0 * self.v_w**3))
        self.l_sc = self.v_w * self.t_sw
        self.m_sc = self.Mdot_w * self.t_sw
        self.E_sc = self.m_sc * self.v_w**2
        self.P_sc = self.m_sc / (self.l_sc * self.t_sw**2)

        self._sol = None

    # ------------------------------------------------------------------ ODEs
    # Dimensionless variables: v_w = Mdot = 1, L = 1/2, rho0 = 3/(4 pi).
    def _shock_speed(self, x):
        """Solve the reverse-shock condition for w = v_w - dR1/dt.

        x = P_b / rho_w(R1)  (dimensionless, v_w = 1).  Energy-consistent
        condition (per unit shocked mass, see module docstring / assumptions):

            1/2 = gamma/(gamma-1) P_b/rho_eff + P_b (1 - w)/(rho_w w)

        with rho_eff = rho_s (P_b/P_s)^(1/gamma), rho_s = (g+1)/(g-1) rho_w,
        P_s = 2/(g+1) rho_w w^2.  f(w) is convex with f(0+) = +inf and
        f(w_max) <= 0 at w_max = sqrt((g+1) x / 2) (P_s = P_b), so the root is
        unique; found by brentq (scalar) or vectorised bisection in log w.
        """
        g = self.gamma
        x = np.asarray(x, dtype=float)

        def f(w):
            Ps_over_Pb = 2.0 / (g + 1.0) * w**2 / x
            term1 = g / (g + 1.0) * x * Ps_over_Pb ** (1.0 / g)
            return term1 + x * (1.0 - w) / w - 0.5

        hi = np.log(np.sqrt(0.5 * (g + 1.0) * x))
        if x.ndim == 0:                       # fast scalar path (ODE RHS)
            lw = brentq(lambda lw: f(np.exp(lw)), float(hi) - 60.0, float(hi),
                        xtol=1e-14, rtol=1e-15)
            return np.exp(lw)
        lo = hi - 60.0
        for _ in range(64):
            mid = 0.5 * (lo + hi)
            pos = f(np.exp(mid)) > 0
            lo = np.where(pos, mid, lo)
            hi = np.where(pos, hi, mid)
        return np.exp(0.5 * (lo + hi))

    def _derived(self, y):
        Rc, V, R1, E = y
        g = self.gamma
        Vb = 4.0 * np.pi / 3.0 * (Rc**3 - R1**3)
        P = (g - 1.0) * E / Vb
        rho_w = 1.0 / (4.0 * np.pi * R1**2)
        w = self._shock_speed(P / rho_w)
        dR1 = 1.0 - w
        return P, Vb, dR1

    def _rho_eff(self, P, R1, dR1):
        """Density of freshly shocked wind after isentropic adjustment to P_b."""
        g = self.gamma
        rho_w = 1.0 / (4.0 * np.pi * R1**2)
        w = 1.0 - dR1
        rho_s = (g + 1.0) / (g - 1.0) * rho_w
        P_s = 2.0 / (g + 1.0) * rho_w * w**2
        return rho_s * (P / P_s) ** (1.0 / g)

    def _rhs(self, t, y):
        Rc, V, R1, E = y
        P, Vb, dR1 = self._derived(y)
        P0 = self.P0 / self.P_sc
        Ms = Rc**3                                   # (4pi/3) rho0 Rc^3, rho0=3/4pi
        dMs = 3.0 * Rc**2 * V
        dV = (4.0 * np.pi * Rc**2 * (P - P0) - dMs * V) / Ms
        dE = 0.5 * (1.0 - dR1) - 4.0 * np.pi * Rc**2 * P * V
        return [V, dV, dR1, dE]

    def solve(self, t_end, t_start_frac=1e-4, R1_frac0=0.9, rtol=1e-9,
              atol=1e-14):
        """Integrate the global evolution from ~0 to ``t_end`` [s].

        The integration starts at t0 = t_start_frac * t_sw in the
        free-expansion phase with Rc = v_w t0, V = v_w, R1 = R1_frac0 * Rc and
        a bubble pressure chosen so that dR1/dt = R1_frac0 v_w (self-similar
        early phase).  The memory of these initial conditions decays rapidly
        (checked in tests/).
        """
        g = self.gamma
        t0 = t_start_frac
        Rc0 = t0
        R10 = R1_frac0 * Rc0
        # P_b such that dR1/dt = R1_frac0, i.e. w = 1 - R1_frac0
        rho_w0 = 1.0 / (4.0 * np.pi * R10**2)
        w0 = 1.0 - R1_frac0
        x = brentq(lambda lx: np.log(self._shock_speed(np.exp(lx)) / w0), -60, 60)
        P0b = np.exp(x) * rho_w0
        E0 = P0b * 4.0 * np.pi / 3.0 * (Rc0**3 - R10**3) / (g - 1.0)
        y0 = [Rc0, 1.0, R10, E0]
        self._t0 = t0
        self._M0_bubble = t0 - R10           # wind mass injected but not in free wind
        self._P0_bubble = P0b
        self._V0_bubble = 4.0 * np.pi / 3.0 * (Rc0**3 - R10**3)

        tau_end = t_end / self.t_sw
        sol = solve_ivp(self._rhs, (t0, tau_end), y0, method="Radau",
                        rtol=rtol, atol=atol, dense_output=True)
        if not sol.success:
            raise RuntimeError(f"ODE integration failed: {sol.message}")
        self._sol = sol
        self._tau_end = tau_end
        return self

    # --------------------------------------------------------------- helpers
    def _check(self, t):
        if self._sol is None:
            raise RuntimeError("call .solve(t_end) first")
        tau = np.asarray(t) / self.t_sw
        if np.any(tau > self._tau_end * (1 + 1e-12)) or np.any(tau < self._t0):
            raise ValueError("requested time outside solved range")
        return tau

    def _shell_structure(self, Rc, V, P_b):
        """Isothermal shell between Rc and R2 (CGS in, CGS out).

        Pressure varies linearly from P_b at Rc to the post-shock value at R2;
        R2 is fixed by requiring the shell to contain all swept-up mass
        (4 pi/3) rho0 R2^3.
        """
        c2 = self.cs2_shell
        # isothermal shock jump: P_out = P0 + rho0 V^2 (1 - rho0/rho_out), rho_out = P_out/c2
        P_out = self.P0 + self.rho0 * V**2
        for _ in range(50):
            rho_out = max(P_out / c2, self.rho0)
            P_new = self.P0 + self.rho0 * V**2 * (1.0 - self.rho0 / rho_out)
            if abs(P_new - P_out) < 1e-12 * P_out:
                break
            P_out = P_new
        P_out = max(P_out, self.P0)

        def mass_excess(R2):
            x, w = np.polynomial.legendre.leggauss(16)
            r = 0.5 * (R2 - Rc) * x + 0.5 * (R2 + Rc)
            P = P_b + (P_out - P_b) * (r - Rc) / (R2 - Rc)
            m_shell = 0.5 * (R2 - Rc) * np.sum(w * 4 * np.pi * r**2 * P / c2)
            return m_shell - 4.0 * np.pi / 3.0 * self.rho0 * R2**3

        rho_min = min(P_b, P_out) / c2
        if rho_min <= self.rho0:
            warnings.warn("isothermal shell would not be denser than the ambient "
                          "medium (P_b < rho0 c_shell^2): thin-shell structure "
                          "invalid; R2 set to NaN", RuntimeWarning)
            return np.nan, np.nan
        R2 = brentq(mass_excess, Rc * (1 + 1e-12), Rc * 10.0, xtol=1e-14 * Rc)
        return R2, P_out

    # ------------------------------------------------------- global quantities
    def state(self, t):
        """Global bubble properties at time ``t`` [s] as a :class:`BubbleState`."""
        tau = float(self._check(t))
        y = self._sol.sol(tau)
        Rc, V, R1, E = y
        P, Vb, dR1 = self._derived(y)
        Rc_cgs, V_cgs = Rc * self.l_sc, V * self.v_w
        P_cgs = P * self.P_sc
        R2, P_out = self._shell_structure(Rc_cgs, V_cgs, P_cgs)
        return BubbleState(
            t=t, R1=R1 * self.l_sc, Rc=Rc_cgs, R2=R2, V=V_cgs,
            dR1dt=dR1 * self.v_w, P_b=P_cgs, E_b=E * self.E_sc,
            M_b=(tau - R1) * self.m_sc,
            M_s=4.0 * np.pi / 3.0 * self.rho0 * R2**3,
            P_shell_out=P_out,
        )

    def history(self, t):
        """Vectorised global quantities (CGS) at times ``t`` (array)."""
        tau = self._check(t)
        y = self._sol.sol(tau)
        Rc, V, R1, E = y
        P, Vb, dR1 = self._derived(y)
        return dict(t=np.asarray(t), R1=R1 * self.l_sc, Rc=Rc * self.l_sc,
                    V=V * self.v_w, dR1dt=dR1 * self.v_w, P_b=P * self.P_sc,
                    E_b=E * self.E_sc, M_b=(tau - R1) * self.m_sc)

    # ------------------------------------------------------- radial profiles
    def _shocked_wind_lagrangian(self, tau, n_elem=20000):
        """Lagrangian reconstruction of the shocked wind at time tau.

        Returns (dimensionless) arrays ordered from R1 outward: radius r,
        density rho, velocity v, together with the leftover volume belonging
        to the gas that was in the bubble at t0 (initial condition gas).
        """
        g = self.gamma
        sol = self._sol
        # shock times s from tau back to t0, log-spaced (dense near tau as well)
        s_log = np.geomspace(self._t0, tau, n_elem)
        s_lin = tau - np.geomspace(1e-10 * tau, tau - self._t0, n_elem)[::-1]
        s = np.unique(np.concatenate([s_log, s_lin, [tau]]))
        s = s[(s >= self._t0) & (s <= tau)]
        ys = sol.sol(s)
        P_s, _, dR1_s = self._derived(ys)
        R1_s = ys[2]
        mdot_sh = 1.0 - dR1_s                                # mass flux through shock
        rho_post = self._rho_eff(P_s, R1_s, dR1_s)           # at pressure P_b(s)

        yt = sol.sol(tau)
        Pt, Vbt, dR1t = self._derived(yt)
        Rc_t, V_t, R1_t, E_t = yt
        rho_now = rho_post * (Pt / P_s) ** (1.0 / g)

        # volume of gas shocked between s and tau (s descending -> outward)
        integrand = mdot_sh / rho_now
        cum = cumulative_trapezoid(integrand[::-1], s[::-1], initial=0.0)
        vol = -cum                          # s decreasing => positive volume
        r = (R1_t**3 + 3.0 / (4.0 * np.pi) * vol) ** (1.0 / 3.0)
        rho = rho_now[::-1]

        # Eulerian velocity from continuity (exact for isobaric adiabatic flow)
        dE = 0.5 * (1.0 - dR1t) - 4.0 * np.pi * Rc_t**2 * Pt * V_t
        dVb = 4.0 * np.pi * (Rc_t**2 * V_t - R1_t**2 * dR1t)
        dlnP = dE / E_t - dVb / Vbt
        rho_in = self._rho_eff(Pt, R1_t, dR1t)
        def vel(rr):
            return (R1_t**2 * dR1t + (1.0 - dR1t) / (4.0 * np.pi * rho_in)
                    - dlnP / (3.0 * g) * (rr**3 - R1_t**3)) / rr**2
        v = vel(r)

        # gas from the initial condition fills the rest of the bubble
        V_ic_expected = self._V0_bubble * (self._P0_bubble / Pt) ** (1.0 / g)
        V_ic_residual = Vbt - 4.0 * np.pi / 3.0 * (r[-1]**3 - R1_t**3)
        return dict(r=r, rho=rho, v=v, vel=vel, P=Pt,
                    V_ic_expected=V_ic_expected, V_ic_residual=V_ic_residual,
                    V_b=Vbt, R1=R1_t, Rc=Rc_t)

    def profiles(self, t, r=None, n_r=2000, r_max_factor=1.5):
        """Radial profiles of the fluid variables at time ``t`` [s].

        Parameters
        ----------
        t : float
            Time since the wind turned on [s].
        r : array_like, optional
            Radii [cm].  Default: ``n_r`` points from 0 to
            ``r_max_factor * R2``, with the discontinuities R1, Rc, R2 inserted
            explicitly (each appears twice, once with the inner and once with
            the outer limiting value).

        Returns
        -------
        dict with CGS arrays ``r, rho, n_H, P, T, v`` and an integer array
        ``region`` (0 free wind, 1 shocked wind, 2 shell, 3 ambient), plus the
        :class:`BubbleState` under key ``state`` and consistency diagnostics.
        """
        tau = float(self._check(t))
        st = self.state(t)
        if not np.isfinite(st.R2):
            raise ValueError("thin-shell structure invalid at this time; cannot "
                             "build profiles (see warning)")
        lag = self._shocked_wind_lagrangian(tau)
        l = self.l_sc

        if r is None:
            r = np.linspace(0.0, r_max_factor * st.R2, n_r)[1:]
            eps = 1e-9
            special = []
            for R in (st.R1, st.Rc, st.R2):
                special += [R * (1 - eps), R * (1 + eps)]
            r = np.sort(np.concatenate([r, special]))
        r = np.asarray(r, dtype=float)

        rho = np.empty_like(r)
        P = np.empty_like(r)
        T = np.empty_like(r)
        v = np.empty_like(r)
        region = np.empty(r.shape, dtype=int)

        # region 0: free wind
        m0 = r < st.R1
        rho[m0] = self.Mdot_w / (4.0 * np.pi * r[m0]**2 * self.v_w)
        v[m0] = self.v_w
        T[m0] = self.T_wind
        P[m0] = rho[m0] * u.k_B * self.T_wind / (self.mu_ion * u.m_p)
        region[m0] = 0

        # region 1: shocked wind
        m1 = (r >= st.R1) & (r < st.Rc)
        rr = r[m1] / l
        r_lag, rho_lag = lag["r"], lag["rho"]
        # gas from initial conditions: uniform density filling remaining volume
        V_ic = max(lag["V_ic_residual"], 1e-300)
        rho_ic = self._M0_bubble / V_ic
        rho1 = np.where(rr <= r_lag[-1],
                        np.interp(rr, r_lag, rho_lag),
                        rho_ic)
        rho[m1] = rho1 * self.m_sc / l**3
        P[m1] = st.P_b
        v[m1] = lag["vel"](rr) * self.v_w
        T[m1] = st.P_b * self.mu_ion * u.m_p / (rho[m1] * u.k_B)
        region[m1] = 1

        # region 2: isothermal shell, linear pressure Rc -> R2
        m2 = (r >= st.Rc) & (r < st.R2)
        f = (r[m2] - st.Rc) / (st.R2 - st.Rc)
        P[m2] = st.P_b + (st.P_shell_out - st.P_b) * f
        rho[m2] = P[m2] / self.cs2_shell
        rho_out = st.P_shell_out / self.cs2_shell
        u_out = st.V * (1.0 - self.rho0 / rho_out)
        v[m2] = st.V + (u_out - st.V) * f
        T[m2] = self.T_shell
        region[m2] = 2

        # region 3: ambient medium
        m3 = r >= st.R2
        rho[m3] = self.rho0
        P[m3] = self.P0
        v[m3] = 0.0
        T[m3] = self.T_ambient
        region[m3] = 3

        n_H = rho / (self.mu_H * u.m_p)

        # ---- diagnostics
        diag = {
            # fraction of bubble volume whose gas came from the initial condition
            "V_ic_fraction": lag["V_ic_residual"] / lag["V_b"],
            # independent adiabatic prediction for that fraction
            "V_ic_fraction_expected": lag["V_ic_expected"] / lag["V_b"],
            # mismatch between Lagrangian velocity at Rc and shell velocity
            "v_at_Rc_over_V": lag["vel"](lag["Rc"]) * self.v_w / st.V,
        }
        # radiative-loss diagnostic (adiabatic dynamics assumed!)
        rs = r_lag * l
        n_lag = rho_lag * self.m_sc / l**3 / (self.mu_H * u.m_p)
        T_lag = st.P_b * self.mu_ion * u.m_p / (rho_lag * self.m_sc / l**3 * u.k_B)
        hot = T_lag >= 1e5                          # validity range of the fit
        integrand = np.where(hot, 4 * np.pi * rs**2 * n_lag**2 * cooling_function(T_lag), 0.0)
        Lcool = np.trapezoid(integrand, rs)
        dm = 4 * np.pi * rs**2 * n_lag * self.mu_H * u.m_p
        diag["L_cool_bubble"] = Lcool
        diag["L_cool_over_L_w"] = Lcool / self.L_w
        diag["shocked_wind_mass_fraction_below_1e5K"] = (
            np.trapezoid(np.where(T_lag < 1e5, dm, 0.0), rs) / np.trapezoid(dm, rs))
        n1, T1 = n_lag[0], T_lag[0]
        diag["t_cool_postshock"] = 1.5 * (n1 * self.mu_H / self.mu_ion) * u.k_B * T1 / (
            n1**2 * cooling_function(T1))
        diag["T_postshock"] = T1

        return dict(r=r, rho=rho, n_H=n_H, P=P, T=T, v=v, region=region,
                    state=st, diagnostics=diag)

    # --------------------------------------------------------------- summary
    def weaver_comparison(self, t):
        """Ratios of this model to the Weaver et al. (1977) self-similar values."""
        st = self.state(t)
        # Weaver's thin-shell radius corresponds to Rc here (R2 adds the
        # finite thickness of the isothermal shell)
        return dict(
            Rc=st.Rc / weaver.outer_radius(self.L_w, self.rho0, t),
            R2=st.R2 / weaver.outer_radius(self.L_w, self.rho0, t),
            V=st.V / weaver.shell_velocity(self.L_w, self.rho0, t),
            P_b=st.P_b / weaver.bubble_pressure(self.L_w, self.rho0, t),
            E_b=st.E_b / weaver.bubble_thermal_energy(self.L_w, t),
            R1=st.R1 / weaver.reverse_shock_radius(self.L_w, self.Mdot_w, self.rho0, t,
                                                   self.gamma),
        )
