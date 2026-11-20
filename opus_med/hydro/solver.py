"""1D finite-volume Euler solver in spherical (or planar) geometry.

Conserved variables per unit volume U = (rho, rho v, E).  Spatial order 1
(piecewise constant) or 2 (piecewise-linear primitive reconstruction with an
MC limiter); time order 1 (forward Euler) or 2 (SSP-RK2 / Heun).  Fluxes are
computed with the HLLC approximate Riemann solver.  All operations are
vectorised over cells with numpy.
"""
import numpy as np

NG = 2  # ghost cells per side


class Grid:
    """Uniform grid of N cells on [r_min, r_max]."""

    def __init__(self, N, r_max, r_min=0.0, geometry="spherical"):
        self.N = N
        self.geometry = geometry
        self.rf = np.linspace(r_min, r_max, N + 1)          # faces
        self.dr = self.rf[1] - self.rf[0]
        self.rc = 0.5 * (self.rf[1:] + self.rf[:-1])        # centres
        if geometry == "spherical":
            self.area = self.rf**2                           # per steradian
            self.vol = (self.rf[1:]**3 - self.rf[:-1]**3) / 3.0
        elif geometry == "planar":
            self.area = np.ones_like(self.rf)
            self.vol = np.diff(self.rf)
        else:
            raise ValueError(geometry)
        self.dA = self.area[1:] - self.area[:-1]             # for P source


def prim_to_cons(rho, v, P, gamma):
    return np.array([rho, rho * v, P / (gamma - 1) + 0.5 * rho * v**2])


def cons_to_prim(U, gamma, rho_floor=1e-300, P_floor=1e-300):
    rho = np.maximum(U[0], rho_floor)
    v = U[1] / rho
    P = np.maximum((gamma - 1) * (U[2] - 0.5 * rho * v**2), P_floor)
    return rho, v, P


def hllc_flux(rhoL, vL, PL, rhoR, vR, PR, gamma):
    """HLLC flux (Toro 2009) with Davis/Einfeldt wave speed estimates."""
    cL = np.sqrt(gamma * PL / rhoL)
    cR = np.sqrt(gamma * PR / rhoR)
    SL = np.minimum(vL - cL, vR - cR)
    SR = np.maximum(vL + cL, vR + cR)
    EL = PL / (gamma - 1) + 0.5 * rhoL * vL**2
    ER = PR / (gamma - 1) + 0.5 * rhoR * vR**2
    num = PR - PL + rhoL * vL * (SL - vL) - rhoR * vR * (SR - vR)
    den = rhoL * (SL - vL) - rhoR * (SR - vR)
    Sm = num / den

    FL = np.array([rhoL * vL, rhoL * vL**2 + PL, (EL + PL) * vL])
    FR = np.array([rhoR * vR, rhoR * vR**2 + PR, (ER + PR) * vR])
    UL = np.array([rhoL, rhoL * vL, EL])
    UR = np.array([rhoR, rhoR * vR, ER])

    def star(rho, v, P, E, S):
        fac = rho * (S - v) / (S - Sm)
        return np.array([fac,
                         fac * Sm,
                         fac * (E / rho + (Sm - v) * (Sm + P / (rho * (S - v))))])

    UsL = star(rhoL, vL, PL, EL, SL)
    UsR = star(rhoR, vR, PR, ER, SR)
    F = np.where(SL >= 0, FL,
        np.where(Sm >= 0, FL + SL * (UsL - UL),
        np.where(SR > 0, FR + SR * (UsR - UR), FR)))
    return F


def mc_limiter(dL, dR):
    """Monotonised-central slope limiter."""
    s = np.where(dL * dR > 0, 1.0, 0.0) * np.sign(dL)
    return s * np.minimum(np.minimum(2 * np.abs(dL), 2 * np.abs(dR)),
                          0.5 * np.abs(dL + dR))


class Hydro:
    """Solver object.

    Parameters
    ----------
    grid : Grid
    gamma : adiabatic index
    order : 1 or 2 (same order in space and time)
    inner_bc, outer_bc : 'reflect' or 'outflow'
    source : callable(grid, t) -> (S_rho, S_mom, S_E) per unit volume per time,
             or None
    cfl : Courant number
    """

    def __init__(self, grid, gamma=5.0 / 3.0, order=2, inner_bc="reflect",
                 outer_bc="outflow", source=None, cfl=0.4,
                 rho_floor=1e-35, P_floor=1e-25):
        self.g = grid
        self.gamma = gamma
        self.order = order
        self.inner_bc = inner_bc
        self.outer_bc = outer_bc
        self.source = source
        self.cfl = cfl
        self.rho_floor = rho_floor
        self.P_floor = P_floor
        self.U = None
        self.t = 0.0

    # ------------------------------------------------------------------
    def set_initial(self, rho, v, P):
        self.U = prim_to_cons(np.asarray(rho, float) * np.ones(self.g.N),
                              np.asarray(v, float) * np.ones(self.g.N),
                              np.asarray(P, float) * np.ones(self.g.N),
                              self.gamma)

    def primitives(self, U=None):
        U = self.U if U is None else U
        return cons_to_prim(U, self.gamma, self.rho_floor, self.P_floor)

    # ------------------------------------------------------------------
    def _ghost(self, q, odd):
        """Pad a primitive array with NG ghost cells. odd=True flips sign at
        reflecting boundaries (velocity)."""
        left = q[:NG][::-1].copy()
        if self.inner_bc == "reflect":
            if odd:
                left = -left
        else:
            left[:] = q[0]
        right = q[-NG:][::-1].copy()
        if self.outer_bc == "reflect":
            if odd:
                right = -right
        else:
            right[:] = q[-1]
        return np.concatenate([left, q, right])

    def _rhs(self, U, t):
        g = self.g
        gam = self.gamma
        rho, v, P = cons_to_prim(U, gam, self.rho_floor, self.P_floor)
        W = [self._ghost(rho, False), self._ghost(v, True), self._ghost(P, False)]
        if self.order == 1:
            # faces i-1/2 for i = 0..N : left state cell NG-1+i, right NG+i
            WL = [w[NG - 1:NG + g.N] for w in W]
            WR = [w[NG:NG + g.N + 1] for w in W]
        else:
            WL, WR = [], []
            for k, w in enumerate(W):
                d = mc_limiter(w[1:-1] - w[:-2], w[2:] - w[1:-1])  # cells 1..-2
                wm = w[1:-1] - 0.5 * d   # left edge of cells 1..n-2
                wp = w[1:-1] + 0.5 * d   # right edge
                # face between padded cells j and j+1 (j = NG-1 .. NG+N-1)
                # in sliced arrays (offset 1): index j-1
                WL.append(wp[NG - 2:NG - 2 + g.N + 1])
                WR.append(wm[NG - 1:NG - 1 + g.N + 1])
            WL[0] = np.maximum(WL[0], self.rho_floor)
            WR[0] = np.maximum(WR[0], self.rho_floor)
            WL[2] = np.maximum(WL[2], self.P_floor)
            WR[2] = np.maximum(WR[2], self.P_floor)
        F = hllc_flux(WL[0], WL[1], WL[2], WR[0], WR[1], WR[2], gam)
        AF = F * g.area
        dUdt = -(AF[:, 1:] - AF[:, :-1]) / g.vol
        # geometric pressure source (exactly balances P * dA for uniform P)
        dUdt[1] += P * g.dA / g.vol
        if self.source is not None:
            S = self.source(g, t)
            dUdt += S
        return dUdt

    def timestep(self):
        rho, v, P = self.primitives()
        c = np.sqrt(self.gamma * P / rho)
        return self.cfl * self.g.dr / np.max(np.abs(v) + c)

    def _fix(self, U):
        rho, v, P = cons_to_prim(U, self.gamma, self.rho_floor, self.P_floor)
        return prim_to_cons(rho, v, P, self.gamma)

    def step(self, dt):
        U0, t0 = self.U, self.t
        if self.order == 1:
            U1 = U0 + dt * self._rhs(U0, t0)
            self.U = self._fix(U1)
        else:
            U1 = self._fix(U0 + dt * self._rhs(U0, t0))
            U2 = 0.5 * (U0 + U1 + dt * self._rhs(U1, t0 + dt))
            self.U = self._fix(U2)
        self.t = t0 + dt

    def evolve(self, t_end, callback=None, max_steps=10**8):
        """Advance to t_end exactly. Returns number of steps."""
        n = 0
        while self.t < t_end * (1 - 1e-14) and n < max_steps:
            dt = min(self.timestep(), t_end - self.t)
            self.step(dt)
            n += 1
            if callback is not None:
                callback(self)
        return n


def injection_source(R_inj, Mdot, L):
    """Uniform mass (Mdot) and energy (L) injection inside R_inj in spherical
    geometry; cells partially inside receive the volume-overlap fraction.
    Rates are per steradian internally (grid volumes are per steradian)."""
    V_inj = R_inj**3 / 3.0            # per steradian
    qm = Mdot / (4 * np.pi) / V_inj
    qe = L / (4 * np.pi) / V_inj
    cache = {}

    def S(grid, t):
        key = id(grid)
        if key not in cache:
            lo = np.minimum(grid.rf[:-1], R_inj)
            hi = np.minimum(grid.rf[1:], R_inj)
            frac = (hi**3 - lo**3) / 3.0 / grid.vol
            out = np.zeros((3, grid.N))
            out[0] = qm * frac
            out[2] = qe * frac
            cache[key] = out
        return cache[key]
    return S
