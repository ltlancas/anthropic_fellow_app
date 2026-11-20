"""1D finite-volume Euler solver in spherical (or Cartesian) geometry.

Conserved variables U = (rho, rho v, E), E = P/(gamma-1) + rho v^2/2.
  order=1 : piecewise-constant reconstruction + forward Euler.
  order=2 : piecewise-linear (MC-limited, primitive variables) reconstruction + SSP-RK2 (Heun).
Riemann solver: HLLC. Optional volumetric mass/energy source (wind injection).
"""
import numpy as np
from numba import njit

FLOOR_RHO = 1e-40
FLOOR_P = 1e-40
SRC_FRAC = 0.3


@njit(cache=True)
def _hllc(rL, vL, pL, rR, vR, pR, g):
    cL = np.sqrt(g * pL / rL)
    cR = np.sqrt(g * pR / rR)
    EL = pL / (g - 1.0) + 0.5 * rL * vL * vL
    ER = pR / (g - 1.0) + 0.5 * rR * vR * vR
    SL = min(vL - cL, vR - cR)
    SR = max(vL + cL, vR + cR)
    fL0 = rL * vL; fL1 = rL * vL * vL + pL; fL2 = (EL + pL) * vL
    fR0 = rR * vR; fR1 = rR * vR * vR + pR; fR2 = (ER + pR) * vR
    if SL >= 0.0:
        return fL0, fL1, fL2
    if SR <= 0.0:
        return fR0, fR1, fR2
    SM = (pR - pL + rL * vL * (SL - vL) - rR * vR * (SR - vR)) / (rL * (SL - vL) - rR * (SR - vR))
    if SM >= 0.0:
        fac = rL * (SL - vL) / (SL - SM)
        u0 = fac; u1 = fac * SM
        u2 = fac * (EL / rL + (SM - vL) * (SM + pL / (rL * (SL - vL))))
        return fL0 + SL * (u0 - rL), fL1 + SL * (u1 - rL * vL), fL2 + SL * (u2 - EL)
    fac = rR * (SR - vR) / (SR - SM)
    u0 = fac; u1 = fac * SM
    u2 = fac * (ER / rR + (SM - vR) * (SM + pR / (rR * (SR - vR))))
    return fR0 + SR * (u0 - rR), fR1 + SR * (u1 - rR * vR), fR2 + SR * (u2 - ER)


@njit(cache=True)
def _mc(a, b):
    if a * b <= 0.0:
        return 0.0
    s = 1.0 if a > 0 else -1.0
    return s * min(2.0 * abs(a), 2.0 * abs(b), 0.5 * abs(a + b))


@njit(cache=True)
def cons_to_prim(U, g, W):
    n = U.shape[1]
    for i in range(n):
        rho = max(U[0, i], FLOOR_RHO)
        v = U[1, i] / rho
        p = (g - 1.0) * (U[2, i] - 0.5 * rho * v * v)
        if p < FLOOR_P:
            p = FLOOR_P
            U[2, i] = p / (g - 1.0) + 0.5 * rho * v * v
        U[0, i] = rho
        W[0, i] = rho; W[1, i] = v; W[2, i] = p


@njit(cache=True)
def prim_to_cons(W, g):
    U = np.empty_like(W)
    U[0] = W[0]
    U[1] = W[0] * W[1]
    U[2] = W[2] / (g - 1.0) + 0.5 * W[0] * W[1] ** 2
    return U


@njit(cache=True)
def rhs(U, area, vol, darea, g, order, inner_reflect, src_rho, src_E, dU):
    """Compute dU/dt; returns maximum signal speed."""
    n = U.shape[1]
    W = np.empty((3, n + 4))
    Wi = np.empty((3, n))
    cons_to_prim(U, g, Wi)
    for k in range(3):
        for i in range(n):
            W[k, i + 2] = Wi[k, i]
        # outer: zero-gradient outflow
        W[k, n + 2] = Wi[k, n - 1]
        W[k, n + 3] = Wi[k, n - 1]
        # inner: reflecting (spherical centre) or outflow
        W[k, 1] = Wi[k, 0]
        W[k, 0] = Wi[k, 1]
    if inner_reflect:
        W[1, 1] = -Wi[1, 0]
        W[1, 0] = -Wi[1, 1]
    # slopes (cells 1..n+2 in padded indexing)
    S = np.zeros((3, n + 4))
    if order == 2:
        for i in range(1, n + 3):
            for k in range(3):
                S[k, i] = _mc(W[k, i] - W[k, i - 1], W[k, i + 1] - W[k, i])
            # positivity fallback
            if (W[0, i] - 0.5 * abs(S[0, i]) <= 0.0) or (W[2, i] - 0.5 * abs(S[2, i]) <= 0.0):
                S[0, i] = 0.0; S[1, i] = 0.0; S[2, i] = 0.0
    F = np.empty((3, n + 1))
    smax = 0.0
    for f in range(n + 1):
        iL = f + 1
        iR = f + 2
        rL = W[0, iL] + 0.5 * S[0, iL]; vL = W[1, iL] + 0.5 * S[1, iL]; pL = W[2, iL] + 0.5 * S[2, iL]
        rR = W[0, iR] - 0.5 * S[0, iR]; vR = W[1, iR] - 0.5 * S[1, iR]; pR = W[2, iR] - 0.5 * S[2, iR]
        a, b, c = _hllc(rL, vL, pL, rR, vR, pR, g)
        F[0, f] = a; F[1, f] = b; F[2, f] = c
    for i in range(n):
        dU[0, i] = -(area[i + 1] * F[0, i + 1] - area[i] * F[0, i]) / vol[i] + src_rho[i]
        dU[1, i] = -(area[i + 1] * F[1, i + 1] - area[i] * F[1, i]) / vol[i] + Wi[2, i] * darea[i] / vol[i]
        dU[2, i] = -(area[i + 1] * F[2, i + 1] - area[i] * F[2, i]) / vol[i] + src_E[i]
        s = abs(Wi[1, i]) + np.sqrt(g * Wi[2, i] / Wi[0, i])
        if s > smax:
            smax = s
    return smax


@njit(cache=True)
def advance(U, t, t_target, dx_min, area, vol, darea, g, order, inner_reflect,
            src_rho, src_E, cfl, max_steps):
    """Advance U in place from t to t_target. Returns (t, nsteps)."""
    n = U.shape[1]
    dU = np.empty_like(U)
    U1 = np.empty_like(U)
    nsteps = 0
    while t < t_target and nsteps < max_steps:
        smax = rhs(U, area, vol, darea, g, order, inner_reflect, src_rho, src_E, dU)
        dt = cfl * dx_min / smax
        # source-term limiter: injected mass/energy per step <= SRC_FRAC of cell content
        for i in range(n):
            if src_E[i] > 0.0:
                dt = min(dt, SRC_FRAC * U[2, i] / src_E[i])
            if src_rho[i] > 0.0:
                dt = min(dt, SRC_FRAC * U[0, i] / src_rho[i])
        if t + dt > t_target:
            dt = t_target - t
        if order == 1:
            for k in range(3):
                for i in range(n):
                    U[k, i] += dt * dU[k, i]
        else:
            for k in range(3):
                for i in range(n):
                    U1[k, i] = U[k, i] + dt * dU[k, i]
            rhs(U1, area, vol, darea, g, order, inner_reflect, src_rho, src_E, dU)
            for k in range(3):
                for i in range(n):
                    U[k, i] = 0.5 * (U[k, i] + U1[k, i] + dt * dU[k, i])
        # enforce floors
        for i in range(n):
            if U[0, i] < FLOOR_RHO:
                U[0, i] = FLOOR_RHO
        t += dt
        nsteps += 1
    return t, nsteps


class Grid:
    """Uniform grid on [r_min, r_max] with N cells. geometry: 'spherical' or 'cartesian'."""

    def __init__(self, N, r_max, r_min=0.0, geometry="spherical"):
        self.N = N
        self.geometry = geometry
        self.rf = np.linspace(r_min, r_max, N + 1)
        self.dr = np.diff(self.rf)
        if geometry == "spherical":
            self.area = 4 * np.pi * self.rf ** 2
            self.vol = 4 * np.pi / 3 * (self.rf[1:] ** 3 - self.rf[:-1] ** 3)
            # volume-weighted cell centre
            self.rc = 0.75 * (self.rf[1:] ** 4 - self.rf[:-1] ** 4) / (self.rf[1:] ** 3 - self.rf[:-1] ** 3)
        elif geometry == "cartesian":
            self.area = np.ones(N + 1)
            self.vol = self.dr.copy()
            self.rc = 0.5 * (self.rf[1:] + self.rf[:-1])
        else:
            raise ValueError(geometry)
        self.darea = self.area[1:] - self.area[:-1]

    def sphere_overlap(self, R):
        """Fraction of each cell's volume lying inside radius R (spherical)."""
        lo = np.minimum(self.rf[:-1], R)
        hi = np.minimum(self.rf[1:], R)
        return (hi ** 3 - lo ** 3) / (self.rf[1:] ** 3 - self.rf[:-1] ** 3)


class Simulation:
    """Container for a run.

    W0: initial primitive array (3, N) [rho, v, P]. Sources: mdot [g/s] and L [erg/s] deposited
    uniformly within radius R_inj (spherical only).
    """

    def __init__(self, grid, W0, gamma=5.0 / 3.0, order=2, cfl=0.4,
                 mdot=0.0, L=0.0, R_inj=0.0):
        self.grid, self.g, self.order, self.cfl = grid, gamma, order, cfl
        self.U = prim_to_cons(np.ascontiguousarray(W0, dtype=float), gamma)
        self.t = 0.0
        self.nsteps = 0
        self.src_rho = np.zeros(grid.N)
        self.src_E = np.zeros(grid.N)
        if (mdot > 0 or L > 0) and R_inj > 0:
            frac = grid.sphere_overlap(R_inj)
            Vinj = 4 * np.pi / 3 * R_inj ** 3
            # rate per unit volume, distributed by overlapping volume
            self.src_rho = mdot / Vinj * frac
            self.src_E = L / Vinj * frac
        self.inner_reflect = grid.geometry == "spherical"

    def run_to(self, t_target, max_steps=10 ** 9):
        g = self.grid
        t, n = advance(self.U, self.t, t_target, float(g.dr.min()), g.area, g.vol, g.darea,
                       self.g, self.order, self.inner_reflect, self.src_rho, self.src_E,
                       self.cfl, max_steps)
        self.t = t
        self.nsteps += n
        return t

    def primitives(self):
        W = np.empty_like(self.U)
        cons_to_prim(self.U, self.g, W)
        return W
