"""Finite-volume solver for the 1D Euler equations with source terms.

Equations (geometry factor alpha = 2 for spherical, 0 for planar):

  d/dt (rho)  + 1/r^a d/dr (r^a rho u)        = S_rho
  d/dt (rho u)+ 1/r^a d/dr (r^a (rho u^2)) + dP/dr = 0
  d/dt (E)    + 1/r^a d/dr (r^a (E+P) u)      = S_E

Discretisation:  dU_i/dt = -(A_{i+1/2} F_{i+1/2} - A_{i-1/2} F_{i-1/2}) / V_i
                           + [0, P_i (A_{i+1/2}-A_{i-1/2})/V_i, 0] + S_i
with A = r^a, V = (r_{i+1/2}^{a+1} - r_{i-1/2}^{a+1})/(a+1).  Writing the
momentum flux as rho u^2 + P and adding the geometric source P dA/V makes a
uniform pressure an exact equilibrium.

order = 1 : piecewise-constant reconstruction + forward Euler.
order = 2 : piecewise-linear reconstruction of primitive variables with a
            limited slope (MC or minmod) + 2-stage SSP Runge-Kutta (Heun).
"""
import numpy as np

from .riemann import hllc

NG = 2  # ghost cells per side


def _minmod(a, b):
    return np.where(a * b > 0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)


def _mc(a, b):
    return np.where(a * b > 0,
                    np.sign(a) * np.minimum(np.minimum(2 * np.abs(a), 2 * np.abs(b)),
                                            0.5 * np.abs(a + b)), 0.0)


class Hydro1D:
    """1D hydro solver.

    Parameters
    ----------
    r_in, r_out : domain boundaries (code units)
    n           : number of cells
    order       : 1 or 2
    geometry    : 'spherical' or 'planar'
    gamma       : adiabatic index
    cfl         : Courant number
    limiter     : 'mc' or 'minmod' (order 2 only)
    inner_bc, outer_bc : 'reflect', 'outflow', or a callable
                  f(r_ghost, t) -> (rho, u, P) giving the ghost-cell state
    source      : optional callable f(solver, U, t) -> dU/dt array (3, n)
    rho_floor, P_floor : positivity floors (code units)
    """

    def __init__(self, r_in, r_out, n, order=2, geometry="spherical",
                 gamma=5.0 / 3.0, cfl=0.4, limiter="mc",
                 inner_bc="reflect", outer_bc="outflow", source=None,
                 rho_floor=1e-12, P_floor=1e-20):
        self.n, self.order, self.g, self.cfl = n, order, gamma, cfl
        self.alpha = 2 if geometry == "spherical" else 0
        self.limiter = _mc if limiter == "mc" else _minmod
        self.inner_bc, self.outer_bc, self.source = inner_bc, outer_bc, source
        self.rho_floor, self.P_floor = rho_floor, P_floor

        self.rf = np.linspace(r_in, r_out, n + 1)          # faces
        self.dr = np.diff(self.rf)
        self.r = 0.5 * (self.rf[1:] + self.rf[:-1])        # centres
        dr0 = self.dr[0]
        self.rg = np.concatenate([self.r[0] - dr0 * np.arange(NG, 0, -1),
                                  self.r,
                                  self.r[-1] + self.dr[-1] * np.arange(1, NG + 1)])
        a = self.alpha
        self.A = self.rf**a
        self.V = (self.rf[1:]**(a + 1) - self.rf[:-1]**(a + 1)) / (a + 1)
        self.dA = (self.A[1:] - self.A[:-1])
        self.U = np.zeros((3, n))
        self.t = 0.0
        self.nstep = 0

    # ------------------------------------------------------------------ #
    def set_primitive(self, rho, u, P):
        self.U = self.prim_to_cons(np.asarray(rho, float) * np.ones(self.n),
                                   np.asarray(u, float) * np.ones(self.n),
                                   np.asarray(P, float) * np.ones(self.n))

    def prim_to_cons(self, rho, u, P):
        return np.array([rho, rho * u, P / (self.g - 1) + 0.5 * rho * u**2])

    def cons_to_prim(self, U):
        rho = np.maximum(U[0], self.rho_floor)
        u = U[1] / rho
        P = (self.g - 1) * (U[2] - 0.5 * rho * u**2)
        P = np.maximum(P, self.P_floor)
        return rho, u, P

    def primitive(self):
        return self.cons_to_prim(self.U)

    # ------------------------------------------------------------------ #
    def _ghost(self, W, t):
        """Return primitive arrays including NG ghost cells on each side."""
        rho, u, P = W
        out = []
        for q, sgn in ((rho, 1), (u, -1), (P, 1)):
            out.append(np.concatenate([np.zeros(NG), q, np.zeros(NG)]))
        rho_g, u_g, P_g = out
        # inner
        if callable(self.inner_bc):
            rr = self.rg[:NG]
            rb, ub, Pb = self.inner_bc(rr, t)
            rho_g[:NG], u_g[:NG], P_g[:NG] = rb, ub, Pb
        elif self.inner_bc == "reflect":
            rho_g[:NG] = rho[NG - 1::-1]
            u_g[:NG] = -u[NG - 1::-1]
            P_g[:NG] = P[NG - 1::-1]
        else:
            rho_g[:NG], u_g[:NG], P_g[:NG] = rho[0], u[0], P[0]
        # outer
        if callable(self.outer_bc):
            rr = self.rg[-NG:]
            rb, ub, Pb = self.outer_bc(rr, t)
            rho_g[-NG:], u_g[-NG:], P_g[-NG:] = rb, ub, Pb
        elif self.outer_bc == "reflect":
            rho_g[-NG:] = rho[:-NG - 1:-1]
            u_g[-NG:] = -u[:-NG - 1:-1]
            P_g[-NG:] = P[:-NG - 1:-1]
        else:
            rho_g[-NG:], u_g[-NG:], P_g[-NG:] = rho[-1], u[-1], P[-1]
        return rho_g, u_g, P_g

    def _interface_states(self, Wg):
        """Left/right states at the n+1 physical faces."""
        if self.order == 1:
            L = [q[NG - 1:-NG] for q in Wg]
            R = [q[NG:-NG + 1 if NG > 1 else None] for q in Wg]
            return L, R
        L, R = [], []
        for q in Wg:
            dq = np.diff(q)                       # q_{i+1}-q_i, length n+2NG-1
            s = self.limiter(dq[:-1], dq[1:])     # slope for cells 1..n+2NG-2
            # s[k] is slope (per cell) of cell k+1 in ghosted indexing
            qc = q[1:-1]
            qL_face = qc + 0.5 * s                # right edge of cell
            qR_face = qc - 0.5 * s                # left edge of cell
            # physical faces j = 0..n lie between ghosted cells NG-1+j and NG+j
            L.append(qL_face[NG - 2:NG - 2 + self.n + 1])
            R.append(qR_face[NG - 1:NG - 1 + self.n + 1])
        # positivity safeguard: fall back to first order where needed
        bad = (L[0] <= 0) | (R[0] <= 0) | (L[2] <= 0) | (R[2] <= 0)
        if np.any(bad):
            for k, q in enumerate(Wg):
                L[k] = np.where(bad, q[NG - 1:NG + self.n], L[k])
                R[k] = np.where(bad, q[NG:NG + self.n + 1], R[k])
        return L, R

    def _rhs(self, U, t):
        W = self.cons_to_prim(U)
        Wg = self._ghost(W, t)
        L, R = self._interface_states(Wg)
        F, smax = hllc(L[0], L[1], L[2], R[0], R[1], R[2], self.g)
        AF = F * self.A
        dU = -(AF[:, 1:] - AF[:, :-1]) / self.V
        dU[1] += W[2] * self.dA / self.V
        if self.source is not None:
            dU += self.source(self, U, t)
        return dU

    def timestep(self):
        """CFL time step, including the boundary (ghost) states so that
        fast inflow through a boundary is accounted for."""
        rho, u, P = self._ghost(self.primitive(), self.t)
        c = np.sqrt(self.g * P / rho)
        s = np.abs(u) + c
        dr = np.concatenate([np.full(NG, self.dr[0]), self.dr, np.full(NG, self.dr[-1])])
        return self.cfl * np.min(dr / s)

    def _fix(self, U):
        rho, u, P = self.cons_to_prim(U)
        return self.prim_to_cons(rho, u, P)

    def step(self, dt):
        t = self.t
        if self.order == 1:
            self.U = self._fix(self.U + dt * self._rhs(self.U, t))
        else:
            U1 = self._fix(self.U + dt * self._rhs(self.U, t))
            self.U = self._fix(0.5 * self.U + 0.5 * (U1 + dt * self._rhs(U1, t + dt)))
        self.t += dt
        self.nstep += 1

    def evolve(self, t_stop, callback=None, max_steps=10**9):
        """Advance to exactly t_stop."""
        while self.t < t_stop * (1 - 1e-14) and self.nstep < max_steps:
            dt = min(self.timestep(), t_stop - self.t)
            self.step(dt)
            if callback is not None:
                callback(self)
        return self


def uniform_injection_source(r_inj, mdot, edot):
    """Source term for uniform mass/energy injection inside r < r_inj.

    Uses exact volume overlap of each (spherical) cell with the injection
    sphere.  mdot, edot in code units (total rates).
    """
    def src(solver, U, t):
        rf = solver.rf
        a = solver.alpha
        lo = np.minimum(rf[:-1], r_inj)
        hi = np.minimum(rf[1:], r_inj)
        Vov = (hi**(a + 1) - lo**(a + 1)) / (a + 1)
        Vinj = r_inj**(a + 1) / (a + 1)
        frac = Vov / solver.V / Vinj
        if a == 2:   # rates are for the full 4 pi sphere; V omits the 4 pi
            frac = frac / (4 * np.pi)
        dU = np.zeros_like(U)
        dU[0] = mdot * frac
        dU[2] = edot * frac
        return dU
    return src
