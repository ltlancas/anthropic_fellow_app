"""1D finite-volume Euler solver (planar or spherical), first or second
order in space and time, HLLC Riemann solver, optional uniform mass/energy
injection inside a radius R_inj.

State arrays have shape (3, n): primitives W = (rho, v, P),
conserved U = (rho, rho v, E).  Everything is vectorised over cells.
"""
import numpy as np

NG = 2  # ghost cells per side


class Hydro:
    def __init__(self, N, r_max, r_min=0.0, gamma=5.0 / 3.0, order=2,
                 geometry="spherical", cfl=0.4, inner_bc="reflect",
                 outer_bc="outflow", rho_floor=1e-12, P_floor=1e-14,
                 Mdot=0.0, L=0.0, R_inj=0.0):
        assert order in (1, 2)
        assert geometry in ("spherical", "cartesian")
        self.N, self.gamma, self.order, self.cfl = N, gamma, order, cfl
        self.geometry = geometry
        self.inner_bc, self.outer_bc = inner_bc, outer_bc
        self.rho_floor, self.P_floor = rho_floor, P_floor

        self.faces = np.linspace(r_min, r_max, N + 1)
        self.r = 0.5 * (self.faces[1:] + self.faces[:-1])
        self.dr = np.diff(self.faces)
        if geometry == "spherical":
            self.area = 4 * np.pi * self.faces**2
            self.vol = 4 * np.pi / 3 * np.diff(self.faces**3)
        else:
            self.area = np.ones(N + 1)
            self.vol = self.dr.copy()
        self.dA = (self.area[1:] - self.area[:-1]) / self.vol  # geometric source

        # injection: per-volume rates times fraction of cell inside R_inj
        self.src_mass = np.zeros(N)
        self.src_energy = np.zeros(N)
        if R_inj > 0 and (Mdot > 0 or L > 0):
            if geometry == "spherical":
                Vinj = 4 * np.pi / 3 * R_inj**3
                lo = np.minimum(self.faces[:-1], R_inj)
                hi = np.minimum(self.faces[1:], R_inj)
                frac = (4 * np.pi / 3 * (hi**3 - lo**3)) / self.vol
            else:
                Vinj = R_inj - r_min
                lo = np.minimum(self.faces[:-1], R_inj)
                hi = np.minimum(self.faces[1:], R_inj)
                frac = (hi - lo) / self.vol
            self.src_mass = Mdot / Vinj * frac
            self.src_energy = L / Vinj * frac

        self.t = 0.0
        self.U = np.zeros((3, N))
        self.nsteps = 0

    # ---------------------------------------------------------------- utils
    def prim_to_cons(self, W):
        rho, v, P = W
        return np.array([rho, rho * v, P / (self.gamma - 1) + 0.5 * rho * v**2])

    def cons_to_prim(self, U):
        rho = np.maximum(U[0], self.rho_floor)
        v = U[1] / rho
        P = (self.gamma - 1) * (U[2] - 0.5 * rho * v**2)
        P = np.maximum(P, self.P_floor)
        return np.array([rho, v, P])

    def set_primitive(self, rho, v, P):
        W = np.array([np.broadcast_to(rho, self.r.shape),
                      np.broadcast_to(v, self.r.shape),
                      np.broadcast_to(P, self.r.shape)], dtype=float)
        self.U = self.prim_to_cons(W)

    @property
    def W(self):
        return self.cons_to_prim(self.U)

    def _pad(self, W):
        Wg = np.empty((3, self.N + 2 * NG))
        Wg[:, NG:-NG] = W
        inner = W[:, NG - 1::-1].copy()          # mirror first NG cells
        if self.inner_bc == "reflect":
            inner[1] *= -1
        elif self.inner_bc == "outflow":
            inner[:] = W[:, :1]
        Wg[:, :NG] = inner
        outer = W[:, :-NG - 1:-1].copy()
        if self.outer_bc == "reflect":
            outer[1] *= -1
        else:
            outer[:] = W[:, -1:]
        Wg[:, -NG:] = outer
        return Wg

    # ------------------------------------------------------- reconstruction
    @staticmethod
    def _minmod(a, b):
        return np.where(a * b > 0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)

    def _face_states(self, Wg):
        """Left/right states at the N+1 faces bounding the active cells."""
        if self.order == 1:
            WL = Wg[:, NG - 1:-NG]
            WR = Wg[:, NG:-NG + 1] if NG > 1 else Wg[:, NG:]
            return WL, WR
        dW = np.diff(Wg, axis=1)                  # n-1 differences
        slope = self._minmod(dW[:, :-1], dW[:, 1:])  # cells 1..n-2 (half-slope*2)
        # slope index j corresponds to padded cell j+1
        Wc = Wg[:, 1:-1]
        Wplus = Wc + 0.5 * slope                  # right edge of cell
        Wminus = Wc - 0.5 * slope                 # left edge of cell
        # faces between padded cells NG-1..NG+N-1 and NG..NG+N
        WL = Wplus[:, NG - 2:NG - 2 + self.N + 1]
        WR = Wminus[:, NG - 1:NG - 1 + self.N + 1]
        WL = np.array([np.maximum(WL[0], self.rho_floor), WL[1], np.maximum(WL[2], self.P_floor)])
        WR = np.array([np.maximum(WR[0], self.rho_floor), WR[1], np.maximum(WR[2], self.P_floor)])
        return WL, WR

    # -------------------------------------------------------- Riemann solver
    def _hllc(self, WL, WR):
        g = self.gamma
        rL, uL, pL = WL
        rR, uR, pR = WR
        cL = np.sqrt(g * pL / rL)
        cR = np.sqrt(g * pR / rR)
        EL = pL / (g - 1) + 0.5 * rL * uL**2
        ER = pR / (g - 1) + 0.5 * rR * uR**2
        SL = np.minimum(uL - cL, uR - cR)
        SR = np.maximum(uL + cL, uR + cR)
        Ss = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (
            rL * (SL - uL) - rR * (SR - uR))

        FL = np.array([rL * uL, rL * uL**2 + pL, (EL + pL) * uL])
        FR = np.array([rR * uR, rR * uR**2 + pR, (ER + pR) * uR])
        UL = np.array([rL, rL * uL, EL])
        UR = np.array([rR, rR * uR, ER])

        def ustar(r, u, p, E, S):
            fac = r * (S - u) / (S - Ss)
            return fac * np.array([np.ones_like(r), Ss,
                                   E / r + (Ss - u) * (Ss + p / (r * (S - u)))])

        FsL = FL + SL * (ustar(rL, uL, pL, EL, SL) - UL)
        FsR = FR + SR * (ustar(rR, uR, pR, ER, SR) - UR)
        F = np.where(SL >= 0, FL, np.where(Ss >= 0, FsL, np.where(SR > 0, FsR, FR)))
        return F

    # ---------------------------------------------------------------- RHS
    def rhs(self, U):
        W = self.cons_to_prim(U)
        Wg = self._pad(W)
        WL, WR = self._face_states(Wg)
        F = self._hllc(WL, WR)
        AF = F * self.area
        dU = -(AF[:, 1:] - AF[:, :-1]) / self.vol
        dU[1] += W[2] * self.dA                   # geometric pressure source
        dU[0] += self.src_mass
        dU[2] += self.src_energy
        return dU

    def timestep(self):
        W = self.W
        c = np.sqrt(self.gamma * W[2] / W[0])
        dt = self.cfl * np.min(self.dr / (np.abs(W[1]) + c))
        m = self.src_energy > 0
        if m.any():  # source limiter: energy injected per step <= 0.3 * cell energy
            dt = min(dt, 0.3 * np.min(self.U[2][m] / self.src_energy[m]))
        return dt

    def step(self, dt):
        if self.order == 1:
            self.U = self.U + dt * self.rhs(self.U)
        else:  # SSP-RK2 (Heun)
            U1 = self.U + dt * self.rhs(self.U)
            U1 = self.prim_to_cons(self.cons_to_prim(U1))
            self.U = 0.5 * (self.U + U1 + dt * self.rhs(U1))
        self.U = self.prim_to_cons(self.cons_to_prim(self.U))  # apply floors
        self.t += dt
        self.nsteps += 1

    def evolve(self, t_end, output_times=(), callback=None, max_steps=10**8):
        """Advance to t_end; call callback(self, k) exactly at each output time."""
        outs = sorted(output_times)
        k = 0
        while self.t < t_end * (1 - 1e-14) and self.nsteps < max_steps:
            dt = self.timestep()
            target = outs[k] if k < len(outs) else t_end
            if self.t + dt >= target:
                dt = target - self.t
            self.step(dt)
            if k < len(outs) and abs(self.t - outs[k]) <= 1e-12 * max(1.0, outs[k]):
                self.t = outs[k]
                if callback is not None:
                    callback(self, k)
                k += 1
        return self

    # --------------------------------------------------------- diagnostics
    def totals(self):
        return (self.U * self.vol).sum(axis=1)
