"""1D spherical finite-volume hydro (HLLE), 1st order (PCM + Euler) or 2nd order (MUSCL minmod + SSP-RK2)."""
import numpy as np

class Hydro:
    def __init__(self, N, rmax, rho0, P0, gamma=5/3, order=2, cfl=0.4,
                 Mdot=0.0, L=0.0, r_inj=0.0):
        self.N, self.g, self.order, self.cfl = N, gamma, order, cfl
        self.rf = np.linspace(0, rmax, N + 1)
        self.r = 0.5 * (self.rf[1:] + self.rf[:-1])
        self.A = 4 * np.pi * self.rf**2
        self.V = 4 / 3 * np.pi * (self.rf[1:]**3 - self.rf[:-1]**3)
        # injection: volume-overlap weighting of a uniform sphere of radius r_inj
        ov = 4 / 3 * np.pi * (np.clip(self.rf[1:], 0, r_inj)**3 - np.clip(self.rf[:-1], 0, r_inj)**3)
        w = ov / ov.sum() if ov.sum() > 0 else ov
        self.qm = Mdot * w / self.V      # mass source per volume per time
        self.qe = L * w / self.V         # energy source per volume per time
        self.U = np.zeros((3, N)); self.U[0] = rho0; self.U[2] = P0 / (gamma - 1)
        self.t = 0.0; self.floor = 1e-12 * P0

    def prim(self, U):
        rho = U[0]; v = U[1] / rho
        P = np.maximum((self.g - 1) * (U[2] - 0.5 * rho * v**2), self.floor)
        return rho, v, P

    def _flux(self, rl, vl, pl, rr, vr, pr):
        g = self.g
        cl, cr = np.sqrt(g * pl / rl), np.sqrt(g * pr / rr)
        sl = np.minimum(vl - cl, vr - cr); sr = np.maximum(vl + cl, vr + cr)
        def F(r, v, p):
            E = p / (g - 1) + 0.5 * r * v**2
            return np.array([r * v, r * v**2 + p, (E + p) * v]), np.array([r, r * v, E])
        Fl, Ul = F(rl, vl, pl); Fr, Ur = F(rr, vr, pr)
        Fh = (sr * Fl - sl * Fr + sl * sr * (Ur - Ul)) / np.where(sr - sl == 0, 1, sr - sl)
        return np.where(sl >= 0, Fl, np.where(sr <= 0, Fr, Fh))

    def rhs(self, U):
        rho, v, P = self.prim(U)
        W = np.array([rho, v, P])
        # ghost cells: reflecting inner, outflow outer (2 each)
        Wg = np.concatenate([W[:, 1::-1] * np.array([[1], [-1], [1]]), W, W[:, -1:], W[:, -1:]], axis=1)
        if self.order == 2:
            dl = Wg[:, 1:-1] - Wg[:, :-2]; dr = Wg[:, 2:] - Wg[:, 1:-1]
            s = np.where(dl * dr > 0, np.sign(dl) * np.minimum(abs(dl), abs(dr)), 0)
        else:
            s = np.zeros_like(Wg[:, 1:-1])
        Wc = Wg[:, 1:-1]  # cells -1..N
        WL = (Wc + 0.5 * s)[:, :-1]; WR = (Wc - 0.5 * s)[:, 1:]  # faces 0..N
        Fx = self._flux(*WL, *WR)
        AF = Fx * self.A
        dU = -(AF[:, 1:] - AF[:, :-1]) / self.V
        dU[1] += P * (self.A[1:] - self.A[:-1]) / self.V
        dU[0] += self.qm; dU[2] += self.qe
        return dU

    def dt(self):
        rho, v, P = self.prim(self.U)
        return self.cfl * np.min((self.rf[1:] - self.rf[:-1]) / (abs(v) + np.sqrt(self.g * P / rho)))

    def _fix(self, U):
        U[0] = np.maximum(U[0], 1e-40); return U

    def step(self, dt):
        if self.order == 1:
            self.U = self._fix(self.U + dt * self.rhs(self.U))
        else:
            U1 = self._fix(self.U + dt * self.rhs(self.U))
            self.U = self._fix(0.5 * (self.U + U1 + dt * self.rhs(U1)))
        self.t += dt

    def evolve_to(self, tend):
        while self.t < tend:
            self.step(min(self.dt(), tend - self.t))
