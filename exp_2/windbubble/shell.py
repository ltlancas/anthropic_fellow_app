"""Quasi-self-similar structure of the adiabatic shocked-ambient shell.

The shell between the contact discontinuity R_c and the forward shock R_2 is
described with the similarity variable xi = r/R_2 and

    rho = rho_0 G(xi),  u = V_2 U(xi),  P = rho_0 V_2^2 Pi(xi),

where V_2 = dR_2/dt.  The only parameter is the deceleration parameter
q = R_2 (d^2R_2/dt^2) / V_2^2 (q = (alpha-1)/alpha for R_2 ~ t^alpha; Weaver:
alpha = 3/5, q = -2/3).  The Euler equations become (W = U - xi):

    W G' + G U'           = -2 U G / xi          (mass)
    W U' + Pi'/G          = -q U                 (momentum)
    W (Pi'/Pi - g G'/G)   = -2 q                 (entropy)

integrated inward from the strong-shock jump at xi = 1 to the contact, where
U = xi (W = 0).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from .units import GAMMA


def _rhs_tau(tau, y, q, g):
    """Regularised system in a pseudo-time tau with d xi/d tau = -W (W^2 - c^2).

    y = (xi, U, lnG, lnPi).  The contact (W -> 0) is approached asymptotically
    as tau -> infinity, which removes the coordinate singularity there.
    """
    xi, U, lnG, lnPi = y
    G, Pi = np.exp(lnG), np.exp(lnPi)
    W = U - xi
    c2 = g * Pi / G
    D = W**2 - c2
    NU = -q * U * W + 2 * q * Pi / G + 2 * c2 * U / xi
    dxi = -W * D
    dU = -W * NU
    dlnG = NU + 2 * U * D / xi
    dlnPi = g * dlnG + 2 * q * D
    return [dxi, dU, dlnG, dlnPi]


class ShellSolution:
    """Self-similar shell profile for a given deceleration parameter q."""

    def __init__(self, q, g=GAMMA, wtol=1e-9):
        self.q = q
        self.g = g
        y0 = [1.0, 2 / (g + 1), np.log((g + 1) / (g - 1)), np.log(2 / (g + 1))]

        def ev(tau, y, *a):
            return (y[1] - y[0]) / y[0] + wtol
        ev.terminal = True

        def ev_sonic(tau, y, *a):
            G, Pi = np.exp(y[2]), np.exp(y[3])
            return (y[1] - y[0])**2 - g * Pi / G
        ev_sonic.terminal = True

        sol = solve_ivp(_rhs_tau, (0.0, 1e6), y0, args=(q, g),
                        events=[ev, ev_sonic], method="DOP853",
                        rtol=1e-11, atol=1e-13, max_step=0.02)
        if sol.status != 1 or len(sol.t_events[0]) == 0:
            raise RuntimeError("shell ODE did not reach contact for q=%g" % q)
        xi, U, lnG, lnPi = sol.y
        # sort ascending in xi
        o = np.argsort(xi)
        xi, U, G, Pi = xi[o], U[o], np.exp(lnG[o]), np.exp(lnPi[o])
        keep = np.concatenate([[True], np.diff(xi) > 1e-14])
        xi, U, G, Pi = xi[keep], U[keep], G[keep], Pi[keep]
        self.xi_c = xi[0]
        self.G_c, self.U_c, self.Pi_c = G[0], U[0], Pi[0]
        self._xi, self._G, self._U, self._Pi = xi, G, U, Pi
        self.mass = np.trapz(4 * np.pi * xi**2 * G, xi)       # should be 4pi/3
        self.kin = np.trapz(4 * np.pi * xi**2 * 0.5 * G * U**2, xi)
        self.th = np.trapz(4 * np.pi * xi**2 * Pi / (g - 1), xi)

    def profile(self, xi):
        """Return (G, U, Pi) at xi in [xi_c, 1] (linear interpolation)."""
        xi = np.clip(xi, self.xi_c, 1.0)
        return (np.interp(xi, self._xi, self._G), np.interp(xi, self._xi, self._U),
                np.interp(xi, self._xi, self._Pi))


class ShellTable:
    """Tabulation of the shell solution over q for fast interpolation."""

    def __init__(self, qmin=-1.2, qmax=0.3, n=151, g=GAMMA):
        self.q = np.linspace(qmin, qmax, n)
        self.sols = [ShellSolution(q, g) for q in self.q]
        self.xi_c = np.array([s.xi_c for s in self.sols])
        self.Pi_c = np.array([s.Pi_c for s in self.sols])
        self._xi_c = interp1d(self.q, self.xi_c, kind="cubic")
        self._Pi_c = interp1d(self.q, self.Pi_c, kind="cubic")

    def xic(self, q):
        return float(self._xi_c(q))

    def pic(self, q):
        return float(self._Pi_c(q))

    def profile(self, q, xi):
        """Profile at arbitrary q: linear interpolation between tabulated q."""
        i = int(np.clip(np.searchsorted(self.q, q) - 1, 0, len(self.q) - 2))
        w = (q - self.q[i]) / (self.q[i + 1] - self.q[i])
        s0, s1 = self.sols[i], self.sols[i + 1]
        # map xi onto each solution's own [xi_c, 1] range so the contact lines up
        xic = (1 - w) * s0.xi_c + w * s1.xi_c
        s = (np.asarray(xi) - xic) / (1 - xic)
        p0 = s0.profile(s0.xi_c + s * (1 - s0.xi_c))
        p1 = s1.profile(s1.xi_c + s * (1 - s1.xi_c))
        return tuple((1 - w) * a + w * b for a, b in zip(p0, p1))
