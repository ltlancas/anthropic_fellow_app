"""Self-similar structure of the adiabatic shocked-ambient shell of a
wind-blown bubble expanding as R2 ∝ t^alpha into a uniform medium.

Variables:  lam = r/R2,  v = R2dot V(lam),  rho = rho0 G(lam),
            P = rho0 R2dot^2 Pi(lam).
See assumptions.md for the derivation of the ODEs.
"""
import numpy as np
from scipy.integrate import solve_ivp

from .units import GAMMA


def _rhs(lam, y, alpha, g):
    G, V, Pi = y
    c1 = (1 - alpha) / alpha
    D = V - lam
    cs2 = g * Pi / G
    dV = (c1 * V * D - (Pi / G) * (2 * c1 - 2 * g * V / lam)) / (D**2 - cs2)
    dG = -G * (dV + 2 * V / lam) / D
    dPi = (Pi / D) * (2 * c1 - g * dV - 2 * g * V / lam)
    return [dG, dV, dPi]


def solve_shell(alpha=0.6, g=GAMMA, n=4000):
    """Integrate from the outer shock (lam=1) inward to the contact
    discontinuity (where V = lam).  Returns dict with profiles and
    integrals."""
    y0 = [(g + 1) / (g - 1), 2 / (g + 1), 2 / (g + 1)]

    def contact(lam, y, *a):
        return (y[1] - lam) + 1e-9
    contact.terminal = True

    sol = solve_ivp(_rhs, (1.0, 0.01), y0, args=(alpha, g), events=contact,
                    rtol=1e-11, atol=1e-13, dense_output=True, method="LSODA")
    lam_c = sol.t[-1]
    lam = np.linspace(1.0, lam_c, n)
    G, V, Pi = sol.sol(lam)
    # energy integral  int lam^2 (G V^2/2 + Pi/(g-1)) dlam  over shell
    integrand = lam**2 * (0.5 * G * V**2 + Pi / (g - 1))
    I_shell = np.trapezoid(integrand[::-1], lam[::-1])
    M_shell = np.trapezoid((lam**2 * G)[::-1], lam[::-1])
    return dict(lam=lam[::-1], G=G[::-1], V=V[::-1], Pi=Pi[::-1],
                lam_c=lam_c, Pi_c=Pi[-1], I_shell=I_shell, M_shell=M_shell,
                alpha=alpha, gamma=g)


_CACHE = {}


def shell_solution(alpha=0.6, g=GAMMA):
    key = (alpha, g)
    if key not in _CACHE:
        _CACHE[key] = solve_shell(alpha, g)
    return _CACHE[key]


if __name__ == "__main__":
    s = shell_solution()
    a, g = s["alpha"], s["gamma"]
    Itot = s["I_shell"] + s["lam_c"]**3 * s["Pi_c"] / (3 * (g - 1))
    xi2 = (1 / (4 * np.pi * a**2 * Itot)) ** 0.2
    print("lam_c=%.5f Pi_c=%.5f I_shell=%.5f M_shell*3=%.5f xi2=%.5f"
          % (s["lam_c"], s["Pi_c"], s["I_shell"], 3 * s["M_shell"], xi2))
