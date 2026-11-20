"""Approximate Riemann solvers (vectorised over interfaces)."""
import numpy as np


def _flux(rho, u, P, E):
    return np.array([rho * u, rho * u * u + P, (E + P) * u])


def hllc(rhoL, uL, PL, rhoR, uR, PR, g):
    """HLLC flux (Toro 2009, Sec. 10.4) with Davis wave-speed estimates."""
    cL = np.sqrt(g * PL / rhoL)
    cR = np.sqrt(g * PR / rhoR)
    EL = PL / (g - 1) + 0.5 * rhoL * uL**2
    ER = PR / (g - 1) + 0.5 * rhoR * uR**2
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)
    dL = rhoL * (SL - uL)
    dR = rhoR * (SR - uR)
    Ss = (PR - PL + uL * dL - uR * dR) / (dL - dR)

    FL = _flux(rhoL, uL, PL, EL)
    FR = _flux(rhoR, uR, PR, ER)

    def ustar(rho, u, P, E, S, d):
        fac = d / (S - Ss)
        return np.array([fac,
                         fac * Ss,
                         fac * (E / rho + (Ss - u) * (Ss + P / d))])

    UL = np.array([rhoL, rhoL * uL, EL])
    UR = np.array([rhoR, rhoR * uR, ER])
    FsL = FL + SL * (ustar(rhoL, uL, PL, EL, SL, dL) - UL)
    FsR = FR + SR * (ustar(rhoR, uR, PR, ER, SR, dR) - UR)

    F = np.where(SL >= 0, FL,
                 np.where(Ss >= 0, FsL,
                          np.where(SR > 0, FsR, FR)))
    return F, np.maximum(np.abs(SL), np.abs(SR))
