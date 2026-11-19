"""
Classic self-similar energy-driven wind-bubble solution
(Castor, McCray & Weaver 1975, ApJL 200, L107; Weaver et al. 1977, ApJ 218, 377).

Valid when (i) the reverse-shock radius R1 << R2, (ii) the swept-up mass
dominates the injected wind mass, and (iii) radiative losses from the hot
bubble are negligible.  All quantities in CGS.
"""
import numpy as np

GAMMA = 5.0 / 3.0

_A_R = (125.0 / (154.0 * np.pi)) ** 0.2           # ~0.763
_A_P = 7.0 / (3850.0 * np.pi) ** 0.4              # ~0.163


def outer_radius(L, rho0, t):
    """Shell radius R2(t) = (125/154pi)^(1/5) (L t^3 / rho0)^(1/5)  [cm]."""
    return _A_R * (L * t**3 / rho0) ** 0.2


def shell_velocity(L, rho0, t):
    """dR2/dt = 3 R2 / (5 t)  [cm/s]."""
    return 0.6 * outer_radius(L, rho0, t) / t


def bubble_pressure(L, rho0, t):
    """Uniform hot-bubble pressure P = 7/(3850 pi)^(2/5) L^(2/5) rho0^(3/5) t^(-4/5)."""
    return _A_P * L**0.4 * rho0**0.6 * t ** (-0.8)


def bubble_thermal_energy(L, t):
    """Thermal energy of the hot bubble, E_th = (5/11) L t."""
    return 5.0 / 11.0 * L * t


def reverse_shock_radius(L, Mdot, rho0, t, gamma=GAMMA):
    """R1 from strong-shock pressure balance P = 2/(gamma+1) rho_w v_w^2 (R1 << R2)."""
    v_w = np.sqrt(2.0 * L / Mdot)
    P = bubble_pressure(L, rho0, t)
    return np.sqrt(2.0 * Mdot * v_w / ((gamma + 1.0) * 4.0 * np.pi * P))
