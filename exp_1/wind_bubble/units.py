"""
Unit system for the wind-bubble model.

All internal calculations are performed in CGS units.  This module provides
physical constants (CGS) and conversion factors so that inputs can be given
in convenient astrophysical units and outputs converted back.

Usage
-----
>>> from wind_bubble import units as u
>>> mdot_cgs = 1e-2 * u.Msun / u.yr          # g/s
>>> t_cgs = 1e5 * u.yr                        # s
>>> R_pc = R_cgs / u.pc                       # cm -> pc

i.e. *multiply* by a unit to convert into CGS, *divide* by a unit to convert
out of CGS.
"""

# --- fundamental / physical constants (CGS; CODATA 2018 / IAU 2015) -------
k_B = 1.380649e-16          # Boltzmann constant            [erg / K]
m_p = 1.67262192369e-24     # proton mass                   [g]
m_H = 1.6735575e-24         # hydrogen-atom mass            [g]
G = 6.67430e-8              # gravitational constant        [cm^3 g^-1 s^-2]

# --- base units expressed in CGS -------------------------------------------
cm = 1.0
g = 1.0
s = 1.0
erg = 1.0
K = 1.0
dyn = 1.0

m = 1.0e2 * cm
km = 1.0e5 * cm
au = 1.495978707e13 * cm
pc = 3.0856775814913673e18 * cm
kpc = 1.0e3 * pc

kg = 1.0e3 * g
Msun = 1.988409870698051e33 * g          # IAU 2015 nominal (GM_sun / G)

yr = 3.15576e7 * s                       # Julian year
kyr = 1.0e3 * yr
Myr = 1.0e6 * yr

km_s = km / s

# --- derived conveniences ---------------------------------------------------
cm3 = cm**3
per_cm3 = 1.0 / cm3
erg_s = erg / s
dyn_cm2 = dyn / cm**2
K_cm3 = k_B                              # P/k_B in K cm^-3 -> dyn cm^-2


def mean_mass_per_particle(mu_H=1.4, x_He=None, state="ionized"):
    """Mean mass per free particle, in units of m_p.

    Parameters
    ----------
    mu_H : float
        Mean mass per hydrogen nucleus in units of m_p (default 1.4).
    x_He : float, optional
        Helium abundance by number relative to H.  If ``None`` it is inferred
        from ``mu_H`` assuming the gas is H + He only: mu_H = 1 + 4 x_He.
    state : {"ionized", "atomic", "molecular"}
        Ionisation/molecular state of the gas.  "ionized" means fully ionised
        H and doubly ionised He.

    Returns
    -------
    float
        Mean mass per particle in units of m_p.
    """
    if x_He is None:
        x_He = (mu_H - 1.0) / 4.0
    if state == "ionized":
        n_part = 2.0 + 3.0 * x_He           # H+ + e-  ;  He++ + 2 e-
    elif state == "atomic":
        n_part = 1.0 + x_He
    elif state == "molecular":
        n_part = 0.5 + x_He
    else:
        raise ValueError(f"unknown state {state!r}")
    return mu_H / n_part
