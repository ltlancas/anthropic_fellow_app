"""Physical constants and unit conversions (CGS internally).

All quantities inside the library are in CGS. Helper functions convert the
astrophysical input units used in the problem specification
(Msun/yr, erg/s, cm^-3, kyr, AU, pc) to and from CGS.
"""

# --- fundamental constants (CODATA 2018 / IAU 2015) -------------------------
M_P = 1.67262192369e-24      # proton mass [g]
K_B = 1.380649e-16           # Boltzmann constant [erg/K]

# --- astronomical units -----------------------------------------------------
MSUN = 1.98847e33            # solar mass [g]
YR = 3.15576e7               # Julian year [s]
KYR = 1.0e3 * YR             # [s]
AU = 1.495978707e13          # astronomical unit [cm]
PC = 3.0856775814913673e18   # parsec [cm]
KM = 1.0e5                   # km [cm]

# --- gas composition ----------------------------------------------------------
# Mean mass per hydrogen nucleus, mu_H = 1.4 (given).  With He/H = 0.1 by number
# this is 1 + 4*0.1 = 1.4.  For the temperature we need the mean mass per free
# particle; we assume fully ionized H and He: per H nucleus there are
# 1 (H+) + 0.1 (He++) + 1.2 (e-) = 2.3 particles, so mu = 1.4/2.3.
MU_H = 1.4
X_HE = 0.1
PARTICLES_PER_H = 1.0 + X_HE + (1.0 + 2.0 * X_HE)   # = 2.3
MU = MU_H / PARTICLES_PER_H                          # ~0.6087

GAMMA = 5.0 / 3.0


def mdot_cgs(mdot_msun_per_yr):
    """Msun/yr -> g/s."""
    return mdot_msun_per_yr * MSUN / YR


def rho_from_nH(n_H, mu_H=MU_H):
    """Hydrogen number density [cm^-3] -> mass density [g cm^-3]."""
    return n_H * mu_H * M_P


def temperature(P, rho, mu=MU):
    """Ideal-gas temperature [K] from pressure [dyn cm^-2], density [g cm^-3]."""
    return mu * M_P * P / (rho * K_B)


def pressure_from_T(T, rho, mu=MU):
    """Ideal-gas pressure [dyn cm^-2] from temperature and density."""
    return rho * K_B * T / (mu * M_P)
