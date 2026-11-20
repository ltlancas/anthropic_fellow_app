"""Physical constants and unit conversions (cgs internally)."""

m_p = 1.67262192e-24      # g
k_B = 1.380649e-16        # erg / K
M_sun = 1.98847e33        # g
yr = 3.15576e7            # s (Julian year)
kyr = 1.0e3 * yr
pc = 3.0856775814913673e18  # cm
AU = 1.495978707e13       # cm
km = 1.0e5                # cm

MU_H = 1.4                # mean mass per H nucleus, in m_p (also used as mean molecular weight)
GAMMA = 5.0 / 3.0


def msun_per_yr_to_cgs(mdot):
    return mdot * M_sun / yr


def nH_to_rho(nH, mu=MU_H):
    return nH * mu * m_p


def rho_to_nH(rho, mu=MU_H):
    return rho / (mu * m_p)


def temperature(rho, P, mu=MU_H):
    """T = mu m_p P / (rho k_B)."""
    return mu * m_p * P / (rho * k_B)


def pressure_from_T(rho, T, mu=MU_H):
    return rho * k_B * T / (mu * m_p)
