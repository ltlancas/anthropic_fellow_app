"""CGS physical constants and unit conversions."""

m_p = 1.67262192e-24      # proton mass [g]
k_B = 1.380649e-16        # Boltzmann constant [erg/K]
M_sun = 1.98847e33        # solar mass [g]
yr = 3.15576e7            # Julian year [s]
kyr = 1.0e3 * yr          # [s]
pc = 3.0856775814913673e18  # parsec [cm]
AU = 1.495978707e13       # astronomical unit [cm]
km = 1.0e5                # [cm]

MU_H = 1.4                # mean molecular weight per H nucleus (task spec)
MU = 1.4                  # mean mass per particle used for T (assumption, see docs)
GAMMA = 5.0 / 3.0

__all__ = ["m_p", "k_B", "M_sun", "yr", "kyr", "pc", "AU", "km",
           "MU_H", "MU", "GAMMA"]
