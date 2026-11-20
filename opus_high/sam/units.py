"""Physical constants (cgs) and a simple code-unit system.

Code units used by the hydro code:  [L] = 1 pc, [t] = 1 kyr, [M] = 1 Msun.
All derived quantities (density, velocity, pressure, energy, luminosity)
follow from these three base units.
"""

# --- physical constants (cgs) ---
M_P = 1.67262192e-24      # proton mass [g]
K_B = 1.380649e-16        # Boltzmann constant [erg/K]
MSUN = 1.98847e33         # solar mass [g]
PC = 3.0856775814913673e18  # parsec [cm]
AU = 1.495978707e13       # astronomical unit [cm]
YR = 3.15576e7            # Julian year [s]
KYR = 1e3 * YR
KMS = 1e5                 # km/s in cm/s

MU_H = 1.4                # mean molecular weight per H nucleus (task spec)
GAMMA = 5.0 / 3.0


class CodeUnits:
    """Conversion between cgs and code units ([L]=pc, [t]=kyr, [M]=Msun)."""

    def __init__(self, length=PC, time=KYR, mass=MSUN):
        self.length = length
        self.time = time
        self.mass = mass
        self.velocity = length / time
        self.density = mass / length**3
        self.pressure = mass / (length * time**2)   # = energy density
        self.energy = mass * length**2 / time**2
        self.luminosity = self.energy / time
        self.mdot = mass / time

    # convenience helpers
    def to_cgs(self, value, kind):
        return value * getattr(self, kind)

    def to_code(self, value, kind):
        return value / getattr(self, kind)


CODE = CodeUnits()


def nH_to_rho(nH):
    """Hydrogen number density [cm^-3] -> mass density [g cm^-3]."""
    return MU_H * M_P * nH


def temperature(P, rho):
    """Temperature [K] from cgs pressure and density, using mu = MU_H."""
    return MU_H * M_P * P / (rho * K_B)


def pressure_from_T(rho, T):
    return rho * K_B * T / (MU_H * M_P)
