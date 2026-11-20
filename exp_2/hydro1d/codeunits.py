"""Code unit system for the hydro solver.

Base units: length = 1 pc, time = 1 kyr, mass density = 1.4 m_p cm^-3
(i.e. density in units of the hydrogen number density for mu_H = 1.4).
Derived: velocity = pc/kyr (~977.8 km/s), pressure = rho_u v_u^2,
energy density = pressure unit, mass rate, luminosity etc.
"""
from windbubble import units as U


class CodeUnits:
    def __init__(self, length=U.PC, time=U.KYR, density=U.MU_H * U.M_P):
        self.length = length
        self.time = time
        self.density = density
        self.velocity = length / time
        self.pressure = density * self.velocity**2
        self.mass = density * length**3
        self.mdot = self.mass / time
        self.luminosity = self.mass * self.velocity**2 / time
        self.energy = self.mass * self.velocity**2

    # cgs -> code
    def to_code(self, value, kind):
        return value / getattr(self, kind)

    # code -> cgs
    def to_cgs(self, value, kind):
        return value * getattr(self, kind)
