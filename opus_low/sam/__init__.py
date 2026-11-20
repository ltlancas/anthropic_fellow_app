"""Semi-analytic model of an adiabatic wind-blown bubble (Weaver et al. 1977, adiabatic phase)."""
import numpy as np

# cgs constants
MP = 1.67262192e-24
KB = 1.380649e-16
MSUN = 1.98847e33
YR = 3.15576e7
KYR = 1e3 * YR
PC = 3.0856775814913673e18
AU = 1.495978707e13
GAMMA = 5.0 / 3.0
MU = 1.4  # mass per H nucleus [m_p]; also used as mean mass per particle for T


class Bubble:
    def __init__(self, mdot_msun_yr, L, nH, r_inj=100 * AU):
        self.Mdot = mdot_msun_yr * MSUN / YR
        self.L = L
        self.rho0 = nH * MU * MP
        self.vw = np.sqrt(2 * L / self.Mdot)
        self.r_inj = r_inj
        self.P0 = 0.0  # cold ambient medium (pressure neglected)

    # Weaver et al. (1977) adiabatic self-similar results
    def R2(self, t):
        return (125.0 / (154.0 * np.pi)) ** 0.2 * (self.L * t**3 / self.rho0) ** 0.2

    def V2(self, t):
        return 0.6 * self.R2(t) / t

    def P_hot(self, t):
        return 7.0 / (3850.0 * np.pi) ** 0.4 * self.L**0.4 * self.rho0**0.6 * t ** (-0.8)

    def rho_wind(self, r):
        return self.Mdot / (4 * np.pi * np.maximum(r, self.r_inj) ** 2 * self.vw)

    def R1(self, t):
        # reverse shock: post-shock pressure (3/4) rho_w vw^2 equals P_hot
        return np.sqrt(0.75 * self.Mdot * self.vw / (4 * np.pi * self.P_hot(t)))

    def Rc(self, t):
        # contact discontinuity: swept mass in shell of density 4 rho0
        return self.R2(t) * 0.75 ** (1.0 / 3.0)

    def profile(self, r, t):
        """Return rho [g/cc], v [cm/s], P [dyn/cm^2], T [K] at radii r [cm], time t [s]."""
        r = np.asarray(r, float)
        R1, Rc, R2, P = self.R1(t), self.Rc(t), self.R2(t), self.P_hot(t)
        rho = np.full_like(r, self.rho0); v = np.zeros_like(r); p = np.full_like(r, self.P0)
        # free wind (inside injection region: uniform-source solution, v ~ linear)
        fw = r < R1
        rho[fw] = self.rho_wind(r[fw])
        v[fw] = self.vw * np.minimum(r[fw] / self.r_inj, 1.0)
        
        p[fw] = 1e-6 * P  # cold supersonic wind (negligible pressure)
        # shocked wind: isobaric; density conserves postshock entropy (adiabatic, uniform P)
        sw = (r >= R1) & (r < Rc)
        rho[sw] = 4 * self.rho_wind(R1)
        p[sw] = P
        # velocity: linear between postshock vw/4 and contact speed V2 * Rc/R2
        vc = self.V2(t) * Rc / R2
        v[sw] = 0.25 * self.vw + (vc - 0.25 * self.vw) * (r[sw] - R1) / (Rc - R1)
        # swept-up shell: strong-shock jump conditions
        sh = (r >= Rc) & (r < R2)
        rho[sh] = 4 * self.rho0
        p[sh] = 0.75 * self.rho0 * self.V2(t) ** 2
        v[sh] = 0.75 * self.V2(t)
        T = p * MU * MP / (rho * KB)
        return rho, v, p, T
