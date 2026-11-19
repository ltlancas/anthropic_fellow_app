# Adiabatic wind-blown bubble model (Weaver et al. 1977, Section 2) together with
# the Chevalier & Clegg (1985) free-wind solution it uses for the unshocked wind.
# Combined from fb_models.py and wind_solutions.py
# author: Lachlan Lancaster

import numpy as np
from astropy import units as u
from astropy import constants as ac
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

#########################################################################################
###########################   CLASSICAL WIND SOLUTION MODELS   ##########################
#########################################################################################

class WindModel:
    def __init__(self, **kwargs):
        self._set_parameters_parent(**kwargs)
        self._check_parameter_units_parent()

    def _set_parameters_parent(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

        # following CC85, set the default wind paramters in Msun/yr and 1e43 erg/s
        # this very roughly follows parameters for M82
        if "Mdot" not in self.__dict__:
            self.Mdot = 1.0 * u.Msun / u.yr
        
        if "Edot" not in self.__dict__:
            self.Edot = 1.0e43 * u.erg / u.s

    def _check_parameter_units_parent(self):
        t1 = u.get_physical_type(self.Edot)=="power"
        t2 = u.get_physical_type(self.Mdot*u.s)=="mass"
        if not(t1):
            raise ValueError("Units of Edot are incorrect")
        if not(t2):
            raise ValueError("Units of Mdot are incorrect")

    def _check_radius_units(self, r):
        t1 = u.get_physical_type(r)=="length"
        if not t1:
            raise ValueError("Units of r are incorrect")

class CC85Wind(WindModel):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._set_parmeters()
        self._check_parameter_units()
        self._set_derived_parameters()

    def _set_parmeters(self):       
        if "gamma" not in self.__dict__:
            self.gamma = 5./3
        if "R" not in self.__dict__:
            # 100 pc based more on CC85 values for M82
            self.R = 100.0 * u.pc

    def _check_parameter_units(self):
        t1 = u.get_physical_type(self.R)=="length"
        t2 = u.get_physical_type(self.gamma)=="dimensionless"
        if not(t1):
            raise ValueError("Units of R are incorrect")
        if not(t2):
            raise ValueError("gamma should be dimensionless")

    def _set_derived_parameters(self):
        self.vinf = np.sqrt(2*self.Edot/self.Mdot).to(u.km/u.s)

    def mach(self, r):
        self._check_radius_units(r)
        g = self.gamma
        # the below represent Equations 4 & 5 of CC85
        # these give M^-2
        y = (r/self.R).to(" ").value
        if np.isscalar(y):
            y = np.array([y])
            was_scalar = True
        else:
            was_scalar = False
        mm2 = np.zeros_like(y)
        for i in range(len(y)):
            if y[i] < 1:
                p1 = -1.*(3*g+1)
                p2 = 0.5*(g + 1)
                p3 = 5*g + 1
                f = lambda x: ((3*g+x)/(1+3*g))**p1 * ((g-1+2*x)/(g+1))**p2 - y[i]**p3
                mm2[i] = brentq(f,1,16*(y[i]**-2))
            else:
                p1 = 0.5*(g + 1)
                p2 = 2*(g - 1)
                f = lambda x: x**-1 * ((g-1+2*x)/(1+g))**p1 - y[i]**p2
                mm2[i] = brentq(f,1e-15,1)
        if was_scalar:
            return mm2[0]**-0.5
        else:
            return mm2**-0.5

    def c(self, r):
        # sound speed profile
        mm = self.mach(r)
        return self.vinf / np.sqrt(mm**2 + 2./(self.gamma - 1))

    def u(self, r):
        # velocity profile
        mm = self.mach(r)
        cc = self.c(r)
        return mm * cc
    
    def rho(self, r):
        uu =self.u(r)
        res = self.Mdot/(4*np.pi*r**2*uu)
        return res.to(u.g/u.cm**3)

    def press(self, r):
        cc = self.c(r)
        return cc**2 * self.rho(r) / self.gamma


#########################################################################################
###########################   ADIABATIC WIND BUBBLE MODEL   #############################
#########################################################################################

class Bubble:
    def __init__(self, **kwargs):
        self._set_parmeters_parent(**kwargs)
        self._check_parameter_units_parent()

    def _set_parmeters_parent(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

        if "rho0" not in self.__dict__:
            self.rho0 = 140*ac.m_p/(u.cm**3)

    def _check_parameter_units_parent(self):
        t1 = u.get_physical_type(self.rho0)=="mass density"
        if not(t1):
            raise ValueError("Units of rho0 are incorrect")
    
    def _check_time_units(self, t):
        t1 = u.get_physical_type(t)=="time"
        if not t1:
            raise ValueError("Units of t are incorrect")

    def _check_radius_units(self, r):
        r1 = u.get_physical_type(r)=="length"
        if not r1:
            raise ValueError("Units of r are incorrect")

class AdiabaticWind(Bubble):
    # Weaver solution Section 2 for an adiabatic wind bubble
    # assumes no radiative losses, even in the shell.
    # We don't treat conduction here either

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._set_parmeters(**kwargs)
        self._check_parameter_units()
        self._ad_shell_solve(-2./3)
        self._set_derived_parameters()

        

        # set free-wind solution
        fw_dict = {"Mdot": self.Mdotw, "Edot": self.Lwind,
                   "R": self.rfb, "gamma":self.gamma}
        self.free_wind = CC85Wind(**fw_dict)

    def _set_parmeters(self, **kwargs):
        # scaling paramter for dimensional analysis solution
        # given after Equation 13 of Weaver et al. (1977)
        self.alpha = 0.88

        if "Lwind" not in self.__dict__:
            self.Lwind = 1e38*u.erg/u.s
        if "Mdotw" not in self.__dict__:
            self.Mdotw = 1e-4*u.Msun/u.yr
        if "rfb" not in self.__dict__:
            self.rfb = 1.0*u.pc
        if "gamma" not in self.__dict__:
            self.gamma = 5./3
    
    def _check_parameter_units(self):
        t1 = u.get_physical_type(self.Lwind)=="power"
        t2 = u.get_physical_type(self.rfb)=="length"
        t3 = u.get_physical_type(self.Mdotw*u.s)=="mass"
        if not(t1):
            raise ValueError("Units of L_wind are incorrect")
        if not(t2):
            raise ValueError("Units of r_fb are incorrect")
        if not(t3):
            raise ValueError("Units of Mdot_w are incorrect")

    def _set_derived_parameters(self):
        # fraction of the shell's outer radius at which the shell's inner radius lies
        # approximate 0.86, but determined here from the numerical solution
        xic = self.ad_shell_sol.t[-1]
        # the dimensionless value of the pressure in the shell at the inner edge of the
        # shell radius. Approx 0.59 but determined here from the numerical solution
        Pxic = self.ad_shell_sol.y[2,-1]
        g = self.gamma
        # the dimensionless pre-factor in the scaling solution. Approx 0.88 but 
        # determined here for general gamma and the numerical solution
        self.alpha = (125*(g - 1)/(12*np.pi*xic**3*Pxic*(9*g - 4)))**0.2
        (self.xic, self.Pxic) = (xic, Pxic)
        self.Vwind = np.sqrt(2*self.Lwind/self.Mdotw).to("km/s")

    #################################################################
    #################   TOP-LINE DEFAULT FUNCTIONS   ################
    #################################################################

    def radius(self, t):
        self._check_time_units(t)
        r_we = self.alpha*(self.Lwind*(t**3)/self.rho0)**(1./5)
        return r_we.to("pc")
    
    def velocity(self, t):
        self._check_time_units(t)
        v_we = 0.6*self.radius(t)/t
        return v_we.to("km/s")
    
    def momentum(self, t):
        self._check_time_units(t)
        pr_we = 4*np.pi*self.rho0*self.radius(t)**3*self.velocity(t)/3
        return pr_we.to("solMass*km/s")
    
    def pressure(self, t):
        self._check_time_units(t)
        g = self.gamma
        prefac = 15*(g-1)/(4*np.pi*(9*g-4)*(self.xic*self.alpha)**3)
        press_we = prefac*(self.Lwind**2 * self.rho0**3 / t**4)**(1./5)
        return (press_we/ac.k_B).to("K/cm3")
    
    def density_profile(self, r, t):
        # returns the density of the bubble at radius r and time t
        # r : radius (generally an array)
        # t : time (should be a scalar)
        self._check_time_units(t)
        self._check_radius_units(r)
        g = self.gamma
        r_rs = self.R_rs(t)
        r_b = self.radius(t)
        r_c = self.xic*r_b
        rho_fw = lambda r: self.free_wind.rho(r)
        rho_sw = lambda r: ((g+1)/(g-1))*self.free_wind.rho(r_rs)
        r_sh_ref = r_b*self.ad_shell_sol.t[::-1]
        rho_sh_ref = self.rho0*self.ad_shell_sol.y[1,::-1]
        rho_sh = lambda r: np.interp(r, r_sh_ref, rho_sh_ref)
        rho_bg = lambda r: self.rho0
        rho = np.piecewise(r, [r<r_rs, (r>=r_rs) & (r<r_c), (r>=r_c) & (r<=r_b), r>r_b],\
                           [rho_fw, rho_sw, rho_sh, rho_bg])
        return rho.to("solMass/pc3")
    
    def velocity_profile(self, r, t):
        # returns the velocity of the bubble at radius r and time t
        # r : radius (generally an array)
        # t : time (should be a scalar)
        self._check_time_units(t)
        self._check_radius_units(r)
        g = self.gamma
        r_rs = self.R_rs(t)
        r_b = self.radius(t)
        r_c = self.xic*r_b
        u_fw = lambda r: self.free_wind.u(r)
        u_sw = lambda r: self._v_sw(r, t)
        r_sh_ref = r_b*self.ad_shell_sol.t[::-1]
        u_sh_ref = self.velocity(t)*self.ad_shell_sol.y[0,::-1]
        u_sh = lambda r: np.interp(r, r_sh_ref, u_sh_ref)
        u_bg = lambda r: 0
        u = np.piecewise(r, [r<r_rs, (r>=r_rs) & (r<r_c), (r>=r_c) & (r<=r_b), r>r_b],\
                           [u_fw, u_sw, u_sh, u_bg])
        return u.to("km/s")
    
    def pressure_profile(self, r, t):
        # returns the pressure of the bubble at radius r and time t
        # r : radius (generally an array)
        # t : time (should be a scalar)
        self._check_time_units(t)
        self._check_radius_units(r)
        g = self.gamma
        r_rs = self.R_rs(t)
        r_b = self.radius(t)
        r_c = self.xic*r_b
        p_fw = lambda r: self.free_wind.press(r)
        p_sw = lambda r: self.pressure(t)*ac.k_B
        r_sh_ref = r_b*self.ad_shell_sol.t[::-1]
        p_sh_ref = self.velocity(t)**2*self.rho0*self.ad_shell_sol.y[2,::-1]
        p_sh = lambda r: np.interp(r, r_sh_ref, p_sh_ref)
        p_bg = lambda r: self.rho0*(1*u.km/u.s)**2
        p = np.piecewise(r, [r<r_rs, (r>=r_rs) & (r<r_c), (r>=r_c) & (r<=r_b), r>r_b],\
                           [p_fw, p_sw, p_sh, p_bg])
        return (p/ac.k_B).to("K/cm3")
    
    # TODO: make a separate file for shell structure solutions
    #       so this is similar to the free-wind implementation here 
    #     - Also, replace the string versions of the 
    #       param.to() unit conversion calls
    #################################################################
    ################ FUNCTIONS FOR INTERNAL STRUCTURE ###############
    #################################################################

    def _ad_shell_solve(self, kappa):
        """
        Solves the structure equation for the dimensionless parameters of the shell
        surrounding an adiabatic wind bubble following section 2 of Weaver et al. (1977).
        Args:
            kappa (float): "second order deceleration parameter" of the shell
                           related to the power-law index of the shell radius in time
                           Equal to -2/3 for the R ~ t^(3/5) solution
        Returns:
            None: sets self.ad_shell_sol to the solution object from solve_ivp
        """
        gamma = self.gamma

        def derivs(xi, ys):
            (U, G, P) = ys
            t1 = kappa*G*(U-xi)/P - 2*gamma/xi - 2*kappa/U
            t2 = gamma - (U-xi)**2 *G/P
            Up = U*(t1/t2)
            t1 = Up + 2*U/xi
            t2 = U-xi
            Gp = -G*(t1/t2)
            Pp = P*(gamma*Gp/G - 2*kappa/(U-xi))
            return (Up, Gp, Pp)

        # stop if density goes to 0
        def event_1(t, ys):
            return ys[1]
        event_1.terminal = True

        U0 = 2./(gamma + 1)
        G0 = (gamma + 1)/(gamma - 1)
        P0 = 2/(gamma + 1)
        if not(hasattr(self, "ad_shell_sol")):
            self.ad_shell_sol = solve_ivp(derivs, (1, 0.75), [U0, G0, P0],\
                                          events=[event_1], dense_output=True,\
                                          rtol=1e-12, atol = 1e-12)
        return None
    
    def _v_sw(self, r, t):
        # gives the radial velocity in the shocked wind region
        self._check_radius_units(r)
        self._check_time_units(t)
        g = self.gamma
        r_c = self.xic*self.radius(t)
        gfac = (9*g-4)/(15*g)
        t1 = (gfac*r_c**3/(r**2*t)).to("km/s")
        t2 = ((4/(15*g))*(r/t)).to("km/s")
        return t1 + t2
    
    def R_rs(self, t):
        # gives the radius of the reverse shock as a function of time in the
        # adiabatic wind bubble solution
        self._check_time_units(t)
        Rc = self.xic*self.radius(t)
        R_ballistic = self.Vwind*t
        g = self.gamma
        gfac = ((g+1)/(g-1)) * ((9*g-4)/(15*g)) * ((g+1)**2/(4*g))**(1/(g-1))
        res = np.sqrt(gfac*Rc**3/R_ballistic)
        return res.to("pc")
