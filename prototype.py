"""Test script for sanity checking.

Units
-----
length : kpc
velocity : km/s
mass : Msun

=> 
time : kpc/km s ~ 0.978 Gyr
specific energy : km2 / s2
specific angular momentum : kpc km/s
"""

import numpy as np
from scipy import integrate, optimize
from astropy import units as u
from astropy import constants as c

G = 4.302e-6 # kpc Msun-1 km2 s-2
Gyr_per_t = 0.9777922216731284 # Gyr

################################################################################
# Rotation matrices
################################################################################

@u.quantity_input(theta=u.radian)
def Rx(theta):
    s = np.sin(theta)
    c = np.cos(theta)
    return np.array([[1, 0, 0],
                     [0, c, -s],
                     [0, s, c]])

@u.quantity_input(theta=u.radian)
def Ry(theta):
    s = np.sin(theta)
    c = np.cos(theta)
    return np.array([[c, 0, s],
                     [0, 1, 0],
                     [-s, 0, c]])

@u.quantity_input(theta=u.radian)
def Rz(theta):
    s = np.sin(theta)
    c = np.cos(theta)
    return np.array([[c, -s, 0],
                     [s, c, 0],
                     [0, 0, 1]])

################################################################################
# Power law density model
################################################################################

class Potential:
    """Potential, density, enclosed mass, acceleration"""

    def __init__(self, name, params, f_Phi, f_rho, f_M):
        """
        name : str, name
        params : iterable of floats, parameters of model
        f_Phi : r, *params -> Phi, potential function
        f_rho : r, *params -> rho, density function
        f_M : r, *params -> M, enclosed mass function
        """
        self.name = name
        self.params = params
        self.f_Phi = f_Phi
        self.f_rho = f_rho
        self.f_M = f_M

    def __repr__(self):
        return "<{}: {}>".format(self.__class__.__name__, self.name)

    @u.quantity_input(r=u.kpc)
    def __call__(self, r):
        return self.Phi(r)

    @u.quantity_input(r=u.kpc)
    def Phi(self, r):
        return self.f_Phi(r, *self.params)

    @u.quantity_input(r=u.kpc)
    def M(self, r):
        return self.f_M(r, *self.params)

    @u.quantity_input(r=u.kpc)
    def rho(self, r):
        return self.f_rho(r, *self.params)

    def D(self, coord, t):
        """ For Gyr timestep units, coords in kpc, kpc/Gyr
        coord -> Dcoord
        coord = [x, y, vx, vy]
        Dcoord = [vx, vy, ax, ay]
        t is ignored for a time-independent potential
        """
        x, y, vx, vy = coord
        r = np.sqrt(x ** 2 + y ** 2) * u.kpc
        F = (-c.G * self.M(r) / r ** 3).to(u.Gyr ** -2).value
        return np.array([vx, vy, F * x, F * y])

@u.quantity_input(r=u.kpc, r0=u.kpc, rho0=(u.Msun / u.kpc ** 3))    
def rho_pl(r, r0, rho0, alpha):
    rho = rho0 * (r / r0) ** alpha
    return rho.to(u.Msun / u.kpc ** 3)

@u.quantity_input(r=u.kpc, r0=u.kpc, rho0=(u.Msun / u.kpc ** 3))    
def M_pl(r, r0, rho0, alpha):
    M = 4 * np.pi * rho0 * r0 ** 3 / (3 + alpha) * (r / r0) ** (3 + alpha)
    return M.to(u.Msun)

@u.quantity_input(r=u.kpc, r0=u.kpc, rho0=(u.Msun / u.kpc ** 3))    
def Phi_pl(r, r0, rho0, alpha):
    phi = -c.G * M_pl(r, r0, rho0, alpha) / r
    return phi.to(u.km ** 2 / u.s ** 2)

class PowerLawPotential(Potential):
    def __init__(self, r0, rho0, alpha):
        self.r0 = r0
        self.rho0 = rho0
        self.alpha = alpha
        super().__init__("Power law", (r0, rho0, alpha), Phi_pl, rho_pl, M_pl)
    def __repr__(self):
        s = "<{}: r0 = {:.2e}, ρ0 = {:.2e}, α = {:.2f}>"
        return s.format(self.__class__.__name__, self.r0, self.rho0, self.alpha)
        

# potential with 1e12 Msun within 10 kpc
M = 1e12 * u.Msun
r0 = 10 * u.kpc
alpha = -2.5
density_unit = u.Msun / u.kpc ** 3
rho0 = optimize.fsolve(lambda rho: np.abs(M - M_pl(r0, r0, rho * density_unit, alpha)), 1e10)[0] * density_unit
potential_pl = PowerLawPotential(r0, rho0, alpha)

################################################################################
# Orbits
################################################################################

class Orbit:

    @u.quantity_input(e=(u.km / u.s) ** 2, j=(u.kpc * u.km / u.s), i=u.radian,
                      Omega=u.radian)
    def __init__(self, e, j, i, Omega, potential):
        """Orbit in a potential.

        Parameters
        ----------
        e : float, specific energy (-infty < e < 0)
        j : float, specific angular momentum
        i : float, inclination in radians
        Omega : float, longitude of ascending node in radians
        potential : Potential instance
        """
        self.e = e
        self.j = j
        f = lambda r: (e - potential(r * u.kpc) - j ** 2 / (r * u.kpc) ** 2 / 2)
        r = optimize.fsolve(f, 1)[0] * u.kpc
        Phi = potential(r)
        v = np.sqrt(2 *  (e - Phi))
        j_unit = u.kpc * u.km / u.s
        e_unit = u.km ** 2 / u.s ** 2
        assert np.isclose(j.to(j_unit).value, (r * v).to(j_unit).value)
        assert np.isclose(e.to(e_unit).value, (Phi + v ** 2 / 2).to(e_unit).value)
        coord0 = np.array([r.to(u.kpc).value, 0, 0, v.to(u.kpc / u.Gyr).value])
        # integrate for 1 Gyr
        t_int = 1
        t_dyn = (c.G * potential.M(r) * 3 / (4 * np.pi * r ** 3)) ** (-0.5)
        t_dyn = t_dyn.to(u.Gyr).value
        n_steps = int(t_int / (t_dyn / 100))
        t = np.linspace(0, t_int, n_steps)
        x, y, vx, vy = integrate.odeint(potential.D, coord0, t).T
        z = np.zeros(x.shape)
        vz = np.zeros(x.shape)
        pos = np.array([x, y, z])
        vel = np.array([vx, vy, vz])
        rotated_pos = Rz(Omega) @ Rx(i) @ pos
        rotated_vel = Rz(Omega) @ Rx(i) @ vel
        self.x, self.y, self.z = rotated_pos
        self.vx, self.vy, self.vz = rotated_vel
    
class OrbitLibrary:

    def __init__(self, ne, nj, nphi, ntheta):
        pass
    
