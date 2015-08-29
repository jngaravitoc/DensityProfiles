import numpy as np
from astropy import constants
from astropy import units

G = constants.G
K = constants.k_B
G = G.to(units.kiloparsec**3 / (units.Msun * units.s**2)) 
K = K.to(units.Msun * units.kpc**2 / (units.s**2 * units.Kelvin))


#++++++++++++++++ PLUMMER ++++++++++++++++++++++++++++++++++++++

def mass_plummer(a, r, M):
    a = a*units.kpc
    Mass = M*r**3 / (a**2 + r**2)**(3/2.)
    return Mass

def rho_plummer(a, r, M):
    a = a*units.kpc
    rho = 3*M / (4 *np.pi * a**3) * (1 + r**2/a**2)**(-5/2)
    return rho

def potential_plummer(a, r, M):
    a = a*units.kpc
    phi =  - G*M / np.sqrt(r**2 + a**2)
    return phi

def vc_plummer(a, r, M):
    a = a*units.kpc
    vc = np.sqrt(G*M*( r**2/(r**2 + a**2)**(3/2.)))
    vc = vc.to(units.km / units.s)
    return vc

#++++++++++++++++ HERNQUIST ++++++++++++++++++++++++++++

def Potential_Hernquist(a, r, M):
    a = a*units.kpc
    phi = -G*M / (r+a)
    return phi

def Density_Hernquist(a, r, M):
    a = a*units.kpc
    rho = M / (2 * np.pi) * a / (r*(r+a)**3)
    return rho

def Mass_Hernquist(a, r, M):
    a = a*units.kpc
    Mass = M * r**2 / (r+a)**2
    return Mass

def vc_Hernquist(a, r, M):
    a = a*units.kpc
    vc = np.sqrt(G*M*r/(r+a)**2)
    vc = vc.to(units.km / units.s)
    return vc

#+++++++++++++++++++ SIS (Singular Isothermal Sphere) ++++++++++++++++++++

def rho_sis(a, v, G, r):
    a = a*units.kpc
    v = v.to(units.kpc / units.s)
    rho = v**2 / (4*np.pi * G*(r**2 + a**2))
    return rho

def mass_sis(v, G, r):
    v = v.to(units.kpc / units.s)
    M = v**2 * r/G
    return M

def phi_sis(v, r):
    v = v.to(units.kpc / units.s)
    phi = v**2 * log(r.value)
    return phi*units.kpc
