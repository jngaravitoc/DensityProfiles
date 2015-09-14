# History: 
# 11/09/15: 
# Implemented: Plummer, Hernquist, SIS, Miyamoto-Nagai
#########################################################


import numpy as np
from astropy import constants
from astropy import units

G = constants.G
K = constants.k_B
G = G.to(units.kiloparsec**3 / (units.Msun * units.s**2)) 
K = K.to(units.Msun * units.kpc**2 / (units.s**2 * units.Kelvin))


#++++++++++++++++ PLUMMER ++++++++++++++++++++++++++++++++++++++

def mass_plummer(a, r, M):
    M = M*units.Msun
    r = r*units.kpc
    a = a*units.kpc
    Mass = M*r**3 / (a**2 + r**2)**(3/2.)
    return Mass

def dens_plummer(a, r, M):
    M = M*units.Msun
    r = r*units.kpc
    a = a*units.kpc
    rho = 3*M / (4 *np.pi * a**3) * (1 + r**2/a**2)**(-5/2)
    return rho

def pot_plummer(a, r, M):
    M = M*units.Msun
    r = r*units.kpc
    a = a*units.kpc
    phi =  - G*M / np.sqrt(r**2 + a**2)
    return phi

def vc_plummer(a, r, M):
    a = a*units.kpc
    M = M*units.Msun
    r = r*units.kpc
    vc = np.sqrt(G*M*( r**2/(r**2 + a**2)**(3/2.)))
    vc = vc.to(units.km / units.s)
    return vc

def a_plummer(a, r, M):
    a = a*units.kpc
    M = M*units.Msun
    r = r*units.kpc
    A = - G * M * r / (r**2 + a**2)**(3/2.0)
    A = A.to(units.km / units.s**2)
    return A
	
#++++++++++++++++ HERNQUIST ++++++++++++++++++++++++++++

def pot_hernquist(a, r, M):
    a = a * units.kpc
    r = r * units.kpc
    M = M * units.Msun
    phi = -G*M / (r+a)
    return phi

def dens_hernquist(a, r, M):
    a = a * units.kpc
    r = r * units.kpc
    M = M * units.Msun
    rho = M / (2 * np.pi) * a / (r*(r+a)**3)
    return rho

def mass_hernquist(a, r, M):
    a = a*  units.kpc
    r = r * units.kpc
    M = M * units.Msun
    Mass = M * r**2 / (r+a)**2
    return Mass

def vc_hernquist(a, r, M):
    a = a*units.kpc
    r = r * units.kpc
    M = M * units.Msun
    vc = np.sqrt(G*M*r/(r+a)**2)
    vc = vc.to(units.km / units.s)
    return vc

def a_hernquist(a, r, M):
    a = a * units.kpc
    r = r * units.kpc
    M = M * units.Msun
    A =  - G * M  / (r + a)**2
    A = A.to(units.km / units.s**2) 
    return A

#+++++++++++++++++++ SIS (Singular Isothermal Sphere) ++++++++++++++++++++

def dens_sis(a, r, v):
    a = a * units.kpc
    r = r * units.kpc
    v = v * units.km / units.s
    v = v.to(units.kpc / units.s)
    rho = v**2 / (4 * np.pi * G * (r**2 + a**2))
    return rho

def mass_sis(a, r, v):
    a = a * units.kpc
    r = r * units.kpc
    v = v * units.km / units.s
    v = v.to(units.kpc / units.s)
    M = v**2 * r/G
    return M

def pot_sis(a, r, v):
    a = a * units.kpc
    r = r * units.kpc
    v = v * units.km / units.s
    v = v.to(units.kpc / units.s)
    phi = v**2 * np.log(r.value + a.value)
    return phi

def vc_sis(a, r, v):
    a = a * units.kpc
    r = r * units.kpc
    v = v * units.km / units.s
    v = v.to(units.kpc / units.s)
    V = v * np.sqrt(r + a) / np.sqrt(r)
    V = V.to(units.km / units.s )
    return V

def a_sis(a, r, v):
    a = a * units.kpc
    r = r * units.kpc
    v = v * units.km / units.s
    v = v.to(units.kpc / units.s)
    A = -v**2 / (r + a)
    A = A.to(units.km / units.s**2)
    return A


#+++++++++++++++++++++++ Miyamoto-Nagai +++++++++++++++++++++++++++

def pot_mn(a, b, z, r, M):
    z = z*units.kpc
    a = a*units.kpc
    b = b*units.kpc
    r = r*units.kpc
    phi = - G*M / (np.sqrt(r**2 + ( a + np.sqrt( z**2 + b**2 ))**2 ) )
    return phi

def dens_mn(a, b, z, r, M):
    z = z*units.kpc
    a = a*units.kpc
    b = b*units.kpc
    R = r*units.kpc
    rho = (b**2 * M / (4*np.pi)) * (a*R**2 + ( a + 3*(np.sqrt(z**2 + b**2)))*( a + np.sqrt(z**2 + b**2))**2 ) /( ( (R**2 + (a + np.sqrt(z**2 + b**2))**2)**(5./2.) * (z**2 + b**2)**(3./2.)) )
    return rho.value


def vc_mn(a, b, z, R, M):
    z = z*units.kpc
    a = a*units.kpc
    b = b*units.kpc
    R = r*units.kpc
    M = M * units.Msun
    vc = r*np.sqrt(G*M / ( (R**2 + (a + b)**2)**(3/2.0) ))
    vc = vc.to(units.km / units.s)	
    return vc

def mass_mn(a, b, z, R, M):
    z = z*units.kpc
    a = a*units.kpc
    b = b*units.kpc
    R = r*units.kpc
    M = M * units.Msun
    v = vc_mn(a, b, z, R, M)
    mass = v**2 * R / G
    return mass
    
def a_mn(a, b, z, R, M):
    z = z*units.kpc
    a = a*units.kpc
    b = b*units.kpc
    R = r*units.kpc
    M = M * units.Msun
    Ar = - G *  M * R / (R**2 + ( a + np.sqrt( z**2 + b**2))**2)**(3.0/2.0)
    Az = - G * M * z * (a + np.sqrt(z**2 + b**2)) / ( (R**2 + (a + np.sqrt(z**2 + b**2))**2)**(3.0/2.0) * np.sqrt(z**2 + b**2)  )
    Ar = Ar.to(units.km / units.s**2)
    Az = Az.to(units.km / units.s**2)
    return Ar, Az



