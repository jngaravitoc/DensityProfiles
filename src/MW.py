import numpy as np 
from astropy import units
import astropy as apy
from profiles import *
from math import erf


def acceleration(x, y, z):
    M_bulge = 1E10
    M_disk = 5.5E10
    M_halo = 1E12
    abulge = a_hernquist(0.7, x, y, z, M_bulge)
    adisk = a_mn(6.5, 0.6, x, y, z, M_disk)
    ahalo = a_NFW(11.0, x, y, z, M_halo)
    #print abulge, adisk, ahalo
    ax = abulge[0] + adisk[0] + ahalo[0]
    ay = abulge[1] + adisk[1] + ahalo[1]
    az = abulge[2] + adisk[2] + ahalo[2]
    ax = ax.to(units.kpc/units.Gyr**2) # / G.value / M_bulge
    ay = ay.to(units.kpc/units.Gyr**2) #/ G.value / M_bulge
    az = az.to(units.kpc/units.Gyr**2) #/ G.value / M_bulge
    return ax.value, ay.value, az.value

