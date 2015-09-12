import numpy as np 
import matplotlib.pyplot as plt
from profiles import *
font = {'size':15, 'family':'serif'}
plt.matplotlib.rc('font', **font)

def plots_plummer(a, r, M):
    Pot = pot_plummer(a, r, M)
    Density = dens_plummer(a, r, M)
    Mass = mass_plummer(a, r, M)
    vc = vc_plummer(a, r, M)
    a = a_plummer(a, r, M)
    plt.figure(figsize=(18, 9))
    plt.subplot(2, 3, 1)
    plt.plot(r, Pot, lw=2)
    plt.ylabel("$\Phi$", fontsize=20)
    plt.subplot(2, 3, 2)
    plt.plot(r, Density, lw=2)
    plt.ylabel(r"$\rho (M_{\odot} / Kpc^3)$", fontsize=20)
    plt.subplot(2, 3, 3)
    plt.plot(r, Mass, lw=2)
    plt.xlabel("$r(Kpc)$", fontsize=20)
    plt.ylabel("$M(M_{'odot})$", fontsize=20)
    plt.subplot(2, 3, 4)
    plt.plot(r, vc, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=20)
    plt.ylabel("$v_c(km/s)$", fontsize=20)
    plt.subplot(2, 3, 5)
    plt.plot(r, a, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=20)
    plt.ylabel("$a(km/s^2)$", fontsize=20)


