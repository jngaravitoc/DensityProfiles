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
    plt.figure(figsize=(22, 11))
    plt.subplot(2, 3, 1)
    plt.plot(r, Pot, lw=2)
    plt.ylabel("$\Phi$", fontsize=30)
    plt.subplot(2, 3, 2)
    plt.plot(r, Density, lw=2)
    plt.ylabel(r"$\rho (M_{\odot} / Kpc^3)$", fontsize=30)
    plt.subplot(2, 3, 3)
    plt.plot(r, Mass, lw=2)
    plt.xlabel("$r(Kpc)$", fontsize=30)
    plt.ylabel("$M(M_{'odot})$", fontsize=30)
    plt.subplot(2, 3, 4)
    plt.plot(r, vc, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$v_c(km/s)$", fontsize=30)
    plt.subplot(2, 3, 5)
    plt.plot(r, a, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$a(km/s^2)$", fontsize=30)

def plots_hernquist(a, r, M):
    Pot = pot_hernquist(a, r, M)
    Density = dens_hernquist(a, r, M)
    Mass = mass_hernquist(a, r, M)
    vc = vc_hernquist(a, r, M)
    a = a_hernquist(a, r, M)
    plt.figure(figsize=(22, 11))
    plt.subplot(2, 3, 1)
    plt.plot(r, Pot, lw=2)
    plt.ylabel("$\Phi$", fontsize=30)
    plt.subplot(2, 3, 2)
    plt.plot(r, Density, lw=2)
    plt.ylabel(r"$\rho (M_{\odot} / Kpc^3)$", fontsize=30)
    plt.subplot(2, 3, 3)
    plt.plot(r, Mass, lw=2)
    plt.xlabel("$r(Kpc)$", fontsize=30)
    plt.ylabel("$M(M_{'odot})$", fontsize=30)
    plt.subplot(2, 3, 4)
    plt.plot(r, vc, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$v_c(km/s)$", fontsize=30)
    plt.subplot(2, 3, 5)
    plt.plot(r, a, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$a(km/s^2)$", fontsize=30)


def plots_sis(a, r, v):
    Pot = pot_sis(a, r, v)
    Density = dens_sis(a, r, v)
    Mass = mass_sis(a, r, v)
    vc = vc_sis(a, r, v)
    a = a_sis(a, r, v)
    plt.figure(figsize=(22, 11))
    plt.subplot(2, 3, 1)
    plt.plot(r, Pot, lw=2)
    plt.ylabel("$\Phi$", fontsize=30)
    plt.subplot(2, 3, 2)
    plt.plot(r, Density, lw=2)
    plt.ylabel(r"$\rho (M_{\odot} / Kpc^3)$", fontsize=30)
    plt.subplot(2, 3, 3)
    plt.plot(r, Mass, lw=2)
    plt.xlabel("$r(Kpc)$", fontsize=30)
    plt.ylabel("$M(M_{'odot})$", fontsize=30)
    plt.subplot(2, 3, 4)
    plt.plot(r, vc, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$v_c(km/s)$", fontsize=30)
    plt.subplot(2, 3, 5)
    plt.plot(r, a, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$a(km/s^2)$", fontsize=30)

def plots_mn(a, b, z, R, M):
    Pot = pot_mn(a, b, z, R, M)
    Density = dens_mn(a, b, z, R, M)
    Mass = mass_mn(a, b, z, R, M)
    vc = vc_mn(a, b, z, R, M)
    ar, az = a_mn(a, b, z, R, M)
    plt.figure(figsize=(22, 11))
    plt.subplot(2, 3, 1)
    plt.plot(R, Pot, lw=2)
    plt.ylabel("$\Phi$", fontsize=30)
    plt.subplot(2, 3, 2)
    plt.plot(R, Density, lw=2)
    plt.ylabel(r"$\rho (M_{\odot} / Kpc^3)$", fontsize=30)
    plt.subplot(2, 3, 3)
    plt.plot(R, Mass, lw=2)
    plt.xlabel("$r(Kpc)$", fontsize=30)
    plt.ylabel("$M(M_{'odot})$", fontsize=30)
    plt.subplot(2, 3, 4)
    plt.plot(R, vc, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$v_c(km/s)$", fontsize=30)
    plt.subplot(2, 3, 5)
    plt.plot(R, ar, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$a_R(km/s^2)$", fontsize=30)
    plt.subplot(2, 3, 6)
    plt.plot(R, az, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$a_z(km/s^2)$", fontsize=30)

def plots_log(Rc, q, z, R, v):
    Pot = pot_log(Rc, q, z, R, v)
    Density = dens_log(Rc, q, z, R, v)
    Mass = mass_log(Rc, q, z, R, v)
    vc = vc_log(Rc, q, z, R, v)
    ar, az  = a_log(Rc, q, z, R, v)
    plt.figure(figsize=(22, 11))
    plt.subplot(2, 3, 1)
    plt.plot(R, Pot, lw=2)
    plt.ylabel("$\Phi$", fontsize=30)
    plt.subplot(2, 3, 2)
    plt.plot(R, Density, lw=2)
    plt.ylabel(r"$\rho (M_{\odot} / Kpc^3)$", fontsize=30)
    plt.subplot(2, 3, 3)
    plt.plot(R, Mass, lw=2)
    plt.xlabel("$r(Kpc)$", fontsize=30)
    plt.ylabel("$M(M_{'odot})$", fontsize=30)
    plt.subplot(2, 3, 4)
    plt.plot(R, vc, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$v_c(km/s)$", fontsize=30)
    plt.subplot(2, 3, 5)
    plt.plot(R, ar, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$a_R(km/s^2)$", fontsize=30)
    plt.subplot(2, 3, 6)
    plt.plot(R, az, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$a_z(km/s^2)$", fontsize=30)

def plots_LMJ(r_h, q1, q2, qz, phi, x, y, z, v):
    Pot = pot_LMJ(r_h, q1, q2, qz, phi, x, y, z, v)
    Density = dens_log(Rc, q, z, R, v)
    Mass = mass_log(Rc, q, z, R, v)
    vc = vc_log(Rc, q, z, R, v)
    ar, az  = a_log(Rc, q, z, R, v)
    plt.figure(figsize=(22, 11))
    plt.subplot(2, 3, 1)
    plt.plot(R, Pot, lw=2)
    plt.ylabel("$\Phi$", fontsize=30)
    plt.subplot(2, 3, 2)
    plt.plot(R, Density, lw=2)
    plt.ylabel(r"$\rho (M_{\odot} / Kpc^3)$", fontsize=30)
    plt.subplot(2, 3, 3)
    plt.plot(R, Mass, lw=2)
    plt.xlabel("$r(Kpc)$", fontsize=30)
    plt.ylabel("$M(M_{'odot})$", fontsize=30)
    plt.subplot(2, 3, 4)
    plt.plot(R, vc, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$v_c(km/s)$", fontsize=30)
    plt.subplot(2, 3, 5)
    plt.plot(R, ar, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$a_R(km/s^2)$", fontsize=30)
    plt.subplot(2, 3, 6)
    plt.plot(R, az, lw=2)
    plt.xlabel("$r(kpc)$", fontsize=30)
    plt.ylabel("$a_z(km/s^2)$", fontsize=30)

