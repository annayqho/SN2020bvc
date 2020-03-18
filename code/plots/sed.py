""" 
Plot the spectrum from radio to X-ray frequencies 

Choose Day 13 ish
because it has the best sampling.
"""

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import numpy as np
from astropy.cosmology import Planck15
import sys
from astropy.modeling.blackbody import blackbody_nu



bigsize=16
smallsize=14

d = Planck15.luminosity_distance(z=0.02507).cgs.value


def mjy_to_lum(f):
    """ Convert flux from mJy to erg/s/Hz """
    lum = f * 1e-3 * 1e-23 * 4 * np.pi * d**2
    return lum


def plot_xray_spindex(x, y, col):
    """ Dashed line indicating X-ray spectral index """
    Gamma = 1.54
    # f_nu \propto nu^{-alpha}
    # Gamma = photon index = alpha + 1
    # 0.3 -- 10 keV
    # 7.2E16 Hz -- 2.4E18 Hz
    xvals = np.linspace(2.4E16, 2.4E18)
    yvals = y * (xvals/x)**(0.46)
    plt.plot(xvals, yvals, ls='-', c=col)


def plot_xray(lum):
    """ OK so the approach is to take the integrated lum across 0.5-7 keV,
    L = int_nu1^nu2 L_nu dnu
    and L_nu = A nu^{-1}

    In this case nu^{-1}, so that means nu L nu ~ constant

    Parameters
    ----------
    x: frequency
    y: flux
    c: color of the point/line

    0.5 keV = 1.21E17 Hz
    7 keV = 1.69E18 Hz
    """
    nu2 = 1.69E18
    nu1 = 1.21E17
    alpha = 1
    xplot = np.linspace(nu1, nu2, 1000)
    yplot = [lum]*1000
    plt.plot(xplot, yplot, c='k', ls='-')

# set up the plot
fig = plt.figure(figsize=(7,4.7), dpi=100)

# Plot the optical blackbody 
R = 12E14
T = 7.83E3
c = 3E18
low_wl = 1928 # Angstroms, UVW2
hi_wl = 21900 # Angstroms, K-band
x = np.logspace(np.log10(c/hi_wl), np.log10(c/low_wl))
y = blackbody_nu(x, T)
plt.plot(x, 4*np.pi**2*R**2*x*y, c='k')

# Plot the radio point
x = 10E9
y = x*mjy_to_lum(66E-3)
plt.scatter(
        x, y, edgecolor='k', facecolor='k', marker='D', lw=0.5, s=50)

# Plot the Chandra data
# geometric mean of the frequencies
x = np.sqrt(1.21E17*1.69E18)
y = 1E40
plt.scatter(x, y, edgecolor='k', facecolor='k', marker='D', lw=0.5, s=50)
#plt.text(2E17,3E39,r'$\nu^{-1}$', fontsize=smallsize)

# Power-law with index 0.22...
nuplot = np.linspace(1E10, 1E18)
# for nu^{-0.78}:
#lplot = (1E37) * (nuplot/1E10)**(0.22)
# for nu^{-0.5} = -(p-1)/2 where p=2
lplot = (1E37) * (nuplot/1E10)**(0.5)
plt.plot(nuplot,lplot,lw=0.5,ls='--',c='k')
plt.text(5E11,2E38,r'$L_\nu \propto \nu^{-0.5}$',fontsize=14)
lplot = (1E37) * (nuplot/1E10)**(0.25)
plt.plot(nuplot,lplot,lw=0.5,ls='--',c='k')
plt.text(1E14,4E37,r'$L_\nu \propto \nu^{-0.75}$',fontsize=14)

# Model for Epoch 2 is a power-law:
# norm = 9E-4 # photons/keV/cm2/s at 1 keV
# f = norm * 10**(-12)  # maybe in units 1E-12 erg/cm2/s
# alpha = 1.39
# nu = np.linspace(7.25E17, 1.93E19, 1000)
# energy = nu * 6.626E-27
# y = f * energy**(-alpha) * 4 * np.pi * d**2
#plt.plot(nu, y, ls=':', c='k')

plt.xlim(1E9, 5E18)
plt.ylim(1E36, 5E43)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Frequency [Hz]", fontsize=bigsize)
plt.ylabel(r"$\nu L_\nu$ (erg/s)", fontsize=bigsize)
plt.xticks(fontsize=bigsize)
plt.yticks(fontsize=bigsize)

plt.tight_layout()

plt.savefig("sed.png", dpi=300)
plt.close()
#plt.show()
