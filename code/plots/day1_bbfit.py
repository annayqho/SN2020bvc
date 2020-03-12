""" Fit blackbody SED but also spectrum on the first day """

import numpy as np
from astropy.io import ascii
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.cosmology import Planck15
from scipy.optimize import curve_fit

fig,ax = plt.subplots(1,1,figsize=(7,4))

def bb_func(nu,T,R):
    d = Planck15.luminosity_distance(z=0.02507).cgs.value
    h = 6.63E-27
    c = 3E10
    k = 1.38E-16
    Inu = (2*h*nu**3/c**2) * (1/(np.exp(h*nu/(k*T)-1)))
    fnu = Inu * np.pi * R**2 / d**2
    return fnu 


datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"
spec = np.loadtxt(datadir + "/ZTF20aalxlis_20200204_P60_v1.ascii")
wl = spec[:,0]
f_lam = spec[:,1] # erg/cm2/s/AA, I think
ef_lam = spec[:,2]
hz = 3E10/(wl*1E-8)
f_nu = wl * f_lam / hz
popt, pcov = curve_fit(bb_func, hz, f_nu, p0=[5000,1E14])
xplt = np.linspace(3500,9500)
xplt_nu = 3E10/(xplt*1E-8)
yplt_nu = bb_func(xplt_nu, popt[0], popt[1])
yplt = xplt_nu * yplt_nu / xplt

T = str(np.round(popt[0]/1E3,1))
R = str(np.round(popt[1]/1E14,1))
L = str(np.round((4 * np.pi * popt[1]**2 * 5.67E-5 * popt[0]**4)/1E42,1))

plt.plot(wl, f_lam/1E-15, c='k', drawstyle='steps-mid', lw=0.5, ls='-', alpha=1)
plt.plot(xplt,yplt/1E-15,lw=1,alpha=1,c='k')
plt.text(
    0.95, 0.95, "$\Delta t=0.67\,$d", fontsize=14, transform=ax.transAxes,
    horizontalalignment='right', verticalalignment='top')
plt.text(
    0.95, 0.85, 
    r"$T=%s \times 10^{3} $\,K, $R=%s\times10^{14}$ cm, $L=%s\times10^{42}\,$erg/sec" %(T,R,L),
    horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,
    fontsize=14)

plt.ylabel(r"Flux ($10^{-15}$ erg\,s$^{-1}$\,cm$^{-2}\,$\AA$^{-1}$)", fontsize=16)
plt.xlabel("Wavelength (Angstrom)", fontsize=16)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
#plt.show()
plt.savefig("bb_first_epoch.png", dpi=200)
