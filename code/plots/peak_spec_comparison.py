""" Compare the peak-light spectrum with other GRB-SNe """

import numpy as np
from astropy.io import ascii
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.cosmology import Planck15
from astropy.io import fits as pyfits
from scipy.optimize import curve_fit

fig,ax = plt.subplots(1,1,figsize=(7,4))


datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"
spec = np.loadtxt(datadir + "/ZTF20aalxlis_20200216_NOT_v1.ascii")
wl = spec[:,0]
f_lam = spec[:,1] # erg/cm2/s/AA, I think
hz = 3E10/(wl*1E-8)
f_nu = wl * f_lam / hz

plt.plot(
        wl/1.0252, f_lam/1E-15, c='k', 
        drawstyle='steps-mid', lw=0.5, ls='-', alpha=1, zorder=1)
plt.text(
    0.95, 0.95, "$\Delta t=12.5\,$d", fontsize=14, transform=ax.transAxes,
    horizontalalignment='right', verticalalignment='top')


# Indicate Ca II 
v = 30000
caii = np.array([8498,8542,8662])*(1.0252-v/3E5)
plt.scatter(caii, [0.3]*3, marker='|', c='k')
plt.text(caii[0]/1.05, 0.32, 'CaII (30,000 km/s)', fontsize=12,
        horizontalalignment='left', verticalalignment='bottom')

# Indicate Fe II
feii = 5169*(1.0252-v/3E5)
plt.scatter(feii, 0.4, marker='|', c='k')
plt.text(feii, 0.42, 'FeII', fontsize=12,
        horizontalalignment='center', verticalalignment='bottom')

# Indicate O I
# oi = 8446*(1.0252-v/3E5)
# plt.scatter(oi, 0.4, marker='|', c='k')
# plt.text(oi, 0.42, 'OI', fontsize=12,
#         horizontalalignment='center', verticalalignment='bottom')

# Indicate Si II 
si = 6355*(1.0252-v/3E5)
plt.scatter(si, 0.4, marker='|', c='k')
plt.text(si, 0.41, 'SiII', fontsize=12,
        horizontalalignment='center', verticalalignment='bottom')

# Plot the peak-light spectrum of SN2006aj for comparison
comp = ascii.read(datadir + "/sn2006aj-20060303-mmt.flm")
wl = comp['col1']
f_lam = comp['col2'] # erg/cm2/s/AA, I think
plt.plot(
        wl/1.0335, f_lam-0.5, c='#f6c746', 
        drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
plt.text(wl[-1]/1.0335, f_lam[-1]-0.5, '06aj at 12.85d', fontsize=12)

plt.ylabel(r"Flux ($10^{-15}$ erg\,s$^{-1}$\,cm$^{-2}\,$\AA$^{-1}$)", fontsize=16)
plt.xlabel("Rest Wavelength ($\AA$)", fontsize=16)
plt.ylim(-0.25, 0.66)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
#plt.show()
plt.savefig("spec_peaklight.png", dpi=200)