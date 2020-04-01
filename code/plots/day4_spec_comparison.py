""" Fit blackbody SED but also spectrum on the first day """

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
spec = np.loadtxt(datadir + "/ZTF20aalxlis_20200207_P60_v1.ascii")
wl = spec[:,0]
f_lam = spec[:,1] # erg/cm2/s/AA, I think
ef_lam = np.sqrt(spec[:,2])
hz = 3E10/(wl*1E-8)
f_nu = wl * f_lam / hz
ef_nu = wl * ef_lam / hz

plt.plot(wl/1.0252, f_lam/1E-15, c='k', drawstyle='steps-mid', lw=0.5, ls='-', alpha=1)
plt.text(
    0.95, 0.95, "$\Delta t=3.7\,$d", fontsize=14, transform=ax.transAxes,
    horizontalalignment='right', verticalalignment='top')

# Indicate Ca II at v=60,000 km/s
v = 51000
caii = np.array([8498,8542,8662])*(1.0252-v/3E5)
plt.scatter(caii, [0.3]*3, marker='|', c='k')
plt.text(caii[0]/1.05, 0.32, 'CaII (51,000 km/s)', fontsize=12,
        horizontalalignment='left', verticalalignment='bottom')


# Indicate Si II at v=higher velocity...
# v = 80000
# si = 6355*(1.0252-v/3E5)
# plt.scatter(si, 0.8, marker='|', c='k')
# plt.text(si/1.05, 0.81, 'SiII (80,000 km/s)', fontsize=12,
#         horizontalalignment='left', verticalalignment='bottom')
# 

# Add a 4d spectrum of 06aj for comparison
comp = ascii.read(datadir + "/sn2006aj-20060222-fast.flm")
wl = comp['col1']
f_lam = comp['col2'] # erg/cm2/s/AA, I think
plt.plot(
        wl/1.0335, f_lam, c='#f6c746',
        drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
plt.text(wl[-1]/1.0335, f_lam[-1], '06aj at 3.97d', fontsize=12)


plt.ylabel(r"Flux ($10^{-15}$ erg\,s$^{-1}$\,cm$^{-2}\,$\AA$^{-1}$)", fontsize=16)
plt.xlabel("Rest Wavelength ($\AA$)", fontsize=16)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
#plt.show()
plt.savefig("spec_day4.png", dpi=200)
