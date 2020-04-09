""" Fit blackbody SED but also spectrum on the first day """

import numpy as np
from astropy.io import ascii
from astropy.io import fits
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

plt.plot(wl/1.0252, f_lam, c='k', drawstyle='steps-mid', lw=0.5, ls='-', alpha=1)
plt.text(
    0.95, 0.95, "20bvc at $\Delta t=3.7\,$d", fontsize=14, transform=ax.transAxes,
    horizontalalignment='right', verticalalignment='top')

# Indicate Ca II at v=60,000 km/s
v = 51000
caii = np.array([8498,8542,8662])*(1-v/3E5)
plt.scatter(caii, np.array([0.3]*3)*1E-15, marker='|', c='k')
plt.text(caii[0]/1.05, 0.32*1E-15, 'CaII (51,000 km/s)', fontsize=12,
        horizontalalignment='left', verticalalignment='bottom')

# Add a spectrum of 17iuk for comparison
comp = fits.open(datadir + "/GTC_OSIRIS_GRB171205A/1D_GRB171205A_171207_R1000B.fits")[0].data
wl = np.arange(3629.598, 3629.598+len(comp)*2.071432195122, 2.071432195122)
plt.plot(
        wl/1.0368, comp, c='#84206b',
        drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
plt.text(wl[-1]/1.03, comp[-1], '17iuk (2d)', fontsize=12)

# Indicate Si II
v = 75000
si = 6355*(1-v/3E5)
plt.scatter(si, 8E-17, marker='|', c='k')
plt.text(si/1.05, 8E-17, 'SiII (75,000 km/s)', fontsize=12,
        horizontalalignment='left', verticalalignment='bottom')

# Indicate Ca II at v=65,000 km/s
v = 105000
caii = np.array([8498,8542,8662])*(1-v/3E5)
plt.scatter(caii, [5.7E-17]*len(caii), marker='|', c='k')
plt.text(caii[0], 5.7E-17, 'CaII (105,000 km/s)', fontsize=12,
        horizontalalignment='left', verticalalignment='bottom')

# Add a 4d spectrum of 06aj for comparison
comp = ascii.read(datadir + "/sn2006aj-20060222-fast.flm")
wl = comp['col1']
f_lam = comp['col2'] # erg/cm2/s/AA, I think
plt.plot(
        wl/1.0335, f_lam*1E-15, c='#f6c746',
        drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
plt.text(wl[-1]/1.0335, f_lam[-1]*1E-15, '06aj (4d)', fontsize=12)


plt.yscale('log')
plt.ylim(1E-17, 2E-15)
plt.ylabel(r"Flux (erg\,s$^{-1}$\,cm$^{-2}\,$\AA$^{-1}$)", fontsize=16)
plt.xlabel("Rest Wavelength ($\AA$)", fontsize=16)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
#plt.show()
plt.savefig("spec_day4.png", dpi=200)
