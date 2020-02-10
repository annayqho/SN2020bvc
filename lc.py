""" Plot the light curve of this SN and compare it to SN 2006aj """

import numpy as np
from astropy.io import ascii
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.cosmology import Planck15

dm = Planck15.distmod(z=0.033).value-Planck15.distmod(z=0.025235).value

# SN 2006aj
dat = ascii.read("lc_060218_swift.dat")
t = dat['t']
mag = dat['mag']
emag = dat['emag']
plt.plot(t, mag-dm, linestyle='-', lw=1, color='k', label="2006aj (UVOT/B)")

# Our SN
t0 = 2458882.5
t = np.array([2458883.9763, 2458884.9212, 2458885.9290, 2458886.9768, 2458886.9941, 2458887.9290])
mag = np.array([16.90, 17.48, 17.48, 17.37, 17.37, 17.15])
emag = np.array([0.04, 0.06, 0.05, 0.04, 0.05, 0.04])
plt.errorbar(t-t0, mag, yerr=emag, fmt='o', c='black', label="SN2020bvc (ZTF/g)")
plt.scatter(2458881.9944-t0, 20.77, marker='o', c='black')
plt.arrow(
        2458881.9944-t0, 20.77, 0, 0.4, length_includes_head=True, 
        head_length=0.2, head_width=0.2, color='black')

# Visibility to GBM
times = Time(np.linspace(2458881.9944, 2458881.9944+1.0, 20), format='jd')
plt.text(2458881.9944+1.0-t0, 20, "Fermi/GBM visibility", horizontalalignment='left')
vis = [1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1]
for ii,time in enumerate(times):
    if vis[ii] == 0:
        plt.axvline(x=time.jd-t0, lw=2.5, c='red', alpha=0.2, zorder=0)
    else:
        plt.axvline(x=time.jd-t0, lw=2.5, c='blue', alpha=0.2, zorder=0)


plt.legend(fontsize=12)
plt.gca().invert_yaxis()
plt.tick_params(labelsize=14, axis='both')
plt.xlabel("Days since GRB 060218", fontsize=16)
plt.ylabel("Apparent Mag ($z=0.025$)", fontsize=16)
plt.ylim(21.5, 16.5)
plt.tight_layout()
#plt.show()
plt.savefig("lc.png", dpi=200)
