""" Plot the bolometric light curve 
and show the best-fit explosion parameters """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import sys
rc("font", family="serif")
rc("text", usetex=True)
from fit_arnett import *


# Load the bolometric light curve
dat = np.loadtxt("bol_lc.txt")
dt = dat[:,0]
lbol = dat[:,1]
lbol_hi = dat[:,2]
lbol_lo = dat[:,3]
elbol = []
# Make the uncertainty symmetric for curve_fit
for ii in np.arange(len(lbol_hi)):
    elbol.append(max(lbol_hi[ii], np.abs(lbol_lo[ii])))
elbol = np.array(elbol)

# Fit the nickel-decay model
popt, pcov = fit(dt, lbol, elbol, plot=True)
mni = popt[0]
tdiff = popt[1]

# Plot the light curve
fig,ax = plt.subplots(1,1,figsize=(6,4))
plt.errorbar(dt/1.02507, lbol, yerr=[lbol_lo, lbol_hi], c='k', fmt='o',
        mec='k', mfc='grey', ms=10)

# Shock-cooling model
Re = 3E12
ve = 0.1*3E10
te = Re/ve
tp = 0.3*86400
Me = 1E-2 * 2E33
Ee = Me*ve**2

tplot = np.linspace(0.01,40,1000)
t = tplot*86400
lplot_sc = ve*Re*Me/(4*t**2)

plt.plot(
        tplot, lplot_sc, lw=1, ls=':', c='#e55c30', 
        label="Post-shock cooling (SC)")

# Radioactive-decay model
tplot = np.linspace(0.01,40,1000)
lplot_rd = lph(tplot, mni, tdiff)
plt.plot(tplot, lplot_rd, lw=1, ls='--', c='#e55c30',
        label="Radioactive decay (RD)")

# Sum
lplot = lplot_sc+lplot_rd
plt.plot(tplot, lplot, lw=1, ls='-', c='k', label="SC+RD")

plt.xlabel("Rest-frame days since last non-detection", fontsize=16)
plt.ylabel("$L_\mathrm{bol}$ (erg s$^{-1}$)", fontsize=16)
plt.tick_params(axis='both', labelsize=16)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1E42,1E43)
plt.xlim(0.5,40)
plt.legend(loc='upper right',fontsize=12)
plt.tight_layout()
#plt.show()
plt.savefig("bol_lc.eps", format='eps', dpi=300)
