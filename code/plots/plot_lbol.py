""" Plot the bolometric light curve 
and show the best-fit explosion parameters """
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)


fig,ax = plt.subplots(1,1,figsize=(5,3))
dt = np.array(
        [0.67, 1.36, 2.88, 3.81, 4.61, 6.27, 9.09, 10.86, 15.49, 26.51, 29.48])
lbol = np.array(
        [5.6, 4.97, 2.26, 1.98, 2.62, 2.95, 3.50, 3.67, 2.85, 2.19, 1.85])
elbol = np.array(
        [0.2, 0.15, 0.12, 0.07, 0.15, 0.15, 0.13, 0.09, 0.08, 0.06, 0.05])
rph = np.array(
        [5.1, 4.61, 8.25, 8.60, 9.80, 11.49, 14.18, 12.06, 19.32, 17.67, 17.05])
#tplot = np.linspace(0.7, 17)
#yplot = tplot[1]*86400 + 0.05*3E10*tplot*86400
#plt.plot(tplot, yplot/1E14, c='k', lw=0.1)

#tplot = np.linspace(0.7, 3)
#yplot = l(tplot*86400, 100, 1E52, 0.67*86400)
#plt.plot(tplot, yplot/1E42, c='k', lw=0.1)

plt.errorbar(dt/1.02507, lbol*1E42, yerr=elbol*1E42, c='k', fmt='o')

# make the first measurement a lower limit
plt.arrow(
        0.67/1.02507, 5.6E42, 0, 1E42, 
        length_includes_head=True, head_length=5E41, head_width=0.1, color='k')

t0 = 1.36
L0 = 4.97E42
tplot = np.linspace(0.5,4)
lplot = L0 * (tplot/t0)**(-2)
plt.plot(tplot, lplot, lw=0.5, ls='--', c='#e55c30')

plt.xlabel("Rest-frame days since last non-detection", fontsize=14)
plt.ylabel("$L_\mathrm{bol}$ (erg s$^{-1}$)", fontsize=14)
plt.tick_params(axis='both', labelsize=14)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1E42,1E43)
plt.tight_layout()
#plt.show()
plt.savefig("bol_lc.png", dpi=200)
