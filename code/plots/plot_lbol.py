""" Plot the bolometric light curve 
and show the best-fit explosion parameters """
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import sys
rc("font", family="serif")
rc("text", usetex=True)
sys.path.append("/Users/annaho/Dropbox/Projects/Research/SN2020bvc/code")
from arnett import lph,lph_khatami


fig,ax = plt.subplots(1,1,figsize=(6,4))
dt = np.array(
        [0.67, 1.36, 2.88, 3.81, 4.61, 6.27, 9.09, 10.86, 15.49, 26.51, 29.48])
lbol = np.array(
        [5.6, 4.97, 2.26, 1.98, 2.62, 2.95, 3.50, 3.67, 2.85, 2.19, 1.85])
elbol = np.array(
        [0.2, 0.15, 0.12, 0.07, 0.15, 0.15, 0.13, 0.09, 0.08, 0.06, 0.05])
rph = np.array(
        [5.1, 4.61, 8.25, 8.60, 9.80, 11.49, 14.18, 12.06, 19.32, 17.67, 17.05])
plt.errorbar(dt/1.02507, lbol*1E42, yerr=elbol*1E42, c='k', fmt='o',
        mec='k', mfc='grey', ms=10)

# make the first measurement a lower limit
plt.arrow(
        0.67/1.02507, 5.6E42, 0, 2E42, 
        length_includes_head=True, head_length=5E41, head_width=0.1, color='k')

# Shock-cooling model
re = 1E13
ve = 0.1*3E10
te = re/ve
tp = 0.6*86400
Re = 4E12
Me = 1E-2 * 2E33
Ee = Me*ve**2

tplot = np.linspace(0.01,40,1000)
t = tplot*86400
lplot_sc = ve*Re*Me/(4*t**2)

#lplot_sc = (te*Ee/(tp**2))*np.exp(-t*(t+2*te)/(2*tp**2))

plt.plot(
        tplot, lplot_sc, lw=1, ls=':', c='#e55c30', 
        label="Post-shock cooling (SC)")

# Radioactive-decay model
tplot = np.linspace(0.01,40,1000)
lplot_rd_kh = np.array([lph_khatami(tval) for tval in tplot])
plt.plot(tplot, lplot_rd_kh, lw=1, ls='--', c='#e55c30',
        label="Radioactive decay (RD)")

# Sum
lplot = lplot_sc+lplot_rd_kh
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
plt.show()
#plt.savefig("bol_lc.png", dpi=300)
