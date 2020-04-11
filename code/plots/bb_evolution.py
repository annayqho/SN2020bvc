""" Evolution of blackbody parameters, compared to other GRB-SNe """

import numpy as np
from astropy.io import ascii
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.cosmology import Planck15
from astropy.table import Table

fig,axarr = plt.subplots(3,1,figsize=(4,6), sharex=True)
dt = np.array([0.67, 1.36, 2.88, 3.81, 4.61, 6.27, 9.09, 10.86, 15.49, 26.51, 29.48])


def llgrbs(ax):
    """ Plot bolometric LC of LLGRBs """
    ddir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/bol_lc"
    dat = Table.read(ddir + "/sn1998bw.dat", format='ascii.fast_no_header')
    dt = dat['col1']
    lum = dat['col2']
    #ax.scatter(dt, lum, marker='.', c=col)
    ax.plot(dt, lum, c='k', ls='-', lw=1, alpha=0.5, label="Ic-BL")
    ax.text(dt[-1], lum[-1], '98bw', fontsize=10)

    dat = Table.read(ddir + "/sn2010bh.dat", format='ascii.fast_no_header')
    dt = dat['col1']
    lum = dat['col2']
    #ax.scatter(dt, lum, marker='.', c=col)
    ax.plot(dt, lum, c='k', ls='-', lw=1, alpha=0.5, label='_nolegend_')
    ax.text(dt[14], lum[14]/1.1, '10bh', fontsize=10, verticalalignment='top')

    dat = Table.read(ddir + "/sn2006aj.dat", format='ascii.fast_no_header')
    dt = dat['col1']
    order = np.argsort(dt)
    dt = dt[order]
    lum = dat['col2'][order]
    ax.plot(dt, lum, c='k', ls='-', lw=1, alpha=0.5, label='_nolegend_')
    ax.text(dt[-1], lum[-1], '06aj', fontsize=10)

    # 2017iuk
    dt = [1.5/24, 0.975, 1.947, 7, 10.952, 14.936]
    lum = [1.7E45, 4.2E42, 3.7E42, 2.4E42, 2.4E42, 1.7E42]
    ax.plot(dt, lum, c='k', ls='-', lw=1, alpha=0.5, label='_nolegend_')
    ax.text(dt[-1], lum[-1], '17iuk', fontsize=10, verticalalignment='top')


# Luminosity panel
ax = axarr[0]
lbol = np.array(
        [5.6, 4.97, 2.26, 1.98, 2.62, 2.95, 3.50, 3.67, 2.85, 2.19, 1.85])
elbol = np.array(
        [0.2, 0.15, 0.12, 0.07, 0.15, 0.15, 0.13, 0.09, 0.08, 0.06, 0.05])
ax.errorbar(dt/1.02507, lbol*1E42, yerr=elbol*1E42, c='k', fmt='o',
        mec='k', mfc='grey', ms=10)
ax.arrow(
        0.67/1.02507, 5.6E42, 0, 2E42, 
        length_includes_head=True, head_length=5E41, head_width=1, color='k')
llgrbs(ax)
ax.set_yscale('log')
ax.set_ylim(1E42, 1E43)
ax.set_ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)


# Radius panel
ax = axarr[1]
rph = np.array(
        [5.1, 4.61, 8.25, 8.60, 9.80, 11.49, 
        14.18, 12.06, 19.32, 17.67, 17.05])
erph = np.array([0.1,0.26,0.55,0.68,0.67,1.06,0.86,1.06,0.78,0.59,0.54])
xvals = np.linspace(0.5,20)
ax.errorbar(dt/1.02507, rph, yerr=erph, c='k', fmt='o',
        mec='k', mfc='grey', ms=10)
ax.arrow(
        0.67/1.02507, 5.1, 0, 1, 
        length_includes_head=True, head_length=1, head_width=0.1, color='k')
yvals = 5E14 + 0.04 * (3E10) * xvals * 86400
ax.plot(xvals, yvals/1E14, ls='--', lw=0.5, c='grey')
#ax.text(26, 6.8, 'v=0.1c', fontsize=14, rotation=20)
ax.text(5, 18, 'v=0.04c', fontsize=14, rotation=0)
ax.set_ylim(1,21)

ax.set_ylabel(r'$R_\mathrm{ph}$ ($10^{14}$ cm)', fontsize=16)

# Temperature panel
ax = axarr[2]
teff = 1000*np.array(
    [13.2,12.54,7.73,7.59,7.43,7.18,6.91,7.83,5.72,5.60,5.47])
eteff = np.array(
    [0.27,0.37,0.26,0.31,0.26,0.25,0.18,0.40,0.09,0.1,0.09])
ax.errorbar(dt/1.02507, teff, yerr=eteff, c='k', fmt='o',
        mec='k', mfc='grey', ms=10)
axarr[2].set_ylabel(r'$T_\mathrm{eff}$ (K)', fontsize=16)
ax.set_yscale('log')
ax.set_ylim(4000, 15000)
ax.axhline(y=5000, c='grey', ls='--', lw=0.5)
ax.text(0.5, 4600, "5000 K", fontsize=14, verticalalignment='top')

axarr[0].xaxis.label.set_visible(False)
axarr[1].xaxis.label.set_visible(False)


axarr[0].tick_params(axis='y', labelsize=16)
axarr[1].tick_params(axis='y', labelsize=16)
axarr[2].tick_params(axis='y', labelsize=16)
axarr[2].tick_params(axis='x', labelsize=16)

axarr[2].set_xlabel(r'Rest-frame days after first light', fontsize=16)

plt.tight_layout()
plt.show()
#plt.savefig("bb_evolution.png", dpi=300)
