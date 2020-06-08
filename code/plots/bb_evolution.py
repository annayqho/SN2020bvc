""" Evolution of blackbody parameters, compared to other GRB-SNe """

import numpy as np
import matplotlib
from astropy.io import ascii
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.cosmology import Planck15
from astropy.table import Table

fig,axarr = plt.subplots(4,1,figsize=(4,9), sharex=False)
dt = np.array(
        [0.67, 1.36, 2.88, 3.81, 4.61, 6.27, 9.09, 10.86, 15.49, 26.51, 29.48])


def llgrbs(ax):
    """ Plot bolometric LC of LLGRBs
    These all come from Cano et al. 2017 (scraped from Fig 9),
    except for the LC of SN2017iuk which comes from 
    """
    ddir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/bol_lc"
    dat = np.loadtxt(ddir + "/sn1998bw_UBVRI_cano2013.dat", delimiter=',')
    dt = dat[:,0]
    lum = dat[:,1]
    ax.plot(dt, lum, c='Goldenrod', ls='-', lw=1, alpha=0.5, label="98bw")
    ax.text(30, 4E42, '98bw', fontsize=12)

    dat = np.loadtxt(ddir + "/sn2010bh_BVRI_cano2013.dat", delimiter=',')
    dt = dat[:,0]
    lum = dat[:,1]
    ax.scatter(
            dt, lum, edgecolor='grey', facecolor='white', 
            marker='D', lw=0.5, label='10bh')
    ax.text(10, 1.2E42, '10bh', fontsize=12)

    dat = np.loadtxt(ddir + "/sn2006aj_UBVRI_cano2013.dat", delimiter=',')
    dt = dat[:,0]
    lum = dat[:,1]
    dt = np.append(dt, 1.4)
    lum = np.append(lum, 2.6E44)
    order = np.argsort(dt)
    ax.plot(dt[order], lum[order], c='k', ls='-', lw=1, alpha=1, label='06aj')
    ax.text(30, 1.1E42, '06aj', fontsize=12)

    # 2017iuk
    lsol = 4E33
    dt = [0.003, 0.06, 0.17, 0.55, 1.04, 2.01, 7, 10.952, 14.936]
    lum = [2.1E45, 4.6E42, 3.9E42, 2.7E42, 2.7E42, 1.8E42, 10**9*lsol, 10**9.23*lsol, 10**9.22*lsol]
    ax.plot(dt, lum, c='#84206b', ls=':', lw=1, alpha=1, label='17iuk')
    ax.text(10, 8E42, '17iuk', fontsize=12)

    #ax.legend(loc='upper right', ncol=2)


def plot_lbol(ax):
    """ Plot the bolometric luminosity over time """
    lbol = np.array(
            [5.6, 4.97, 2.26, 1.98, 2.62, 2.95, 3.50, 3.67, 2.85, 2.19, 1.85])
    elbol = np.array(
            [0.2, 0.15, 0.12, 0.07, 0.15, 0.15, 0.13, 0.09, 0.08, 0.06, 0.05])
    ax.errorbar(dt/1.02507, lbol*1E42, yerr=elbol*1E42, c='k', fmt='o',
            mec='k', mfc='grey', ms=10)
    ax.set_ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)
    ax.set_ylim(1E42, 1E43)


# Luminosity panel
# linear
ax = axarr[0]
plot_lbol(ax)
llgrbs(ax)
ax.set_yscale('log')
ax.set_xlim(-1, 35)
ax.tick_params(axis='both', labelsize=16)
ax.arrow(
        0.67/1.02507, 5.6E42, 0, 2E42, 
        length_includes_head=True, head_length=5E41, head_width=1, color='k')

# Luminosity in log-log space without comparison
ax = axarr[1]
plot_lbol(ax)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(0.4,50)
ax.set_xticks([1,10])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.tick_params(axis='both', labelsize=16)
ax.arrow(
        0.67/1.02507, 5.6E42, 0, 2E42, 
        length_includes_head=True, head_length=5E41, head_width=0.1, color='k')

# Radius panel
ax = axarr[2]
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
ax.text(4, 18, 'v=0.04c', fontsize=14, rotation=0)
ax.set_ylim(0,25)
ax.set_xlim(-1,31)
ax.set_ylabel(r'$R_\mathrm{ph}$ ($10^{14}$ cm)', fontsize=16)
ax.scatter(1.4, 3.3, marker='x', zorder=5, c='k')
ax.text(
        2, 3.3, '06aj', fontsize=12, 
        verticalalignment='top', horizontalalignment='left')
ax.tick_params(axis='both', labelsize=16)

# Temperature panel
ax = axarr[3]
teff = 1000*np.array(
    [13.2,12.54,7.73,7.59,7.43,7.18,6.91,7.83,5.72,5.60,5.47])
eteff = np.array(
    [0.27,0.37,0.26,0.31,0.26,0.25,0.18,0.40,0.09,0.1,0.09])
ax.errorbar(dt/1.02507, teff, yerr=eteff, c='k', fmt='o',
        mec='k', mfc='grey', ms=10)
ax.set_ylabel(r'$T_\mathrm{eff}$ (K)', fontsize=16)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(4000, 15000)
ax.axhline(y=5000, c='grey', ls='--', lw=0.5)
ax.text(0.5, 5000, "5000 K", fontsize=14, verticalalignment='bottom')
ax.set_xscale('log')
ax.set_xlim(0.4,50)
ax.set_xlabel(r'Rest-frame days after first light', fontsize=16)
ax.set_xticks([1,10])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.tick_params(axis='both', labelsize=16)

plt.tight_layout()
#plt.show()
plt.savefig("bb_evolution.png", dpi=300)
