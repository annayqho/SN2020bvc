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

fig,axarr = plt.subplots(2,2,figsize=(8,6), sharex=False)
dat = np.loadtxt("bol_lc.txt")
dt = dat[:,0]
lbol = dat[:,1]
lbol_hi = dat[:,2]
lbol_lo = dat[:,3]
teff = dat[:,4]
teff_hi = dat[:,5]
teff_lo = dat[:,6]
rph = dat[:,7]
rph_hi = dat[:,8]
rph_lo = dat[:,9]

# data directory
ddir = "/Users/annaho/Dropbox/astro/papers/papers_complete/SN2020bvc/data/bol_lc"

def llgrbs(ax):
    """ Plot bolometric LC of LLGRBs
    These all come from Cano et al. 2017 (scraped from Fig 9),
    except for the LC of SN2017iuk which comes from 
    """
    dat = np.loadtxt(ddir + "/sn1998bw_UBVRI_cano2013.dat", delimiter=',')
    dt = dat[:,0]
    lum = dat[:,1]
    ax.plot(dt, lum, c='Goldenrod', ls='-', lw=1, alpha=0.5, label="98bw")
    ax.text(26, 4E42, '98bw', fontsize=12)

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
    ax.text(26, 1.1E42, '06aj', fontsize=12)

    # 2017iuk
    lsol = 4E33
    dt = [0.003, 0.06, 0.17, 0.55, 1.04, 2.01, 7, 10.952, 14.936]
    lum = [2.1E45, 4.6E42, 3.9E42, 2.7E42, 2.7E42, 1.8E42, 10**9*lsol, 10**9.23*lsol, 10**9.22*lsol]
    ax.plot(dt, lum, c='#84206b', ls=':', lw=1, alpha=1, label='17iuk')
    ax.text(10, 7E42, '17iuk', fontsize=12)

    #ax.legend(loc='upper right', ncol=2)


def plot_lbol(ax):
    """ Plot the bolometric luminosity over time """
    ax.errorbar(dt/1.02507, lbol, yerr=[lbol_lo,lbol_hi], c='k', 
        fmt='o', mec='k', mfc='grey', ms=10)
    ax.set_ylabel(r'$L_\mathrm{bol}$ (erg/s)', fontsize=16)
    ax.set_ylim(1E42, 1E44)


# Luminosity panel
# linear
ax = axarr[0,0]
plot_lbol(ax)
llgrbs(ax)
ax.set_yscale('log')
ax.set_xlim(-1, 31)
ax.tick_params(axis='both', labelsize=16)

# Luminosity in log-log space without comparison
ax = axarr[0,1]
plot_lbol(ax)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(0.4,50)
ax.set_xticks([1,10])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.tick_params(axis='both', labelsize=16)

# Radius panel
ax = axarr[1,0]
xvals = np.linspace(0.5,20)
ax.errorbar(dt/1.02507, rph/1E14, yerr=[rph_lo/1E14,rph_hi/1E14], c='k', 
        fmt='o', mec='k', mfc='grey', ms=10)
yvals = 9E14 + 0.06 * (3E10) * xvals * 86400
ax.plot(xvals, yvals/1E14, ls='--', lw=0.5, c='grey')
ax.text(6, 29, '$v=0.06c$', fontsize=14, rotation=0)
ax.set_ylim(0,40)
ax.set_xlim(-1,31)
ax.set_ylabel(r'$R_\mathrm{ph}$ ($10^{14}$ cm)', fontsize=16)
ax.tick_params(axis='both', labelsize=16)
ax.set_xlabel(r'Rest-frame days after first light', fontsize=16)

# Temperature panel
ax = axarr[1,1]
ax.errorbar(dt/1.02507, teff, yerr=[teff_lo,teff_hi], c='k', fmt='o',
        mec='k', mfc='grey', ms=10)
ax.set_ylabel(r'$T_\mathrm{eff}$ (K)', fontsize=16)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(4000, 20000)
ax.axhline(y=5000, c='grey', ls='--', lw=0.5)
ax.text(0.5, 5000, "5000 K", fontsize=14, verticalalignment='bottom')
ax.set_xscale('log')
ax.set_xlim(0.4,50)
ax.set_xlabel(r'Rest-frame days after first light', fontsize=16)
ax.set_xticks([1,10])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.tick_params(axis='both', labelsize=16)

plt.tight_layout()
plt.show()
#plt.savefig("bb_evolution.eps", dpi=300)
