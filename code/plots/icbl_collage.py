""" A collage of early light curves, for the introduction of the paper
I think this should be about GRB-SNe specifically
"""

import matplotlib.pyplot as plt
import numpy as np
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from astropy.io import ascii
from astropy.cosmology import Planck15
from astropy.time import Time
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/code")
from load_lc import get_forced_phot, get_lc


# I think I should focus on GRB-SNe here, and not get to the
# complications of early SN light curves.
# I'll also add in SN2018gep

def plot_980425(ax):
    """ Plot GRB 980425 early (t < 4 days) light curve """
    dm = Planck15.distmod(z=0.0085).value
    dat = ascii.read(
        "../../data/lc_980425.txt", format="fixed_width", delimiter="&")
    jd_raw = np.array(dat['JD'])
    delta = jd_raw[0] - 17/24 # first obs is 17 hrs after burst
    jd = jd_raw - delta
    band = 'Rc'
    lc = dat[band]
    bad = np.array(['nodata' in val for val in lc])
    mag = np.array([float(val.split("$\\pm$")[0]) for val in lc[~bad]])
    mag_err = np.array(
            [float(val.split("$\\pm$")[1]) for val in lc[~bad]])
    t = jd[~bad]
    ax.plot(t,mag-dm,c='k', linestyle='--')
    ax.errorbar(
            t,mag-dm,yerr=mag_err,ecolor='k', fmt='o',
            mfc='k',mec='k',label=band,ms=5)
    ax.text(0.01, 0.9, "GRB\,980425/SN\,98bw (Rc)", transform=ax.transAxes,
            fontsize=11)


def plot_2006aj(ax):
    """ Plot SN2006aj early (t < 4 days) light curve """
    dat = ascii.read("../../data/lc_060218_swift.dat")
    t = dat['t']
    dm = Planck15.distmod(z=0.033).value
    mag = dat['mag']
    emag = dat['emag']
    ax.errorbar(t, mag-dm, fmt='o', c='k')
    ax.plot(
        t, mag-dm, linestyle='--', lw=1, color='k')
    ax.text(
            1,0.1,"GRB\,060218/SN\,06aj (B)", 
            fontsize=11,transform=ax.transAxes,
            horizontalalignment='right')


def plot_171205A(ax):
    """ Plot from d'Elia 2018 """
    dm = Planck15.distmod(z=0.0368).value
    dat = ascii.read("../../data/lc_171205a.txt")
    t0 = Time("2017-12-05T07:20:43.9", format='isot')
    dt = dat['col1'].data / 86400
    mag = dat['col3'].data-dm
    mag_u = dat['col4'].data
    mag_l = dat['col5'].data
    filt = dat['col12'].data
    choose = filt == 'B'
    ax.errorbar(
            dt[choose], mag[choose], yerr=[-mag_l[choose], mag_u[choose]], 
            c='k', fmt='o')
    ax.text(
            1,0.9,"GRB\,171205A/SN\,17iuk (B)", 
            fontsize=11,transform=ax.transAxes,
            horizontalalignment='right')


def plot_18gep(ax):
    """ Plot from Ho 2019 """
    dm = Planck15.distmod(z=0.03154).value
    dt,filt,mag,emag = get_lc()
    choose = filt == 'g'
    ax.errorbar(dt[choose], mag[choose]-dm, emag[choose], fmt='o', c='k')
    ax.text(1.0, 0.1, 'SN\,18gep (g)', fontsize=11, transform=ax.transAxes,
    horizontalalignment='right')


#start, let's make a nice plot of two light curves: SN1998bw and SN2006aj
fig,axarr = plt.subplots(2,2,figsize=(7,4), sharex=True, sharey=False)
ax = axarr[0,0]
plot_980425(ax)
ax.set_xlim(-0.5,4)
ax.set_ylim(-16.8,-17.6)

ax = axarr[0,1]
plot_2006aj(ax)
ax.set_ylim(-16.2,-18.5)

ax = axarr[1,0]
plot_171205A(ax)
ax.set_ylim(-15.5,-18)

ax = axarr[1,1]
plot_18gep(ax)
ax.set_ylim(-15.5,-20.5)

for ax in axarr.flatten():
    ax.tick_params(axis='both', labelsize=12)

fig.text(0.5, 0.04, '$\Delta$ t (days)', ha='center', fontsize=14)
fig.text(0.04, 0.5, 'Abs Mag', va='center', rotation='vertical', fontsize=14)
fig.subplots_adjust(
        left=0.16, bottom=0.18, right=None, top=None, wspace=0.3, hspace=0)
plt.savefig("early_lc_collage.png", dpi=300, bbox_inches='tight')

#plt.show()
