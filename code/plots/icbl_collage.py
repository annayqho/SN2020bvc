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
    ax.text(0.01, 0.8, "$z=0.0085$ (40 Mpc)", transform=ax.transAxes,
            fontsize=11)


def plot_2006aj(ax):
    """ Plot SN2006aj early (t < 4 days) light curve """
    dat = ascii.read("../../data/lc_060218_swift.dat")
    t = dat['t']
    dm = Planck15.distmod(z=0.033).value
    mag = dat['mag']
    emag = dat['emag']
    ax.errorbar(t, mag-dm, fmt='o', c='k', ms=5)
    ax.plot(
        t, mag-dm, linestyle='--', lw=1, color='k')
    ax.text(
            1,0.9,"GRB\,060218/SN\,06aj (B)", 
            fontsize=11,transform=ax.transAxes,
            horizontalalignment='right')
    ax.text(
            1.0,0.8,"$z=0.033$ (150 Mpc)", 
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
    ax.plot(
            dt[choose], mag[choose], c='k', ls='--')
    ax.text(
            1,0.9,"GRB\,171205A/SN\,17iuk (B)", 
            fontsize=11,transform=ax.transAxes,
            horizontalalignment='right')
    ax.text(
            1,0.8,"$z=0.0368$ (170 Mpc)", 
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


def plot_14gqr(ax):
    """ Data from De 2018 """
    dm = Planck15.distmod(z=0.063).value
    dt = np.array([0.42, 0.56, 0.63, 1.27, 1.34, 2.32, 3.25, 5.26])
    mag = np.array([19.95, 20.04, 20.12, 20.52, 20.47, 20.36, 20.17, 19.85])
    emag = np.array([0.05, 0.04, 0.05, 0.11, 0.10, 0.04, 0.07, 0.02])
    ax.errorbar(dt, mag-dm, emag, fmt='o', c='k')
    ax.plot(dt, mag-dm, emag, ls='--', c='k')
    ax.text(
            0.95, 0.9, 'Type Ic iPTF\,14gqr (g)', 
            fontsize=11, transform=ax.transAxes,
    horizontalalignment='right')
    ax.text(0.95, 0.8, '$z=0.063$', fontsize=11, transform=ax.transAxes,
    horizontalalignment='right')


def plot_18xas(ax):
    """ Data from Fremling 2019 """
    dm = Planck15.distmod(z=0.058832).value
    dt = np.array([0.42, 0.56, 0.63, 1.27, 1.34, 2.32, 3.25, 5.26])
    mag = np.array([19.95, 20.04, 20.12, 20.52, 20.47, 20.36, 20.17, 19.85])
    emag = np.array([0.05, 0.04, 0.05, 0.11, 0.10, 0.04, 0.07, 0.02])
    ax.errorbar(dt, mag-dm, emag, fmt='o', c='k')
    ax.plot(dt, mag-dm, emag, ls='--', c='k')
    ax.text(
            0.95, 0.9, 'Type IIb ZTF18aalrxas (g)', 
            fontsize=11, transform=ax.transAxes,
    horizontalalignment='right')
    ax.text(
            0.95, 0.8, '$z=0.063$ (292 Mpc)', 
            fontsize=11, transform=ax.transAxes,
            horizontalalignment='right')


def plot_2011dh(ax):
    """ Data from """
    dm = Planck15.distmod(z=0.001638).value
    dat = ascii.read("sn2011dh.txt")
    t = dat['time']
    filt = dat['band']
    mag = dat['magnitude']
    emag = dat['e_magnitude']
    choose = filt == 'V'
    dt = t[choose]-t[choose][0]
    mag = mag[choose]
    emag = emag[choose]
    ax.errorbar(dt, mag-dm, emag, fmt='o', c='k')
    ax.plot(dt, mag-dm, emag, ls='--', c='k')
    ax.text(
            0.95, 0.9, 'Type IIb SN2011dh (V)', 
            fontsize=11, transform=ax.transAxes,
    horizontalalignment='right')
    ax.text(0.95, 0.8, '$z=0.0016$ (7.3 Mpc)', fontsize=11, transform=ax.transAxes,
    horizontalalignment='right')


def plot_2011dh(ax):
    """ Data from """
    dm = Planck15.distmod(z=0.001638).value
    dat = ascii.read("sn2011dh.txt")
    t = dat['time']
    filt = dat['band']
    mag = dat['magnitude']
    emag = dat['e_magnitude']
    choose = filt == 'V'
    dt = t[choose]-t[choose][0]
    mag = mag[choose]
    emag = emag[choose]
    ax.errorbar(dt, mag-dm, emag, fmt='o', c='k')
    ax.plot(dt, mag-dm, emag, ls='--', c='k')
    ax.text(
            0.95, 0.9, 'Type IIb SN2011dh (V)', 
            fontsize=11, transform=ax.transAxes,
    horizontalalignment='right')
    ax.text(0.95, 0.8, '$z=0.0016$ (7.3 Mpc)', fontsize=11, transform=ax.transAxes,
    horizontalalignment='right')


def plot_15dtg(ax):
    """ Data from """
    dm = Planck15.distmod(z=0.0524).value
    t = np.array([333.931, 334.931, 335.931, 337.948, 337.960, 338.584])
    mag = np.array([19.634, 19.720, 19.943, 20.383, 20.412, 20.288])
    emag = np.array([0.162, 0.029, 0.027, 0.082, 0.278, 0.109])
    dt = (t-t[0]-0.4)/(1.0524)
    ax.errorbar(dt, mag-dm, emag, fmt='o', c='k')
    ax.plot(dt, mag-dm, emag, ls='--', c='k')
    ax.text(
            0.95, 0.9, 'Type Ic iPTF15dtg (V)', 
            fontsize=11, transform=ax.transAxes,
    horizontalalignment='right')
    ax.text(0.95, 0.8, '$z=0.0524$ (241 Mpc)', fontsize=11, transform=ax.transAxes,
    horizontalalignment='right')


def top_panel():
    """ GRB-SNe """
    fig,axarr = plt.subplots(1,3,figsize=(11,3), sharex=True, sharey=False)
    ax = axarr[0]
    plot_980425(ax)
    ax.set_xlim(-0.5,4)
    ax.set_ylim(-16.8,-17.6)

    ax = axarr[1]
    plot_2006aj(ax)
    ax.set_ylim(-16.2,-18.5)

    ax = axarr[2]
    plot_171205A(ax)
    ax.set_ylim(-15.5,-18)

    #ax = axarr[1,1]
    #plot_18gep(ax)
    #ax.set_ylim(-15.5,-20.5)

    for ax in axarr:
        ax.tick_params(axis='both', labelsize=12)

    axarr[1].set_xlabel('Days after GRB', fontsize=14)
    axarr[0].set_ylabel('Abs Mag', fontsize=14)
    #plt.tight_layout()
    plt.savefig("early_lc_collage.eps", format='eps', dpi=300, bbox_inches='tight')

def bottom_panel():
    """ 'ordinary' SNe """
    fig,axarr = plt.subplots(1,3,figsize=(10,3), sharex=True, sharey=False)
    ax = axarr[0]
    plot_14gqr(ax)
    ax.set_xlim(-0.5,4)
    ax.set_ylim(-16.6,-17.5)
    
    ax = axarr[1]
    plot_2011dh(ax)
    ax.set_ylim(-14,-16)

    ax = axarr[2]
    plot_15dtg(ax)
    ax.set_ylim(-16,-17.5)

    for ax in axarr:
        ax.tick_params(axis='both', labelsize=12)

    axarr[1].set_xlabel('Days', fontsize=14)
    axarr[0].set_ylabel('Abs Mag', fontsize=14)
    plt.tight_layout()
    plt.savefig("early_lc_collage_sn.png", dpi=300, bbox_inches='tight')
    #plt.show()

if __name__=="__main__":
    top_panel()
    #bottom_panel()
