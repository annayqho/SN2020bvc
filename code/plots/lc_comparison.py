""" Compare the g and r-band light curves to GRB-SNe """

import numpy as np
from astropy.io import ascii
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.cosmology import Planck15

EXT_G = 0.041
EXT_R = 0.029
EXT_I = 0.021
rcol = '#e55c30'
gcol = '#140b34'
icol = 'k'
z = 0.025201
t0 = 2458883.17
msize = 6
data_dir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data"


def sn2020bvc():
    dat = ascii.read("../../data/marshal_lc.txt")
    dt = dat['jdobs']-t0
    filt = dat['filter']
    mag = dat['magpsf']
    emag = dat['sigmamagpsf']
    inst = dat['instrument']
    maglim = dat['limmag']

    # r-band light curve
    ax = axarr[1]
    choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='r', mag < 99))
    ax.errorbar(dt[choose], mag[choose]-EXT_R, emag[choose], ms=msize, fmt='o', 
            mfc=rcol, mec=rcol, label='20bvc $r$', c=rcol, zorder=1)

    # r-band upper limits
    choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='r', mag == 99, dt<5))
    ax.scatter(
            dt[choose], maglim[choose]-EXT_R, marker='o', 
            c=rcol, label=None, zorder=1)

    for ii,val in enumerate(dt[choose]):
        ax.arrow(
                val, maglim[choose][ii]-EXT_R, 0, 0.25, 
                color=rcol, head_length=0.2, head_width=1, zorder=1)

    # g-band light curve
    ax = axarr[0]
    choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='g', mag < 99))
    ax.errorbar(dt[choose], mag[choose]-EXT_G, emag[choose], ms=msize, fmt='s', 
            mfc=gcol, mec=gcol, label='20bvc $g$', c=gcol, zorder=2)

    # g-band upper limits
    choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='g', mag == 99, dt<5))
    ax.scatter(
            dt[choose], maglim[choose]-EXT_G, marker='s', c=gcol, 
            label=None, zorder=2)
    for ii,val in enumerate(dt[choose]):
        ax.arrow(
                val, maglim[choose][ii]-EXT_G, 0, 0.25, 
                color=gcol, head_length=0.2, head_width=1, zorder=2)


def sn2006aj():
    dm = Planck15.distmod(z=0.0335).value-Planck15.distmod(z=0.025235).value
    print(dm)
    dat = ascii.read("../../data/lc_060218_full.txt")
    t = dat['time']
    mag = dat['magnitude']
    emag = dat['e_magnitude']
    band = dat['band']

    # g-ish band
    ax = axarr[0]
    dt_plot = np.array([2.7, 2.9, 4.0, 4.9, 6.0, 6.9, 7.9, 8.9, 9.9, 10.9, 
        11.9, 12.9, 13.9, 14.9, 15.9, 17.7, 17.9, 18.8, 19.9, 21.8, 22.9, 23.8])
    m_plot = np.array([18.14, 18.07, 17.87, 17.68, 17.74, 17.56, 17.54, 17.47, 17.43, 17.49, 
        17.57, 17.66, 17.81, 17.95, 18.11, 18.56, 18.64, 18.73, 18.95, 19.16, 19.43, 19.43])
    ax.plot(dt_plot, m_plot-dm, linestyle='--', lw=1, color='k', label="06aj $B$", zorder=0, alpha=1)

    dat = ascii.read("../../data/lc_060218_swift.dat")
    dt_plot = dat['t']
    m_plot = dat['mag']
    ax.plot(dt_plot, m_plot-dm-0.5, linestyle='-', lw=1, color='k', zorder=0, alpha=1, label="_none")

    # r-ish band
    ax = axarr[1]
    dt_plot = np.array([2.7, 2.9, 4.0, 4.9, 5.0, 6.0, 6.7, 6.9, 7.9, 8.7, 
       8.9, 8.9, 9.9, 10.9, 11.9, 12.7, 12.9, 13.9, 14.9, 15.9, 
       16.7, 17.9, 17.9, 18.9, 18.9, 19.9, 19.9, 21.8, 22.9, 23.9, 
       24.9, 25.9])
    #These V-band magnitudes are corrected for host-galaxy flux, 
    #extinction, host-galaxy extinction
    m_plot = np.array(
           [17.86, 17.83, 17.54, 17.41, 17.37, 17.28, 17.19, 17.16, 17.08, 17.05, 
           17.02, 17.03, 17.02, 17.02, 17.04, 17.07, 17.08, 17.14, 17.18, 17.27, 
           17.26, 17.46, 17.47, 17.54, 17.54, 17.65, 17.61, 17.80, 17.89, 17.99, 
           18.12, 18.14])
    ax.plot(dt_plot, m_plot-dm, linestyle='--', lw=1, color='k', label="06aj $V$", zorder=0)


def sn1998bw():
    dm = Planck15.distmod(z=0.0085).value-Planck15.distmod(z=0.025235).value
    print(dm)   
    dat = ascii.read(data_dir + "/sn1998bw.dat", delimiter=';')

    # r-band
    ax = axarr[1]
    jd = dat['JD']
    rband = dat['Rcmag']
    erband = dat['e_Rcmag']
    # Extinction is 0.127 in R-band in this direction
    ax.plot(jd-jd[0], rband-dm-0.127, color='#e55c30', lw=1,
            label="98bw $Rc$")

    # g-band
    ax = axarr[0]
    gband = dat['Bmag']
    egband = dat['e_Bmag']
    # Extinction is 0.212 in B-band in this direction
    ax.plot(jd-jd[0], gband-dm-0.212, color='#e55c30', lw=1, 
            label="98bw $B$")


def sn2017iuk():
    dm = Planck15.distmod(z=0.037).value-Planck15.distmod(z=0.025235).value
    dat = ascii.read(data_dir + "/lc_171205a.txt")
    dt = dat['col1'].data / 86400
    mag = dat['col3'].data-dm
    filt = dat['col12'].data
    ax = axarr[0]
    choose = filt == 'B'
    ax.plot(
            dt[choose], mag[choose], c='#84206b', ls=':', label="17iuk $B$") 
    choose = filt == 'V'
    ax = axarr[1]
    ax.plot(
            dt[choose], mag[choose], c='#84206b', ls=':', label="17iuk $V$") 
    


if __name__=="__main__":
    # Two panels: left is g, right is r
    fig,axarr = plt.subplots(1,2, figsize=(10,4))

    # Show SNe
    sn2020bvc()
    sn2006aj()
    sn1998bw()
    sn2017iuk()

    # Formatting, axis labeling
    axarr[0].set_ylabel(r"Apparent Mag ($z=0.025201$)", fontsize=14)

    for ax in axarr:
        ax2 = ax.twinx()
        y_f = lambda y_i: y_i - Planck15.distmod(z=z).value
        ymin, ymax = ax.get_ylim()
        ax2.set_ylim((y_f(ymin), y_f(ymax)))
        ax2.plot([],[])
        ax2.tick_params(axis='both', labelsize=14)
        ax2.invert_yaxis()

    ax2.set_ylabel(
            "Absolute Mag ($z=0.025201$)",
            fontsize=14, rotation=270, labelpad=15.0)

    for ax in axarr:
        ax.invert_yaxis()
        ax.set_xlabel("Days After First Light", fontsize=16)
        ax.tick_params(labelsize=14)
        ax.tick_params(labelsize=14)
        ax.set_xlim(-2.5,31)
        ax.set_ylim(21.2,15.7)
        ax.legend(loc='lower right', fontsize=14, ncol=2)
    plt.tight_layout()
    plt.savefig("lc_comparison.png", dpi=200)
    #plt.show()
