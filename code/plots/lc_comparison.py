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
rcol = 'Crimson'
gcol = 'Aquamarine'
icol = 'Goldenrod'
t0 = 2458883.17
msize = 6
data_dir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data"


def sn2020bvc():
    z = 0.025201
    dm = Planck15.distmod(z=z).value
    dat = ascii.read("../../data/marshal_lc.txt")
    dt = dat['jdobs']-t0
    filt = dat['filter']
    mag = dat['magpsf']
    emag = dat['sigmamagpsf']
    inst = dat['instrument']
    maglim = dat['limmag']-dm

    # r-band light curve
    ax = axarr[1]
    choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='r', mag < 99))
    ax.errorbar(dt[choose], mag[choose]-dm-EXT_R, emag[choose], ms=msize, fmt='o', 
            mfc=rcol, mec=rcol, label='20bvc $r$', c=rcol, zorder=1)

    # r-band upper limits
    choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='r', mag == 99, dt<5))
    ax.scatter(
            dt[choose], maglim[choose]-dm-EXT_R, marker='o', 
            c=rcol, label=None, zorder=1)

    for ii,val in enumerate(dt[choose]):
        ax.arrow(
                val, maglim[choose][ii]-dm-EXT_R, 0, 0.25, 
                color=rcol, head_length=0.2, head_width=1, zorder=1)

    # g-band light curve
    ax = axarr[0]
    choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='g', mag < 99))
    ax.errorbar(dt[choose], mag[choose]-dm-EXT_G, emag[choose], ms=msize, fmt='s', 
            mfc=gcol, mec=gcol, label='20bvc $g$', c=gcol, zorder=2)

    # g-band upper limits
    choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='g', mag == 99, dt<5))
    ax.scatter(
            dt[choose], maglim[choose]-dm-EXT_G, marker='s', c=gcol, 
            label=None, zorder=2)
    for ii,val in enumerate(dt[choose]):
        ax.arrow(
                val, maglim[choose][ii]-dm-EXT_G, 0, 0.25, 
                color=gcol, head_length=0.2, head_width=1, zorder=2)


def sn2006aj():
    dm = Planck15.distmod(z=0.0335).value
    dat = ascii.read("../../data/lc_060218_full.txt")
    b_plot = dat['band']
    dt_plot = dat['time']-53784.149
    t_plot = dat['time']
    m_plot = dat['magnitude']-dm

    # g-ish band
    ax = axarr[0]
    choose = np.logical_and(b_plot == 'B', dat['upperlimit']=='F')
    ax.plot(
        dt_plot[choose]/1.033, m_plot[choose]-0.527, linestyle='--', lw=1, 
        color='k', label="06aj $B$", zorder=0, alpha=1)

    # r-ish band
    ax = axarr[1]
    choose = np.logical_and(b_plot == 'V', dat['upperlimit']=='F')
    ax.plot(dt_plot[choose]/1.033, m_plot[choose]-0.399, linestyle='--', 
        lw=1, color='k', label="06aj $V$", zorder=0)


def sn1998bw():
    dm = Planck15.distmod(z=0.0085).value
    dat = ascii.read(data_dir + "/sn1998bw.dat", delimiter=';')

    # r-band
    ax = axarr[1]
    jd = dat['JD']
    rband = dat['Rcmag']
    erband = dat['e_Rcmag']
    # Extinction is 0.127 in R-band in this direction
    ax.plot(jd-jd[0], rband-dm-0.127, color='Goldenrod', lw=1,
            label="98bw $Rc$")

    # g-band
    ax = axarr[0]
    gband = dat['Bmag']
    egband = dat['e_Bmag']
    # Extinction is 0.212 in B-band in this direction
    ax.plot(jd-jd[0], gband-dm-0.212, color='Goldenrod', lw=1, 
            label="98bw $B$")


def sn2017iuk():
    dm = Planck15.distmod(z=0.037).value
    dat = ascii.read(data_dir + "/lc_171205a.txt")
    dt = dat['col1'].data / 86400
    mag = dat['col3'].data-dm
    filt = dat['col12'].data
    ax = axarr[0]
    choose = filt == 'B'
    # correcting for MW extinction on my own
    ax.plot(
            dt[choose], mag[choose]-0.055, 
            c='#84206b', ls=':', label="17iuk $B$") 
    choose = filt == 'V'
    ax = axarr[1]
    ax.plot(
            dt[choose], mag[choose]-0.042, 
            c='#84206b', ls=':', label="17iuk $V$") 
    

def sn2010bh():
    dm = Planck15.distmod(z=0.059).value
    dat = ascii.read(data_dir + "/sn2010bh.txt")
    dt = dat['col2']-55271.53113425926
    mag = dat['col3'].data-dm
    filt = dat['col6'].data
    ax = axarr[0]
    choose = filt == 'B'
    ax.scatter(
            dt[choose], mag[choose], edgecolor='grey', facecolor='white', 
            marker='D', lw=0.5, label="10bh $B$") 
    ax = axarr[1]
    choose = filt == 'R_c'
    ax.scatter(
            dt[choose], mag[choose], edgecolor='grey', facecolor='white', 
            marker='D', lw=0.5, label="10bh $Rc$") 


if __name__=="__main__":
    # Two panels: left is g, right is r
    fig,axarr = plt.subplots(1,2, figsize=(10,4))

    # Show SNe
    sn2020bvc()
    sn2006aj()
    sn1998bw()
    sn2017iuk()
    sn2010bh()

    # Formatting, axis labeling
    axarr[0].set_ylabel(r"Absolute Mag", fontsize=14)

    for ax in axarr:
        ax.invert_yaxis()
        ax.set_xlabel("Days After First Light", fontsize=16)
        ax.tick_params(labelsize=14)
        ax.tick_params(labelsize=14)
        ax.set_xlim(-2.5,31)
        ax.set_ylim(-14,-19.5)
        ax.legend(loc='lower right', fontsize=12, ncol=3)
    plt.tight_layout()
    plt.savefig("lc_comparison.png", dpi=200)
    #plt.show()
