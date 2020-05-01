""" Compare the optical light curve of SN2020bvc to FBOTs in the literature """

import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck15
from get_lc import get_opt_lc


def sn20bvc(ax):
    t,mag,emag,maglim,filt,inst = get_opt_lc()
    dm = Planck15.distmod(z=0.0252).value

    choose = np.logical_and.reduce(
            (filt == 'g', mag < 50, inst=='P48+ZTF'))
    t0 = t[choose][0]-1
    ax.errorbar(
            (t[choose]-t0)/1.0252, mag[choose]-dm, emag[choose], 
            fmt='o', c='mediumAquaMarine', label='P48 $g$')

    choose = np.logical_and.reduce(
            (filt == 'r', mag < 50, inst=='P48+ZTF'))
    ax.errorbar(
            (t[choose]-t0)/1.0252, mag[choose]-dm, emag[choose], 
            fmt='o', c='Crimson', label='P48 $r$')
    ax.scatter((t[choose][0]-t0-0.67)/1.0252, 19.4-dm, marker='v', c='Crimson')


def drout(ax):
    # Get the Drout14 light curves
    inputf = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc/drout14.txt"
    dat = np.loadtxt(inputf, delimiter=';', dtype=str)
    name = dat[:,0]
    filt = dat[:,1]
    dt = dat[:,2].astype(float)
    islim = dat[:,3] == '<'
    mag = dat[:,4].astype(float)
    emag =dat[:,5]

    # Plot 10ah
    dm = Planck15.distmod(z=0.074).value
    choose = np.logical_and.reduce((filt=='g_P1', name=='10ah ', ~islim))
    ax.plot(
            dt[choose]/1.074, mag[choose]-dm, c='mediumAquaMarine', 
            lw=1, ls='--', label='PS1-10ah')


def ksn2015k(ax):
    inputf = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc/ksn2015k.txt"
    dat = np.loadtxt(inputf)
    dt = dat[:,0]
    M = dat[:,1]
    ax.plot(
            dt/1.09, M, c='Goldenrod', lw=1, ls='--', label='KSN2015K')


def at2018cow(ax):
    inputf = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/lc/18cow.txt"
    dm = Planck15.distmod(z=0.0141).value
    dat = np.loadtxt(inputf, dtype=str)
    tmjd = dat[:,0].astype(float)
    filt = dat[:,2]
    mag = dat[:,4].astype(float)
    choose = filt == 'g'
    ax.plot(tmjd[choose]-tmjd[choose][0], mag[choose]-dm, c='mediumAquaMarine', lw=0.5, label='18cow')
    # last upper limit, which was also in g-band
    ax.plot([58284.1300-tmjd[choose][0], 0], [18.90-dm, mag[choose][0]-dm], 
            c='mediumAquaMarine', lw=0.5, label='_nolabel', ls=':')


if __name__=="__main__":
    fig,ax = plt.subplots(1,1, figsize=(5,4))
    sn20bvc(ax)
    drout(ax)
    ksn2015k(ax)
    at2018cow(ax)
    ax.invert_yaxis()
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xlabel("Time [days]", fontsize=14)
    ax.set_ylabel("Absolute Mag", fontsize=14)
    ax.set_xlim(-2,30)
    ax.legend(fontsize=11, ncol=3)

    plt.tight_layout()
    plt.savefig("sn20bvc_fbot_comparison.png", dpi=200)
    #plt.show()
