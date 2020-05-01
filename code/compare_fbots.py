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

    # 12brf

    # 10bjp

    # 2015K

    # 12bv

    # 04D4ec

    # 18cow


if __name__=="__main__":
    fig,ax = plt.subplots(1,1, figsize=(5,4))
    sn20bvc(ax)
    drout(ax)
    ax.invert_yaxis()
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xlabel("Time [days]", fontsize=14)
    ax.set_ylabel("Absolute Mag", fontsize=14)
    ax.set_xlim(-2,30)
    ax.legend(fontsize=12)

    plt.tight_layout()
    plt.show()
