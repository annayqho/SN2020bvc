""" Fit blackbody SED but also spectrum on the first day """

import numpy as np
from astropy.io import ascii
from astropy.io import fits
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.cosmology import Planck15
from astropy.io import fits as pyfits
from scipy.optimize import curve_fit

def plot_day4(ax):
    z = 0.0252
    datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"
    spec = np.loadtxt(datadir + "/ZTF20aalxlis_20200207_P60_v1.ascii")
    wl = spec[:,0]
    f_lam = spec[:,1] # erg/cm2/s/AA, I think
    ef_lam = np.sqrt(spec[:,2])
    hz = 3E10/(wl*1E-8)
    f_nu = wl * f_lam / hz
    ef_nu = wl * ef_lam / hz

    ax.plot(
            wl/1.0252, f_lam/1E-15, c='#e55c30', drawstyle='steps-mid', 
            lw=1.5, ls='-', alpha=1)

    # Indicate Ca II at v=60,000 km/s
    v = 60000
    caii = np.array([8498,8542,8662])*(1+z)*np.sqrt((1-v/3E5)/(1+v/3E5))
    ax.scatter(caii[1], 0.3, marker='|', c='k')
    ax.text(caii[0]/1.05, 0.32, 'CaII (60,000 km/s)', fontsize=12,
            horizontalalignment='left', verticalalignment='bottom')

    # Add a spectrum of 17iuk for comparison
    z = 0.0368
    dat = np.loadtxt(
            datadir + "/GTC_OSIRIS_GRB171205A/1D_GRB171205A_171207_R1000B_clipped.ascii")
    comp = dat[:,1]
    wl = dat[:,0]
    ax.plot(
            wl/1.0368, comp*3/1E-15, c='#84206b',
            drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
    ax.text(wl[-1]/1.03, comp[-1]*3/1E-15, '17iuk (2d)', fontsize=12)

    # Indicate Si II
    v = 75000
    si = 6355*(1+z)*np.sqrt((1-v/3E5)/(1+v/3E5))
    ax.scatter(si, 21E-2, marker='|', c='k')
    ax.text(si/1.05, 21E-2, 'SiII', fontsize=12,
            horizontalalignment='left', verticalalignment='bottom')

    # Indicate Ca II at v=65,000 km/s
    v = 105000
    caii = np.array([8498,8542,8662])*(1+z)*np.sqrt((1-v/3E5)/(1+v/3E5))
    ax.scatter(caii[1], 9E-2, marker='|', c='k')
    ax.text(caii[0], 9E-2, 'CaII', fontsize=12,
            horizontalalignment='left', verticalalignment='bottom')

    # Add a 4d spectrum of 06aj for comparison
    comp = ascii.read(
            datadir + "/SN2006aj_2006-02-21_00-00-00_BTA-6_SCORPIO_None_clipped.dat")
    wl = comp['col1']
    f_lam = comp['col2'] # erg/cm2/s/AA, I think
    ax.plot(
            wl/1.0335, f_lam/1E-15, c='#f6c746',
            drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
    ax.text(wl[-1]/1.0335, f_lam[-1]/1E-15, '06aj (3.6d)', fontsize=12,
            verticalalignment='top')

    #plt.show()
