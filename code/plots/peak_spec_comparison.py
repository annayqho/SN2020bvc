""" Compare the peak-light spectrum with other GRB-SNe """

import numpy as np
from astropy.io import ascii
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.cosmology import Planck15
from astropy.io import fits as pyfits
from scipy.optimize import curve_fit

def plot_peak(ax):
    datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"
    spec = np.loadtxt(datadir + "/ZTF20aalxlis_20200216_NOT_v1_clipped.ascii")
    wl = spec[:,0]
    f_lam = spec[:,1] # erg/cm2/s/AA, I think
    hz = 3E10/(wl*1E-8)
    f_nu = wl * f_lam / hz

    ax.plot(
            wl/1.0252, f_lam/1E-15, c='#e55c30', 
            drawstyle='steps-mid', lw=1.5, ls='-', alpha=1, zorder=1)


    # Indicate Ca II 
    v = 30000
    caii = np.array([8498,8542,8662])*(1-v/3E5)
    ax.scatter(caii, [0.23]*3, marker='|', c='k')
    ax.text(caii[0], 0.25, 'CaII', fontsize=12,
            horizontalalignment='center', verticalalignment='bottom')

    # Indicate Fe II
    v = 25000
    feii = 5169*(1-v/3E5)
    ax.scatter(feii, 0.4, marker='|', c='k')
    ax.text(feii, 0.42, 'FeII', fontsize=12,
            horizontalalignment='center', verticalalignment='bottom')

    # Indicate O I
    # oi = 8446*(1.0252-v/3E5)
    # plt.scatter(oi, 0.4, marker='|', c='k')
    # plt.text(oi, 0.42, 'OI', fontsize=12,
    #         horizontalalignment='center', verticalalignment='bottom')

    # Indicate Si II 
    v = 20000
    si = 6355*(1-v/3E5)
    ax.scatter(si, 0.4, marker='|', c='k')
    ax.text(si, 0.41, 'SiII', fontsize=12,
            horizontalalignment='center', verticalalignment='bottom')

    # Plot the peak-light spectrum of 171205A for comparison
    datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"

    comp = pyfits.open(
            datadir + "/GTC_OSIRIS_GRB171205A/1D_GRB171205A_171218_R1000BR2500I.fits")[0].data
    wl = np.arange(3629.598, 3629.598+len(comp)*2.071432195122, 2.071432195122)
    ax.plot(
            wl/1.0368, 0.4 + comp*5/1E-15, c='#84206b',
            drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
    ax.text(wl[-1]/1.03, 0.45 + comp[-1]*5/1E-15, '17iuk (13d)', fontsize=12,
            verticalalignment='bottom', horizontalalignment='center')

    # Plot the peak-light spectrum of SN2006aj for comparison
    comp = ascii.read(datadir + "/sn2006aj-20060303-mmt_clipped.flm")
    wl = comp['col1']
    f_lam = comp['col2'] # erg/cm2/s/AA, I think
    ax.plot(
            wl/1.0335, f_lam-0.5, c='#f6c746', 
            drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
    ax.text(wl[-1]/1.0335, f_lam[-1]-0.5, '06aj at 12.85d', fontsize=12)
