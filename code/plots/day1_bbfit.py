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


def bb_func(nu,R,T):
    d = Planck15.luminosity_distance(z=0.02507).cgs.value
    h = 6.63E-27
    c = 3E10
    k = 1.38E-16
    Inu = (2*h*nu**3/c**2) * (1/(np.exp(h*nu/(k*T)-1)))
    fnu = Inu * np.pi * R**2 / d**2
    return fnu 


def plot(xvals, yvals, eyvals, T):
    nsim = 100
    temps = np.zeros(nsim)
    radii = np.zeros(nsim)
    ysamples = np.zeros((nsim, len(xvals)))
    for jj,val in enumerate(yvals):
        ysamples[:,jj] = np.random.normal(loc=val,scale=eyvals[jj],size=nsim)
    for jj in np.arange(nsim):
        popt, pcov = curve_fit(
                lambda nu, R: bb_func(nu, R, T), xvals, ysamples[jj], p0=[1E14])
        radii[jj] = popt[0]
    choose = radii > 0
    lums = 4*np.pi*radii[choose]**2 * (5.67E-5)*T**4
    R = np.mean(radii)
    eR = np.std(radii)
    L = np.mean(lums)
    eL = np.std(lums)

    print("%s +/- %s" %(R/1E14, eR/1E14))
    print("%s +/- %s" %(L/1E42, eL/1E42))

    xplt = np.linspace(3500,9500)
    xplt_nu = 3E10/(xplt*1E-8)
    yplt_nu = bb_func(xplt_nu, R, T)
    yplt = xplt_nu * yplt_nu / xplt

    return xplt,yplt


def sn17iuk(ax):
    # Now plot the earliest spectra of SN2017iuk, from Izzo+2019
    datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"

    dat = np.loadtxt(
            datadir + "/GTC_OSIRIS_GRB171205A/1D_GRB171205A_171206_R1000B_clipped.ascii")
    comp = dat[:,1]
    wl = dat[:,0]
    ax.plot(
            wl/1.0368, comp/1.1/1E-15, c='#84206b', 
            drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
    ax.text(4000, 1.5E-1, '17iuk/1.1 (0.95d)', fontsize=12)


def floyds(ax):
    # Plot the early FLOYDS spectrum of 20bvc
    datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"

    dat = np.loadtxt(datadir + "/floyds_feb5.txt")
    wl = dat[:,0]
    comp = dat[:,1]
    ax.plot(
            wl/1.0252, comp*6/1E-15, c='grey', 
            drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
    ax.text(7450, 7E-1, '20bvc (x6; FLOYDS 2d)', fontsize=12)
    # line IDs: Ca II
    v = 70000
    z = 0.0252
    caii = np.array([8498,8542,8662])*(1+z)*np.sqrt((1-v/3E5)/(1+v/3E5))
    ax.scatter(caii, np.array([9E-1]*3), marker='|', c='k')
    ax.text(caii[0]/1.05, 1, 'CaII (70,000 km/s)', fontsize=12,
            horizontalalignment='left', verticalalignment='bottom')
    feii = np.array([4924,5018,5169])*(1+z)*np.sqrt((1-v/3E5)/(1+v/3E5))
    ax.scatter(feii, np.array([2.3]*3), marker='|', c='k')
    ax.text(feii[-1]/1.05, 2.3, 'FeII', fontsize=12,
            horizontalalignment='left', verticalalignment='bottom')



def sn06aj(ax):
    datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"
    comp = ascii.read(
            datadir + "/SN2006aj_2006-02-20_00-00-00_BTA-6_SCORPIO_None_clipped.dat")
    wl = comp['col1']
    f_lam = comp['col2'] # erg/cm2/s/AA, I think
    ax.plot(
            wl/1.0335, f_lam/1E-15, c='#f6c746',
            drawstyle='steps-mid', lw=0.5, ls='-', zorder=0, label="_none")
    ax.text(4200, 4E-1, '06aj (2.6d)', fontsize=12, verticalalignment='bottom')


def panel(ax):
    datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"
    spec = np.loadtxt(datadir + "/ZTF20aalxlis_20200204_P60_v1.ascii")
    wl = spec[:,0] 
    f_lam = spec[:,1] # erg/cm2/s/AA, I think
    ef_lam = np.sqrt(spec[:,2])
    hz = 3E10/(wl*1E-8)
    f_nu = wl * f_lam / hz
    ef_nu = wl * ef_lam / hz

    # with the Monte Carlo stuff
    xvals = np.array(hz)
    yvals = np.array(f_nu)
    eyvals = np.array(ef_nu)

    # Plot the spectrum
    ax.plot(
            wl/1.0252, f_lam/1E-15, c='#e55c30', drawstyle='steps-mid', 
            lw=1.5, ls='-', alpha=1)

    # Fit for bb=20,000K
    xplt, yplt = plot(xvals, yvals, eyvals, 20000) # highest reasonable temperature
    ax.plot(xplt,yplt/1E-15,lw=0.5,alpha=1,c='k', ls='--', label="$T=20,000$K")

    # Fit for bb=13,200
    xplt, yplt = plot(xvals, yvals, eyvals, 13200) # highest reasonable temperature
    ax.plot(xplt,yplt/1E-15,lw=1,alpha=1,c='k', ls='-', label="$T=13,200$K")

    T = '13.2'
    R = '5.1'
    L = '5.6'

    ax.text(
        3378, 1E-15, "20bvc (0.7d)", fontsize=14, 
        horizontalalignment='left', verticalalignment='top')

    # Comparison objects
    sn17iuk(ax)
    sn06aj(ax)
    floyds(ax)

    #plt.savefig("bb_first_epoch.png", dpi=200)


def run():
    fig,ax = plt.subplots(1,1,figsize=(7,5))
    ax.set_ylabel(r"Flux (erg\,s$^{-1}$\,cm$^{-2}\,$\AA$^{-1}$)", fontsize=16)
    ax.set_xlabel("Rest Wavelength (\AA)", fontsize=16)
    ax.set_yscale('log')
    ax.set_ylim(4E-17, 3E-15)
    ax.set_xlim(3300, 9300)
    ax.tick_params(axis='both', labelsize=16)
    ax.legend(fontsize=12)
    plt.tight_layout()  
    plt.show()
