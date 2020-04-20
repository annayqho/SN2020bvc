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


def sn17iuk():
    # Now plot the earliest spectra of SN2017iuk, from Izzo+2019
    datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"

    comp = fits.open(
            datadir + "/GTC_OSIRIS_GRB171205A/1D_GRB171205A_171206_R1000B.fits")[0].data
    wl = np.arange(3629.598, 3629.598+len(comp)*2.071432195122, 2.071432195122)
    plt.plot(
            wl/1.0368, comp/2, c='#84206b', 
            drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
    plt.text(wl[-1]/1.03, comp[-1]/2, '17iuk/2 (0.95d)', fontsize=12)


def floyds():
    # Plot the early FLOYDS spectrum of 20bvc
    datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"

    dat = np.loadtxt(datadir + "/floyds_feb5.txt")
    wl = dat[:,0]
    comp = dat[:,1]
    plt.plot(
            wl/1.0252, comp*6, c='grey', 
            drawstyle='steps-mid', lw=0.5, ls='-', zorder=0)
    plt.text(7450, 7E-16, '20bvc (x6; FLOYDS 2d)', fontsize=12)
    # line IDs: Ca II
    v = 70000
    z = 0.0252
    caii = np.array([8498,8542,8662])*(1+z)/(1+v/3E5)
    plt.scatter(caii, np.array([9E-16]*3), marker='|', c='k')
    plt.text(caii[0]/1.05, 1E-15, 'CaII (70,000 km/s)', fontsize=12,
            horizontalalignment='left', verticalalignment='bottom')
    feii = np.array([4924,5018,5169])*(1+z)/(1+v/3E5)
    plt.scatter(feii, np.array([2.3E-15]*3), marker='|', c='k')
    plt.text(feii[-1]/1.05, 2.3E-15, 'FeII', fontsize=12,
            horizontalalignment='left', verticalalignment='bottom')



def sn06aj():
    datadir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"
    comp = ascii.read(datadir + "/SN2006aj_2006-02-20_00-00-00_BTA-6_SCORPIO_None.dat")
    wl = comp['col1']
    f_lam = comp['col2'] # erg/cm2/s/AA, I think
    plt.plot(
            wl/1.0335, f_lam, c='#f6c746',
            drawstyle='steps-mid', lw=0.5, ls='-', zorder=0, label="_none")
    plt.text(wl[-1]/1.0335, f_lam[-1], '06aj (2.6d)', fontsize=12,
            verticalalignment='top')


if __name__=="__main__":
    fig,ax = plt.subplots(1,1,figsize=(7,5))

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
    plt.plot(
            wl/1.0252, f_lam, c='k', drawstyle='steps-mid', 
            lw=0.5, ls='-', alpha=1)

    # Fit for bb=20,000K
    xplt, yplt = plot(xvals, yvals, eyvals, 20000) # highest reasonable temperature
    plt.plot(xplt,yplt,lw=0.5,alpha=1,c='k', ls='--', label="$T=20,000$K")

    # Fit for bb=13,200
    xplt, yplt = plot(xvals, yvals, eyvals, 13200) # highest reasonable temperature
    plt.plot(xplt,yplt,lw=1,alpha=1,c='k', ls='-', label="$T=13,200$K")

    T = '13.2'
    R = '5.1'
    L = '5.6'

    plt.text(
        3378, 1E-15, "20bvc (0.7d)", fontsize=14, 
        horizontalalignment='left', verticalalignment='top')

    # Comparison objects
    sn17iuk()
    sn06aj()
    floyds()

    plt.ylabel(r"Flux (erg\,s$^{-1}$\,cm$^{-2}\,$\AA$^{-1}$)", fontsize=16)
    plt.xlabel("Rest Wavelength (\AA)", fontsize=16)
    plt.yscale('log')
    plt.ylim(4E-17, 3E-15)
    plt.xlim(3300, 9300)
    plt.tick_params(axis='both', labelsize=16)
    plt.legend(fontsize=12)
    plt.tight_layout()  
    plt.show()

    #plt.savefig("bb_first_epoch.png", dpi=200)
