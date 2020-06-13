""" fit blackbodies to the photometry from the marshal """

import numpy as np
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.io import ascii
from scipy.optimize import curve_fit
from astropy.cosmology import Planck15
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/SN2020bvc/code")
from get_lc import get_opt_lc, get_uv_lc

t0 = 2458883.17

# a mapping from filter to central wavelength
bands = {}
bands['V'] = 5468
bands['B'] = 4392
bands['U'] = 3465
bands['UVW1'] = 2600
bands['UVM2'] = 2246
bands['UVW2'] = 1928
# p48
bands['g'] = 4722.7 
bands['r'] = 6339.6
bands['i'] = 7886.1
# lt: assume everything is the same except
bands['u'] = 3513.7
bands['z'] = 8972.9


def toflux(mag,emag):
    """ Convert from magnitude to flux in uJy """
    f = 1E6 * 10**((mag-8.90)/(-2.5))
    ef = f*emag
    return f,ef


def bb_func(wl,T,R):
    """
    Return a blackbody function

    Parameters
    ----------
    wl: wavelength in angstroms
    T: temperature in Kelvin
    R: radius in cm
    """
    d = Planck15.luminosity_distance(z=0.02507).cgs.value
    h = 6.626E-27
    c = 3E10
    k = 1.38E-16
    Blam = (2*h*c**2/wl**5) * (1/(np.exp(h*c/(wl*k*T)-1)))
    flam = Blam * np.pi * R**2 / d**2
    # in units of uJy
    fnu = wl**2 * flam / c
    return fnu / 1E-23 / 1E-6


def get_binned_p48_band(dt,mag,emag,filt,instr,use_filt='r'):
    p48 = instr == 'P48+ZTF'
    band = filt==use_filt
    bins = np.unique(
            np.array([np.round(val,0) for val in dt[np.logical_and(p48, band)]]))
    dt_use = []
    mag_use = []
    emag_use = []

    for ii,b in enumerate(bins):
        use = np.logical_and.reduce((p48, band, mag<50, np.abs(dt-b)<0.5))
        if sum(use) == 1:
            dt_use.append(dt[use][0])
            mag_use.append(mag[use][0])
            emag_use.append(emag[use][0])
        elif sum(use) > 1:
            dt_use.append(np.mean(dt[use]))
            ivar = 1/(emag[use])**2
            mag_use.append(sum(mag[use]*ivar)/sum(ivar))
            emag_use.append(1/sqrt(sum(ivar)))
    dt_use = np.array(dt_use)
    mag_use = np.array(mag_use)
    emag_use = np.array(emag_use)
    return dt_use,mag_use,emag_use


def get_binned_p48(dt,mag,emag,filt,instr):
    # get binned P48 light curves
    rdt,rmag,remag = get_binned_p48_band(dt,mag,emag,filt,instr,'r')
    gdt,gmag,gemag = get_binned_p48_band(dt,mag,emag,filt,instr,'g')
    idt,imag,iemag = get_binned_p48_band(dt,mag,emag,filt,instr,'i')
    return rdt,rmag,remag,gdt,gmag,gemag,idt,imag,iemag


# Ultimately, here are the parameters I will measure
Lbol = []
eLbol = []
Teff = []
eTeff = []
Rph = []
eRph = []

t,mag,emag,maglim,filt,instr = get_opt_lc()
uvdt,uvfilt,uvflux,uveflux = get_uv_lc()
dt = t-t0

fig,axarr = plt.subplots(3, 5, figsize=(8,5), sharex=True, sharey=True)

# Define time bins
# For dt < 11d, use Swift+UVOT epochs 
# For dt > 11d, use epochs of LT photometry
# In the paper, we only present photometry up to 30 days
early_bins = np.logical_and(instr=='Swift+UVOT', dt<11)
late_bins = np.logical_and(instr=='LT+IOO', dt>11)
use = np.logical_or(early_bins, late_bins)
# and round to one decimal place to group observations
dtbins = np.unique(np.array([np.round(val,1) for val in dt[use]]))

# for each bin, plot the UVOT and/or LT photometry
# and interpolate the p48 light curve onto that epoch
# also round the dts to one decimal place
uvdt = np.array([np.round(val,1) for val in uvdt])
dt = np.array([np.round(val,1) for val in dt])

for ii,dtbin in enumerate(dtbins):
    # choose the panel
    ax = axarr.flatten()[ii]
    xvals = []
    yvals = []
    eyvals = []

    if dtbin < 15:
        #Get UVOT photometry
        choose = uvdt == dtbin
        for jj in np.arange(sum(choose)):
            wl = bands[uvfilt[choose][jj]]
            f = uvflux[choose][jj]*1E3
            ef = uveflux[choose][jj]*1E3
            xvals.append(wl)
            yvals.append(f)
            eyvals.append(ef)
    else:
        #LT
        choose = np.logical_and(dt == dtbin, instr=='LT+SPRAT')
        for jj in np.arange(sum(choose)):
            wl = bands[filt[choose][jj]]
            f,ef = toflux(mag[choose][jj],emag[choose][jj])
            xvals.append(wl)
            yvals.append(f)
            eyvals.append(ef)

    #For all dtbins, interpolate P48 photometry
    p48_rdt,p48_rmag,p48_remag,p48_gdt,p48_gmag,p48_gemag,p48_idt,p48_imag,p48_iemag = \
            get_binned_p48(dt,mag,emag,filt,instr)


        ztfmag = np.interp(dtbin, dt[choose], mag[choose])
        # THIS ERROR BAR IS TEMPORARY
        f,ef = toflux(ztfmag,emag=0.1) 
        xvals.append(bands[ztffilt])
        yvals.append(f)
        eyvals.append(ef)

    xvals = np.array(xvals)
    yvals = np.array(yvals)
    eyvals = np.array(eyvals)

    ax.errorbar(xvals, yvals, yerr=eyvals, c='k', fmt='.')
    txt = "$\Delta$t=" + str(np.round(dtbin,2))
    ax.text(0.9, 0.9, txt, transform=ax.transAxes,
            horizontalalignment='right', verticalalignment='top')

    # Fit a blackbody 1000 times
    nsim = 10
    temps = np.zeros(nsim)
    radii = np.zeros(nsim)
    ysamples = np.zeros((nsim, len(xvals)))
    for jj,val in enumerate(yvals):
        ysamples[:,jj] = np.random.normal(loc=val,scale=eyvals[jj],size=nsim)
    for jj in np.arange(nsim):
        popt, pcov = curve_fit(
                bb_func, xvals*1E-8, ysamples[jj], p0=[5000,1E14])
        temps[jj] = popt[0]
        radii[jj] = popt[1]
        xplot = np.linspace(2000,8000)
        yplot = bb_func(xplot*1E-8, popt[0], popt[1])
        ax.plot(xplot,yplot,lw=0.1,alpha=0.1)
    lums = 4*np.pi*radii**2 * (5.67E-5)*temps**4
    T = np.mean(temps)
    Teff.append(T)
    eT = np.std(temps)
    eTeff.append(eT)
    R = np.mean(radii)
    Rph.append(R)
    eR = np.std(radii)
    eRph.append(eR)
    L = np.mean(lums)
    Lbol.append(L)
    eL = np.std(lums)
    eLbol.append(eL)

    print(dtbin)
    print("%s +/- %s" %(T/1E3, eT/1E3))
    print("%s +/- %s" %(R/1E14, eR/1E14))
    print("%s +/- %s" %(L/1E42, eL/1E42))

    # Calculate the chi squared of the final fit
    chisq = sum((yvals-bb_func(xvals*1E-8, T, R))**2/eyvals**2)
    dof = len(yvals)-2
    print(chisq, dof)

    # Plot the blackbody
    #xplot = np.linspace(2000,8000)
    #yplot = bb_func(3E10/(xplot*1E-8), T, R)
    #ax.plot(xplot,yplot,lw=0.5,alpha=1,c='k')

Lbol = np.array(Lbol)
eLbol = np.array(eLbol)
Rph = np.array(Rph)
eRph = np.array(eRph)
Teff = np.array(Teff)
eTeff = np.array(eTeff)

axarr[0,0].set_xlim(1E3, 2E4)
axarr[0,0].set_ylim(10, 5000)
axarr[0,0].set_yscale('log')
axarr[0,0].set_xscale('log')
plt.subplots_adjust(wspace=0,hspace=0)
#fig.text(0.5,0.04,"Wavelength (AA)", ha='center',fontsize=14)
fig.text(0.04,0.5,r'Flux ($\mu$Jy)',fontsize=14,verticalalignment='center',
        horizontalalignment='center',rotation='vertical')
axarr[1,2].set_xlabel(r'Wavelength (\AA)',fontsize=14)
#plt.tight_layout()

plt.show()
#plt.savefig("bbfits.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
