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
            emag_use.append(1/np.sqrt(sum(ivar)))
    dt_use = np.array(dt_use)
    mag_use = np.array(mag_use)
    emag_use = np.array(emag_use)
    return dt_use,mag_use,emag_use


def get_binned_p48(dt,mag,emag,filt,instr):
    # get binned P48 light curves
    rdt,rmag,remag = get_binned_p48_band(dt,mag,emag,filt,instr,'r')
    gdt,gmag,gemag = get_binned_p48_band(dt,mag,emag,filt,instr,'g')
    idt,imag,iemag = get_binned_p48_band(dt,mag,emag,filt,instr,'i')
    dt_all = np.hstack((rdt, gdt, idt))
    mag_all = np.hstack((rmag, gmag, imag))
    emag_all = np.hstack((remag, gemag, iemag))
    filt_all = np.hstack((['r']*len(rdt), ['g']*len(gdt), ['i']*len(idt)))
    return dt_all,mag_all,emag_all,filt_all


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
dtbins = [0.76, 1.36, 1.8, 2.8, 3.8, 4.74, 5.78]
bs = 0.15 # bin size

# for each bin, select photometry within 0.1d

for ii,dtbin in enumerate(dtbins):
    # choose the panel
    ax = axarr.flatten()[ii]
    xvals = []
    yvals = []
    eyvals = []

    #Get UVOT photometry
    choose = np.abs(uvdt-dtbin)<bs
    for jj in np.arange(sum(choose)):
        wl = bands[uvfilt[choose][jj]]
        f = uvflux[choose][jj]*1E3
        ef = uveflux[choose][jj]*1E3
        xvals.append(wl)
        yvals.append(f)
        eyvals.append(ef)
    #LT
    choose = np.logical_and(np.abs(dt-dtbin)<bs, instr=='LT+SPRAT')
    for jj in np.arange(sum(choose)):
        wl = bands[filt[choose][jj]]
        f,ef = toflux(mag[choose][jj],emag[choose][jj])
        xvals.append(wl)
        yvals.append(f)
        eyvals.append(ef)
    #P48
    dt_p48,mag_p48,emag_p48,filt_p48 = get_binned_p48(dt,mag,emag,filt,instr)
    choose = np.abs(dt_p48-dtbin) < 0.5
    for jj in np.arange(sum(choose)):
        wl = bands[filt_p48[choose][jj]]
        f,ef = toflux(mag_p48[choose][jj],emag_p48[choose][jj])
        xvals.append(wl)
        yvals.append(f)
        eyvals.append(ef)

    # Done loading photometry; now plot and measure bb
    xvals = np.array(xvals)
    yvals = np.array(yvals)
    eyvals = np.array(eyvals)

    # Sort in order of wavelength
    order = np.argsort(xvals)
    xvals = xvals[order]
    yvals = yvals[order]
    eyvals = eyvals[order]

    ax.errorbar(xvals, yvals, yerr=eyvals, c='k', fmt='.')
    txt = "$\Delta t \\approx $" + str(np.round(dtbin,1))
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
                bb_func, xvals*1E-8, ysamples[jj], p0=[10000,1E14],
                bounds=([2000,1E12],[50000, 1E17])) # must be positive
        temps[jj] = popt[0]
        radii[jj] = popt[1]
        xplot = np.linspace(2000,8000)
        yplot = bb_func(xplot*1E-8, popt[0], popt[1])
        ax.plot(xplot,yplot,lw=0.1,alpha=0.1)
    lums = 4*np.pi*radii**2 * (5.67E-5)*temps**4

    # Sort results and calculate 16-to-84 percentile range
    print(dtbin)

    nvals = len(temps)
    start_ind = int(0.16*nvals)
    end_ind = int(0.84*nvals)
    temps = np.sort(temps)
    T = np.median(temps)
    low = T-temps[start_ind]
    hi = temps[end_ind]-T
    print("%s +%s -%s" %(T/1E3, hi/1E3, low/1E3))
    radii = np.sort(radii)
    R = np.median(radii)
    low = R-radii[start_ind]
    hi = radii[end_ind]-R
    print("%s +%s -%s" %(R/1E14, hi/1E14, low/1E14))
    lums = np.sort(lums)
    L = np.median(lums)
    low = L-lums[start_ind]
    hi = lums[end_ind]-L
    print("%s +%s -%s" %(L/1E42, hi/1E42, low/1E42))


    # Calculate the chi squared of the final fit
    chisq = sum((yvals-bb_func(xvals*1E-8, T, R))**2/eyvals**2)
    dof = len(yvals)-2
    print(chisq, dof)
    print(chisq/dof)

    xplot = np.linspace(2000,8000)
    yplot = bb_func(xplot*1E-8, T, R)
    ax.plot(xplot,yplot,lw=0.5,alpha=1,c='Crimson')

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
