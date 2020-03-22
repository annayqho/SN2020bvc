""" fit blackbodies to the photometry from the marshal """

import numpy as np
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.io import ascii
import extinction
from scipy.optimize import curve_fit
from astropy.cosmology import Planck15

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

# extinction for each filter
ext = {}
for band in bands.keys():
    ext[band] = extinction.fitzpatrick99(
        np.array([bands[band]]), 0.034, 3.1)[0]

def toflux(mag,emag):
    f = 1E6 * 10**((mag-8.90)/(-2.5))
    ef = f*emag
    return f,ef


def bb_func(nu,T,R):
    d = Planck15.luminosity_distance(z=0.02507).cgs.value
    h = 6.63E-27
    c = 3E10
    k = 1.38E-16
    Inu = (2*h*nu**3/c**2) * (1/(np.exp(h*nu/(k*T)-1)))
    fnu = Inu * np.pi * R**2 / d**2
    # in units of uJy
    return fnu / 1E-23 / 1E-6


dat = ascii.read("../../data/marshal_lc.txt")
uvdat = ascii.read("../../data/UVOT_hostsub.ascii")
uvt = uvdat['MJD']+2400000.5
uvdt = uvt-t0
uvfilt = uvdat['FILTER']
uvflux = uvdat['AB_FNU_mJy']
uveflux = uvdat['AB_FNU_mJy_ERRM']

fig,axarr = plt.subplots(2, 5, figsize=(8,3), sharex=True, sharey=True)

t = dat['jdobs']
dt = t-t0
mag = dat['magpsf']
emag = dat['sigmamagpsf']
instr = dat['instrument']
filt = dat['filter']

# To define time bins, use the u-band observations and ignore P60.
use = np.logical_and(filt=='u', instr!='P60+SEDM')
dtbins = dt[use]

# for each bin, plot all photometry within half an hour of that bin
# and interpolate the p48 light curves onto that epoch

for ii,dtbin in enumerate(dtbins):
    # choose the panel
    ax = axarr.flatten()[ii]
    xvals = []
    yvals = []
    eyvals = []

    #UVOT
    choose = np.abs(uvdt-dtbin)<0.05
    for jj in np.arange(sum(choose)):
        wl = bands[uvfilt[choose][jj]]
        f = uvflux[choose][jj]*1E3
        ef = uveflux[choose][jj]*1E3
        xvals.append(wl)
        yvals.append(f)
        eyvals.append(ef)
     
    #LT
    choose = np.logical_and(np.abs(dt-dtbin)<0.05, instr=='LT+IOO')
    for jj in np.arange(sum(choose)):
        wl = bands[filt[choose][jj]]
        f,ef = toflux(mag[choose][jj]-ext[filt[choose][jj]],emag[choose][jj])
        xvals.append(wl)
        yvals.append(f)
        eyvals.append(ef)
    #P48
    for ztffilt in ['r','g','i']:
        choose = np.logical_and(instr=='P48+ZTF', filt==ztffilt)        
        ztfmag = np.interp(dtbin, dt[choose], mag[choose])
        f,ef = toflux(ztfmag-ext[ztffilt],emag=0.1) 
        xvals.append(bands[ztffilt])
        yvals.append(f)
        eyvals.append(ef)
    ax.errorbar(xvals, yvals, yerr=eyvals, c='k', fmt='.')
    txt = "$\Delta$t=" + str(np.round(dtbin,2))
    ax.text(0.9, 0.9, txt, transform=ax.transAxes,
            horizontalalignment='right', verticalalignment='top')

    # Fit a blackbody 1000 times
    xvals = np.array(xvals)
    yvals = np.array(yvals)
    eyvals = np.array(eyvals)
    nsim = 1000
    temps = np.zeros(nsim)
    radii = np.zeros(nsim)
    ysamples = np.zeros((nsim, len(xvals)))
    for jj,val in enumerate(yvals):
        ysamples[:,jj] = np.random.normal(loc=val,scale=eyvals[jj],size=nsim)
    for jj in np.arange(nsim):
        popt, pcov = curve_fit(
                bb_func, 3E10/(xvals*1E-8), ysamples[jj], p0=[5000,1E14])
        temps[jj] = popt[0]
        radii[jj] = popt[1]
        xplot = np.linspace(2000,8000)
        yplot = bb_func(3E10/(xplot*1E-8), popt[0], popt[1])
        #ax.plot(xplot,yplot,lw=0.1,alpha=0.1)
        plot(xplot,yplot,lw=0.1,alpha=0.1)
    lums = 4*np.pi*radii**2 * (5.67E-5)*temps**4
    T = np.mean(temps)
    eT = np.std(temps)
    R = np.mean(radii)
    eR = np.std(radii)
    L = np.mean(lums)
    eL = np.std(lums)

    print(dtbin)
    print("%s +/- %s" %(T/1E3, eT/1E3))
    print("%s +/- %s" %(R/1E14, eR/1E14))
    print("%s +/- %s" %(L/1E42, eL/1E42))

    # Plot the blackbody
    #xplot = np.linspace(2000,8000)
    #yplot = bb_func(3E10/(xplot*1E-8), T, R)
    #ax.plot(xplot,yplot,lw=0.5,alpha=1,c='k')

axarr[0,0].set_xlim(1E3, 2E4)
axarr[0,0].set_ylim(10, 5000)
axarr[0,0].set_yscale('log')
axarr[0,0].set_xscale('log')
plt.subplots_adjust(wspace=0,hspace=0)
#fig.text(0.5,0.04,"Wavelength (AA)", ha='center',fontsize=14)
fig.text(0.04,0.5,r'Flux ($\mu$Jy)',fontsize=14,verticalalignment='center',
        horizontalalignment='center',rotation='vertical')
axarr[1,2].set_xlabel(r'Wavelength (AA)',fontsize=14)
#plt.tight_layout()

#plt.show()
plt.savefig("bbfits.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
