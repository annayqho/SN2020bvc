""" Make my own bolometric light curve of SN1998bw, to test against the
Cano+2013 bolometric light curve """

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
from scipy.optimize import curve_fit
from astropy.io import ascii
data_dir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data"
dm = Planck15.distmod(z=0.0085).value
dat = ascii.read(data_dir + "/sn1998bw.dat", delimiter=';')
jd_all = (dat['JD'].data).data

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
    d = Planck15.luminosity_distance(z=0.0085).cgs.value
    h = 6.626E-27
    c = 3E10
    k = 1.38E-16
    Blam = (2*h*c**2/wl**5) * (1/(np.exp(h*c/(wl*k*T)-1)))
    flam = Blam * np.pi * R**2 / d**2
    # in units of uJy
    fnu = wl**2 * flam / c
    return fnu / 1E-23 / 1E-6


# Ultimately, here are the parameters I will measure
Lbol = []
Teff = []
Rph = []
Lbol_lo = []
Teff_lo = []
Rph_lo = []
Lbol_hi = []
Teff_hi = []
Rph_hi = []

# Only choose epochs with all five bands (it works to just select U-band)
# And we only want data from the first 30 days, i.e. jd < 962
ind = np.where(np.logical_and(dat['Umag'].data.data>0, jd_all<962))[0]
jd = jd_all[np.logical_and(dat['Umag'].data.data>0, jd_all<962)]

# There are 17 epochs
fig,axarr = plt.subplots(4,5,figsize=(8,8),sharex=True,sharey=True)
xvals = np.array([3465, 4392, 5417, 6500, 8176])

for ii,jd_val in enumerate(jd):
    ax = axarr.flatten()[ii]
    # UVBRI
    yvalues_temp = np.array([dat['Umag'][ind[ii]], dat['Bmag'][ind[ii]], 
        dat['Vmag'][ind[ii]], dat['Rcmag'][ind[ii]], dat['Icmag'][ind[ii]]])
    eyvalues_temp = np.array([dat['e_Umag'][ind[ii]], dat['e_Bmag'][ind[ii]], 
        dat['e_Vmag'][ind[ii]], dat['e_Rcmag'][ind[ii]], 
        dat['e_Icmag'][ind[ii]]])
    yvals,eyvals = toflux(yvalues_temp, eyvalues_temp)
    ax.errorbar(xvals, yvals, yerr=eyvals, fmt='o')

    nsim = 10
    temps = np.zeros(nsim)
    radii = np.zeros(nsim)
    ysamples = np.zeros((nsim, len(xvals)))
    for jj,val in enumerate(yvals):
        ysamples[:,jj] = np.random.normal(loc=val,scale=eyvals[jj],size=nsim)
    for jj in np.arange(nsim):
        popt, pcov = curve_fit(
                bb_func, xvals*1E-8, ysamples[jj], p0=[10000,1E14],
                bounds=([2000,1E12],[20000, 2.5E15]))
        temps[jj] = popt[0]
        radii[jj] = popt[1]
        xplot = np.linspace(1000,20000)
        yplot = bb_func(xplot*1E-8, popt[0], popt[1])
        ax.plot(xplot,yplot,lw=0.1,alpha=0.1)
    lums = 4*np.pi*radii**2 * (5.67E-5)*temps**4

    nvals = len(temps)
    start_ind = int(0.16*nvals)
    end_ind = int(0.84*nvals)
    temps = np.sort(temps)
    T = np.median(temps)
    Teff.append(T)
    low = T-temps[start_ind]
    Teff_lo.append(low)
    hi = temps[end_ind]-T
    Teff_hi.append(hi)
    print("%s +%s -%s" %(T/1E3, hi/1E3, low/1E3))
    radii = np.sort(radii)
    R = np.median(radii)
    Rph.append(R)
    low = R-radii[start_ind]
    Rph_lo.append(low)
    hi = radii[end_ind]-R
    Rph_hi.append(hi)
    print("%s +%s -%s" %(R/1E14, hi/1E14, low/1E14))
    lums = np.sort(lums)
    L = np.median(lums)
    Lbol.append(L)
    low = L-lums[start_ind]
    Lbol_lo.append(low)
    hi = lums[end_ind]-L
    Lbol_hi.append(hi)
    print("%s +%s -%s" %(L/1E42, hi/1E42, low/1E42))

    xplot = np.linspace(1000,20000)
    yplot = bb_func(xplot*1E-8, T, R)
    ax.plot(xplot,yplot,lw=0.5,alpha=1,c='Crimson')

    txt = "$\Delta t \\approx $" + str(np.round(jd_val-928.5,1))
    ax.text(0.95, 0.95, txt, transform=ax.transAxes,
            horizontalalignment='right', verticalalignment='top')

axarr[0,0].set_xlim(2E3, 1.5E4)
axarr[0,0].set_ylim(1000, 20000)
axarr[0,0].set_yscale('log')
axarr[0,0].set_xscale('log')
plt.subplots_adjust(wspace=0,hspace=0)
fig.text(0.04,0.5,r'Flux ($\mu$Jy)',fontsize=14,verticalalignment='center',
    horizontalalignment='center',rotation='vertical')
axarr[1,2].set_xlabel(r'Wavelength (\AA)',fontsize=14)
plt.show()
