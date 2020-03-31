""" Fit blackbody SED but also spectrum on the first day """

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

fig,ax = plt.subplots(1,1,figsize=(7,4))


def bb_func(nu,T,R):
    d = Planck15.luminosity_distance(z=0.02507).cgs.value
    h = 6.63E-27
    c = 3E10
    k = 1.38E-16
    Inu = (2*h*nu**3/c**2) * (1/(np.exp(h*nu/(k*T)-1)))
    fnu = Inu * np.pi * R**2 / d**2
    return fnu 


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
nsim = 1000
temps = np.zeros(nsim)
radii = np.zeros(nsim)
ysamples = np.zeros((nsim, len(xvals)))
for jj,val in enumerate(yvals):
    ysamples[:,jj] = np.random.normal(loc=val,scale=eyvals[jj],size=nsim)
for jj in np.arange(nsim):
    popt, pcov = curve_fit(
            bb_func, xvals, ysamples[jj], p0=[5000,1E14])
    temps[jj] = popt[0]
    radii[jj] = popt[1]
choose = np.logical_and(temps > 0, radii > 0)
lums = 4*np.pi*radii[choose]**2 * (5.67E-5)*temps[choose]**4
T = np.mean(temps[choose])
eT = np.std(temps[choose])
R = np.mean(radii)
eR = np.std(radii)
L = np.mean(lums)
eL = np.std(lums)

print("%s +/- %s" %(T/1E3, eT/1E3))
print("%s +/- %s" %(R/1E14, eR/1E14))
print("%s +/- %s" %(L/1E42, eL/1E42))

xplt = np.linspace(3500,9500)
xplt_nu = 3E10/(xplt*1E-8)
yplt_nu = bb_func(xplt_nu, T, R)
yplt = xplt_nu * yplt_nu / xplt

T = '13.2'
R = '5.1'
L = '5.6'

plt.plot(wl/1.0252, f_lam/1E-15, c='k', drawstyle='steps-mid', lw=0.5, ls='-', alpha=1)
plt.plot(xplt,yplt/1E-15,lw=1,alpha=1,c='k')
plt.text(
    0.95, 0.95, "$\Delta t=0.67\,$d", fontsize=14, transform=ax.transAxes,
    horizontalalignment='right', verticalalignment='top')
plt.text(
    0.95, 0.85, 
    r"$T=%s \times 10^{3} $\,K, $R=%s\times10^{14}$ cm, $L=%s\times10^{42}\,$erg/s" %(T,R,L),
    horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,
    fontsize=14)

# Now plot the earliest spectrum of SN2017iuk, from Izzo+2019
dat = pyfits.open("../../data/171205A_731_SF1.fits")[0].data
wl = np.linspace(3200,5600,len(dat))
plt.plot(wl/1.0368, 2*dat/1E-15, c='k', drawstyle='steps-mid', lw=0.5, ls='-', alpha=1)


plt.ylabel(r"Flux ($10^{-15}$ erg\,s$^{-1}$\,cm$^{-2}\,$\AA$^{-1}$)", fontsize=16)
plt.xlabel("Rest Wavelength ($\AA$)", fontsize=16)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
#plt.show()
plt.savefig("bb_first_epoch.png", dpi=200)
