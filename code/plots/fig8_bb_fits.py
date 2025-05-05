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
sys.path.append("/Users/annaho/Dropbox/astro/papers/papers_complete/SN2020bvc/code")
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
    Blam = (2*h*c**2/wl**5) * (1/(np.exp(h*c/(wl*k*T))-1))
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



if __name__=="__main__":
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
    tot_red_chisq = []
    tot_chisq = []
    dof = []

    # Get non-P48 photometry
    optt,optmag,optemag,optmaglim,optfilt,optinstr = get_opt_lc()
    optdt = optt-t0
    choose = np.logical_and(optinstr != 'P48+ZTF', optemag<50) # only dets
    flux_opt,eflux_opt= toflux(optmag[choose], optemag[choose]) # uJy

    # Get P48 photometry
    dt_p48,mag_p48,emag_p48,filt_p48 = get_binned_p48(
            optdt,optmag,optemag,optfilt,optinstr)
    flux_p48,eflux_p48= toflux(mag_p48, emag_p48) # uJy

    # Get UVOT photometry
    uvt,uvdt,uvfilt,uvflux,uveflux = get_uv_lc() # mJy

    # Merge 
    dt = np.hstack((uvdt,dt_p48,optdt[choose]))
    f = np.hstack((uvflux*1E3,flux_p48,flux_opt))
    filt = np.hstack((uvfilt,filt_p48,optfilt[choose]))
    ufilt = np.unique(filt)

    # Initialize figure
    fig,axarr = plt.subplots(5, 5, figsize=(8,8), sharex=True, sharey=True)

    # Define time bins manually
    dtbins = np.array(
            [0.9, 1.36, 1.8, 2.8, 3.8, 4.74, 5.78, 6.27, 7.8, 9.1, 9.8, 10.75, 
             11.09, 11.77, 12.47, 15.49, 20, 21.75, 23.77, 25.65, 26.5, 28.73, 
             29.48])

    # Err fac
    fac = 0.4 # assign systematic uncertainty that is 50% of the flux

    yvals_len = []
    for ii,dtbin in enumerate(dtbins):
        # choose the panel
        ax = axarr.flatten()[ii]
        xvals = []
        yvals = []

        # Interpolate each filter
        for uf in ufilt:
            choose = filt == uf
            # Check that interpolation is actually appropriate
            if np.logical_and(min(dt[choose])<dtbin, max(dt[choose])>dtbin):
                xvals.append(bands[uf])
                new_f = np.interp(dtbin, dt[choose], f[choose])
                yvals.append(new_f)

        # Done loading photometry; now plot and measure bb
        xvals = np.array(xvals)
        yvals = np.array(yvals)
        eyvals = fac*yvals

        # Sort in order of wavelength
        order = np.argsort(xvals)
        xvals = xvals[order]
        yvals = yvals[order]
        eyvals = eyvals[order]
        yvals_len.append(len(yvals))


        # Ignore UVW2 after 2 days
        if dtbin > 2:
            keep = xvals > bands['UVW2']
            ax.scatter(xvals[~keep], yvals[~keep], c='lightgrey')
            xvals = xvals[keep]
            yvals = yvals[keep]
            eyvals = eyvals[keep]

        ax.errorbar(
                xvals, yvals, yerr=eyvals, fmt='.', 
                mfc='lightgrey', mec='k', ecolor='k', elinewidth=0.3, ms=10)
        txt = "$\Delta t \\approx $" + str(np.round(dtbin,1))
        ax.text(0.95, 0.95, txt, transform=ax.transAxes,
                horizontalalignment='right', verticalalignment='top')

        # Fit a blackbody 600 times
        nsim = 600
        temps = np.zeros(nsim)
        radii = np.zeros(nsim)
        ysamples = np.zeros((nsim, len(xvals)))
        if dtbin < 1:
            t0 = 13000
        else:
            # initialize temperature guess at wien's value
            t0 = 0.0029/(xvals[np.argmax(yvals)]*1E-10)
        # initialize radius guess
        dcm = Planck15.luminosity_distance(z=0.0252).cgs.value
        L0 = max(yvals)*1E-6*1E-23*4*np.pi*dcm**2*3E10/(xvals[np.argmax(yvals)]*1E-8)
        r0 = (L0 / (4 * np.pi * (5.67E-5) * t0**4))**0.5
        for jj,val in enumerate(yvals):
            ysamples[:,jj] = np.random.normal(
                    loc=val,scale=eyvals[jj],size=nsim)
        for jj in np.arange(nsim):
            popt, pcov = curve_fit(
                    bb_func, xvals*1E-8, ysamples[jj], p0=[t0,r0],
                    bounds=([0,0],[30000, np.inf]),maxfev=100000)
            temps[jj] = popt[0]
            radii[jj] = popt[1]
            xplot = np.linspace(1000,20000)
            yplot = bb_func(xplot*1E-8, popt[0], popt[1])
            ax.plot(xplot,yplot,lw=0.5,alpha=0.1,c='lightgrey')
        lums = 4*np.pi*radii**2 * (5.67E-5)*temps**4

        # Sort results and calculate 16-to-84 percentile range
        print(dtbin)

        nvals = len(temps)
        start_ind = int(0.16*nvals)
        end_ind = int(0.84*nvals)
        temps = np.sort(temps)
        radii = np.sort(radii)
        lums = np.sort(lums)
        L = np.median(lums)
        Lbol.append(L)
        low = L-lums[start_ind]
        Lbol_lo.append(low)
        hi = lums[end_ind]-L
        Lbol_hi.append(hi)
        print("%s +%s -%s" %(L/1E42, hi/1E42, low/1E42))
        T = np.median(temps)
        Teff.append(T)
        low = T-temps[start_ind]
        Teff_lo.append(low)
        hi = temps[end_ind]-T
        Teff_hi.append(hi)
        print("%s +%s -%s" %(T/1E3, hi/1E3, low/1E3))
        R = np.median(radii)
        Rph.append(R)
        low = R-radii[start_ind]
        Rph_lo.append(low)
        hi = radii[end_ind]-R
        Rph_hi.append(hi)
        print("%s +%s -%s" %(R/1E14, hi/1E14, low/1E14))

        # Calculate the chi squared of the final fit
        chisq = sum((yvals-bb_func(xvals*1E-8, T, R))**2/eyvals**2)
        tot_chisq.append(chisq)
        dof_val = len(yvals)-1
        dof.append(dof_val)
        red_chisq = chisq/dof_val
        print(chisq,dof_val)
        tot_red_chisq.append(red_chisq)

        xplot = np.linspace(1000,20000)
        yplot = bb_func(xplot*1E-8, T, R)
        ax.plot(xplot,yplot,lw=0.5,alpha=1,c='Crimson')

    Lbol = np.array(Lbol)
    np.savetxt('lbol.txt', np.array([dtbins,Lbol,Lbol_hi,Lbol_lo]).T)
    Rph = np.array(Rph)
    Teff = np.array(Teff)
    tot_red_chisq = np.array(tot_red_chisq)
    tot_chisq = np.array(tot_chisq)
    dof = np.array(dof)

    axarr[0,0].set_xlim(1E3, 2E4)
    axarr[0,0].set_ylim(3, 5000)
    axarr[0,0].set_yscale('log')
    axarr[0,0].set_xscale('log')
    axarr[4,3].axis('off')
    axarr[4,4].axis('off')
    plt.subplots_adjust(wspace=0,hspace=0)
    #fig.text(0.5,0.04,"Wavelength (AA)", ha='center',fontsize=14)
    fig.text(0.04,0.5,r'Flux ($\mu$Jy)',fontsize=14,verticalalignment='center',
            horizontalalignment='center',rotation='vertical')
    axarr[4,2].set_xlabel(r'Wavelength (\AA)',fontsize=14)
    #plt.tight_layout()

    #plt.show()
    plt.savefig("bbfits.eps", dpi=300, bbox_inches='tight')
    plt.close()
