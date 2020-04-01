""" Plot the spectral sequence of ZTF18aaqjovh
This relies on having scripts from the following repository:
    https://github.com/annayqho/Spectra
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import ascii
from ztfquery import query
from ztfquery import marshal
import glob
from astropy.time import Time
sys.path.append("/Users/annaho/Github/Spectra")
from normalize import smooth_spec
from measure_snr import get_snr


z = 0.025235
SPEC_DIR = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/spec"
t0 = 2458883.17-2400000.5


def get_res(tel):
    """ Here, this means the width of a line in Angstroms """
    if tel == 'LT':
        res = 18 # Angstrom, res at central wavelength
        res = 30 # add a couple of Ang?
    elif tel == 'P200':
        res = 10 # determined by eye from the spectrum
        # basically, width of a galaxy emission line is 10 AA
        # and each pixel is 1 AA
    elif tel == 'Keck1':
        res = 7*2 # determined by eye from spectrum
        # width of a line is around 7 pixels
        # and each pixel is 2 Angstroms
    elif tel == 'NOT':
        # width of a line is around 8 pixels
        # and each pixel is around 2.63 Ang
        res = 8*2.63
    elif tel == 'DCT':
        # width of a line is around 7 pixels
        # and each pixel is 2.2 Ang
        res = 7*2.2
    elif tel == 'P60':
        res = 20
    else:   
        print("I don't have this telescope")
    return res


def download_spec():
    """ Download the spectra """
    # connect to databases
    m = marshal.MarshalAccess()
    print("Connected")

    # download light curves
    marshal.download_spectra('ZTF18aaqjovh')

    # I moved all the spectra to the directory SPEC_DIR

    # return filenames
    f = glob.glob(SPEC_DIR+ "/*.ascii")
    print(f)
    return f


def get_files(sind, eind):
    """ start_ind: starting index; end_ind: end index """
    files = np.array(glob.glob(SPEC_DIR + "/*.ascii"))
    dt = np.zeros(len(files))
    tels = []
    cols = np.array([""]*len(dt), dtype='U10')
    # Read in files, pull out the corresponding dates, and sort by date
    for ii,f in enumerate(files):
        tel = f.split("_")[2]
        tels.append(tel)
        alldat = open(f).readlines()
        if tel == 'LT': 
            for line in alldat:
                if 'DATE-OBS' in line:
                    obsdate = line.split(" ")[2].rstrip().strip('\'')
                    print(obsdate)
                    print(Time(obsdate, format='isot'))
                    t = Time(obsdate, format='isot').mjd
                    print("LT")
                    print(t-t0)
                    dt[ii] = t-t0
            cols[ii] = 'magenta'
        elif tel == 'P200':
            for line in alldat:
                if 'UTSHUT' in line:
                    obsdate = line[11:]
                    t = Time(obsdate, format='isot').mjd
                    print("P200")
                    print(t-t0)
                    dt[ii] = t-t0
            cols[ii] = 'lightblue'
        elif tel == 'Keck1':
            for line in alldat:
                if 'DATE_BEG' in line:
                    obsdate = line[13:32]
                    t = Time(obsdate, format='isot').mjd
                    print("Keck")
                    print(t-t0)
                    dt[ii] = t-t0
            cols[ii] = 'red'
        elif tel == 'NOT':
            for line in alldat:
                if '#' in line:
                    t = float(line[1:])
                    print("NOT")
                    print(t-t0)
                    dt[ii] = t-t0
            cols[ii] = 'green'
        elif tel == 'P60':
            for line in alldat:
                if 'UTC' in line:
                    obsdate = line[7:]
                    t = Time(obsdate, format='isot').mjd
                    print("P60")
                    dt[ii] = t-t0
                    print(t-t0)
            cols[ii] = 'black'
        else:
            print("couldn't find telescope")
            print(tel)
    order = np.argsort(dt)
    files_sorted = files[order]
    dt_sorted = dt[order]
    tel_sorted = np.array(tels)[order]
    cols = cols[order]
    return files_sorted[sind:eind], dt_sorted[sind:eind], tel_sorted[sind:eind]


def load_spec(f):
    lc = np.loadtxt(f)
    obs_wl = lc[:,0]
    # shift to rest wl
    wl = obs_wl / (1+z)
    f = lc[:,1]
    if tel == 'Keck':
        eflux = lc[:,3]
    else:
        # estimate uncertainty from scatter in a continuum region
        eflux = np.array([get_snr(wl, f, 6300, 6500)]*len(wl))
    ivar = 1/eflux**2
    return wl, f, ivar


def plot_smoothed_spec(ax, x, y, ivar, tel, epoch, ls='-', lw=0.5, c='black', label=None, text=True):
    """ plot the smoothed spectrum """
    res = get_res(tel)
    smoothed = smooth_spec(x, y, ivar, res*3)
    ax.plot(
            x, smoothed, c=c,
            drawstyle='steps-mid', lw=lw, ls=ls, alpha=1.0, label=label,
            zorder=10)
    dt_str = r"+%s\,d" %(
            str(np.round(epoch, 1)))
    if text:
        ax.text(
                x[-1]+100, smoothed[-1],  s=dt_str,
                horizontalalignment='left', verticalalignment='center',
                fontsize=12)
    return smoothed


def plot_spec(ax, x, y, tel, epoch):
    """ plot the spectrum """
    ax.plot(
            x, y, c='lightgrey',
            drawstyle='steps-mid', lw=1.0, alpha=1.0)
    return ax


# def get_98bw(ind):
#     """ Epochs are -2, +3, +29, +73 """
#     ddir = "/Users/annaho/Dropbox/Projects/Research/ZTF18aaqjovh/data/spectra"
#     inputf = glob.glob(ddir + "/*.txt")
#     dat = ascii.read(inputf[ind])
#     wl = dat['wavelength'] / (1+0.0085)
#     flux = dat['flux']
#     return wl, flux    


if __name__=="__main__":
    fig,ax = plt.subplots(figsize=(6,10))
    files, epochs, tels = get_files(0, 6)
    #files, epochs, tels = get_files(6, 13)
    nfiles = len(files)
    shift = [1, 1.3, 1.8, 2.3, 2.8, 3.3, 3.6, 4.0, 4.5, 5, 5.5, 6, 6.5]
    #bw_shift = [2.1, 3, 4, 6, 7.1]
    for ii,f in enumerate(files):
        tel = tels[ii]
        print(tel)
        dt = epochs[ii]
        print(dt)
        wl, flux, ivar = load_spec(f)
        # Remove tellurics
        if ii >= 3:
            # wl 7150 to 7300
            mask = np.logical_and(wl>7150, wl<7300)
            ivar[mask] = 0
        choose = np.logical_and(wl>3660, wl < 9000)
        scale = flux[wl>4100][0]
        shifted = flux/scale-shift[ii]
        plot_spec(ax, wl[choose], (shifted-shift[ii])[choose], tel, dt)
        plot_smoothed_spec(
                ax, wl[choose], (shifted-shift[ii])[choose], 
                ivar[choose], tel, dt, lw=2)
#     for ii in np.arange(5):
#         wcomp,fcomp = get_98bw(ii)
#         scale = fcomp[wcomp>4100][ii]
#         shifted = fcomp/scale-bw_shift[ii]
#         if ii == 0:
#             ax.plot(
#                     wcomp,shifted,c='k',lw=0.3,
#                     label="98bw at similar phase")
#         else:
#             ax.plot(
#                     wcomp,shifted,c='k',lw=0.3,
#                     label='_nolegend_')
# 
    plt.tick_params(axis='both', labelsize=14)
    plt.xlim(3660, 10140)
    # for the first set
    plt.ylim(-6.3 -0.5)
    # for the second set
    #plt.ylim(-6.5 -0.5)
    plt.xlabel(r"Rest Wavelength (\AA)", fontsize=16)
    plt.ylabel(r"Scaled $F_{\lambda}$ + const.", fontsize=16)
    ax.get_yaxis().set_ticks([])
    plt.tight_layout()
    #plt.show()
    #plt.savefig("spec_sequence_second.png", dpi=500, bbox_inches='tight')
    plt.savefig("spec_sequence_first.png", dpi=500, bbox_inches='tight')
