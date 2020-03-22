""" A collage of early light curves for ZTF objects, 
for the rate section of the paper 
"""

import matplotlib.pyplot as plt
import numpy as np
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from astropy.io import ascii
from astropy.cosmology import Planck15
from astropy.time import Time
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF_fast_transient_search/code")
from forced_phot.run_forced_phot import get_forced_phot


def plot_2006aj(ax):
    """ Plot SN2006aj early (t < 4 days) light curve """
    dat = ascii.read("../../data/lc_060218_swift.dat")
    t = dat['t']
    dm = Planck15.distmod(z=0.033).value
    mag = dat['mag']
    emag = dat['emag']
    ax.plot(
        t, mag-dm, linestyle='-', lw=1, color='grey')


def plot_980425(ax):
    """ Plot GRB 980425 early (t < 4 days) light curve """
    dm = Planck15.distmod(z=0.0085).value
    dat = ascii.read(
        "../../data/lc_980425.txt", format="fixed_width", delimiter="&")
    jd_raw = np.array(dat['JD'])
    delta = jd_raw[0] - 17/24 # first obs is 17 hrs after burst
    jd = jd_raw - delta
    band = 'Rc'
    lc = dat[band]
    bad = np.array(['nodata' in val for val in lc])
    mag = np.array([float(val.split("$\\pm$")[0]) for val in lc[~bad]])
    mag_err = np.array(
            [float(val.split("$\\pm$")[1]) for val in lc[~bad]])
    t = jd[~bad]
    ax.plot(t,mag-dm,c='blue', linestyle='-')

def plot_19aaxfcpq(ax):
    mw_ext_g = 0.043 # g-band
    mw_ext_r = 0.030 # r-band
    z = 0.038
    dm = Planck15.distmod(z=z).value
    t0 = 2458638

    name = 'ZTF19aaxfcpq'
    #ra = 240.861997 
    #dec = 38.184068
    #jdobs = 2458638.7611
    #zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,5])
    #ascii.write(
    #        [filt,jd,mag,emag], 'ZTF19aaxfcpq.txt', 
    #        names=['filt', 'JD', 'mag', 'emag'])
    dat = ascii.read('ZTF19aaxfcpq.txt')
    filt = dat['filt']
    mag = dat['mag']
    emag = dat['emag']
    jd = dat['JD']
    choose = np.logical_and(filt=='g', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose],mag[choose]-mw_ext_g-dm,yerr=emag[choose], 
            c='k', fmt='o')
    choose = np.logical_and(filt=='r', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose], mag[choose]-mw_ext_r-dm, yerr=emag[choose], 
            mec='red', mfc='white', fmt='s', zorder=0)
    ax.text(0.01, 0.9, "ZTF19aaxfcpq (g\&r)", transform=ax.transAxes,
            fontsize=11)

    # The LT point
    ax.errorbar(
            2458639.4464-t0,19.27-mw_ext_g-dm,yerr=0.02, c='k', fmt='o')
    ax.errorbar(
            2458639.4464-t0, 19.38-mw_ext_r-dm, yerr=0.02, 
            mec='red', mfc='white', fmt='s', zorder=0)

    # the last g-band upper limit
    x = 2458637.7621-t0
    y = 20.78-mw_ext_g-dm
    print(y)
    ax.scatter(x, y, marker='o', c='k')
    ax.arrow(x, y, 0, 0.5, length_includes_head=True, 
            head_width=0.2, color='k', head_length=0.2)


def plot_19abqwtfu(ax):
    mw_ext_g = 0.641 # g-band
    mw_ext_r = 0.444
    z = 0.014353
    dm = Planck15.distmod(z=z).value
    t0 = 2458717

    # name = 'ZTF19abqwtfu'
    # ra = 346.829550 
    # dec = 13.855965
    # jdobs = t0
    # zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,5])
    # ascii.write(
    #         [filt,jd,mag,emag], 'ZTF19abqwtfu.txt', 
    #         names=['filt', 'JD', 'mag', 'emag'], overwrite=True)
    dat = ascii.read('ZTF19abqwtfu.txt')
    filt = dat['filt']
    mag = dat['mag']
    emag = dat['emag']
    jd = dat['JD']
    choose = np.logical_and(filt=='g', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose],mag[choose]-mw_ext_g-dm,yerr=emag[choose], 
            c='k', fmt='o', zorder=1)
    choose = np.logical_and(filt=='r', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose], mag[choose]-mw_ext_r-dm, yerr=emag[choose], 
            mec='red', mfc='white', fmt='s', zorder=0)
    ax.text(0.01, 0.9, "ZTF19abqwtfu (r\&g)", transform=ax.transAxes,
            fontsize=11)

    # The P60 points
    ax.errorbar(
            np.array([2458717.66, 2458718.96])-t0,
            np.array([18.74, 18.61])-mw_ext_r-dm,
            [0.17, 0.03], fmt='s', mec='red', mfc='white', zorder=0)


def plot_19abupned(ax):
    name = 'ZTF19abupned'
    mw_ext_g = 0.188
    mw_ext_r = 0.130
    z = 0.05
    dm = Planck15.distmod(z=z).value
    t0 = 2458724.9013-0.8

    # ra = 358.250149 
    # dec = 25.121270
    # zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,t0,[5,5])
    # ascii.write(
    #         [filt,jd,mag,emag], '%s.txt' %name, 
    #         names=['filt', 'JD', 'mag', 'emag'], overwrite=True)
    dat = ascii.read('%s.txt' %name)
    filt = dat['filt']
    mag = dat['mag']
    emag = dat['emag']
    jd = dat['JD']
    choose = np.logical_and(filt=='g', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose],mag[choose]-mw_ext_g-dm,yerr=emag[choose], 
            c='k', fmt='o', zorder=1)
    choose = np.logical_and(filt=='r', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose], mag[choose]-mw_ext_r-dm, yerr=emag[choose], 
            mec='red', mfc='white', fmt='s', zorder=0, c='red')
    ax.text(0.01, 0.9, "%s (r\&g)" %name, transform=ax.transAxes,
            fontsize=11)

    # the last g-band upper limit
    x = 2458723.9464-t0
    y = 20.99-mw_ext_g-dm
    ax.scatter(x, y, marker='o', c='k')
    ax.arrow(x, y, 0, 0.4, length_includes_head=True, 
            head_width=0.2, color='k', head_length=0.2)


def plot_20aaiqiti(ax):
    name = 'ZTF20aaiqiti'
    mw_ext_g = 0.046
    mw_ext_r = 0.032
    z = 0.025
    dm = Planck15.distmod(z=z).value
    t0 = 2458873.1

    # ra = 183.020377 
    # dec = 32.733871
    # zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,t0,[5,5])
    # ascii.write(
    #         [filt,jd,mag,emag], '%s.txt' %name, 
    #         names=['filt', 'JD', 'mag', 'emag'], overwrite=True)
    dat = ascii.read('%s.txt' %name)
    filt = dat['filt']
    mag = dat['mag']
    emag = dat['emag']
    jd = dat['JD']
    choose = np.logical_and(filt=='g', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose],mag[choose]-mw_ext_g-dm,yerr=emag[choose], 
            c='k', fmt='o', zorder=1)
    choose = np.logical_and(filt=='r', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose], mag[choose]-mw_ext_r-dm, yerr=emag[choose], 
            mec='red', mfc='white', fmt='s', zorder=0, c='red')
    ax.text(0.01, 0.9, "%s (r\&g)" %name, transform=ax.transAxes,
            fontsize=11)

    # the last g-band upper limit
    x = 2458872.9224-t0
    y = 20.39-mw_ext_g-dm
    ax.scatter(x, y, marker='o', c='k')
    ax.arrow(x, y, 0, 0.5, length_includes_head=True, 
            head_width=0.2, color='k', head_length=0.2)


def plot_19abqshry(ax):
    name = 'ZTF19abqshry'
    mw_ext_g = 0.044
    mw_ext_r = 0.031
    z = 0.030821
    dm = Planck15.distmod(z=z).value
    t0 = 2458715

    # ra = 249.638354
    # dec = 45.631184
    # zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,t0,[5,5])
    # ascii.write(
    #         [filt,jd,mag,emag], '%s.txt' %name, 
    #         names=['filt', 'JD', 'mag', 'emag'], overwrite=True)
    dat = ascii.read('%s.txt' %name)
    filt = dat['filt']
    mag = dat['mag']
    emag = dat['emag']
    jd = dat['JD']
    choose = np.logical_and(filt=='g', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose],mag[choose]-mw_ext_g-dm,yerr=emag[choose], 
            c='k', fmt='o', zorder=1)
    # MSIP point
    ax.errorbar(
            2458715.6635-t0, 20.73-mw_ext_g-dm, 0.22,
            c='k', fmt='o', zorder=1)
    choose = np.logical_and(filt=='r', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose], mag[choose]-mw_ext_r-dm, yerr=emag[choose], 
            mec='red', mfc='white', fmt='s', zorder=0, c='red')
    ax.text(0.01, 0.9, "%s (r\&g)" %name, transform=ax.transAxes,
            fontsize=11)

    # The P60 points
    ax.errorbar(2458717.7442-t0, 19.97-dm-mw_ext_g, 0.12, fmt='o', c='k')
    ax.errorbar(2458717.7415-t0, 19.88-dm-mw_ext_r, 0.17, 
            mec='red', mfc='white', c='red')

    # the last g-band upper limit
    x = 2458714.7238-t0
    y = 20.18-mw_ext_g-dm
    ax.scatter(x, y, marker='o', c='k')
    ax.arrow(x, y, 0, 0.3, length_includes_head=True, 
            head_width=0.2, color='k', head_length=0.15)


def plot_19ablesob(ax):
    name = 'ZTF19ablesob'
    mw_ext_g = 0.207
    mw_ext_r = 0.143
    z = 0.0558
    dm = Planck15.distmod(z=z).value
    t0 = 2458695.2

    # ra = 358.941433 
    # dec = 21.955484
    # zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,t0,[5,5])
    # ascii.write(
    #         [filt,jd,mag,emag], '%s.txt' %name, 
    #         names=['filt', 'JD', 'mag', 'emag'], overwrite=True)
    dat = ascii.read('%s.txt' %name)
    filt = dat['filt']
    mag = dat['mag']
    emag = dat['emag']
    jd = dat['JD']
    choose = np.logical_and(filt=='g', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose],mag[choose]-mw_ext_g-dm,yerr=emag[choose], 
            c='k', fmt='o', zorder=1)

    choose = np.logical_and(filt=='r', mag<99)
    dt = jd-t0
    ax.errorbar(
            dt[choose], mag[choose]-mw_ext_r-dm, yerr=emag[choose], 
            mec='red', mfc='white', fmt='s', zorder=0, c='red')
    ax.text(0.01, 0.9, "%s (r\&g)" %name, transform=ax.transAxes,
            fontsize=11)

    # MSIP points
    ax.errorbar(
            2458697.98-t0, 19.32-mw_ext_g-dm, 0.22,
            c='k', fmt='o', zorder=1)
    ax.errorbar(
            2458697.94-t0, 19.37-mw_ext_r-dm, yerr=0.12, 
            mec='red', mfc='white', fmt='s', zorder=0, c='red')


    # # The P60 points
    # ax.errorbar(2458717.7442-t0, 19.97-dm-mw_ext_g, 0.12, fmt='o', c='k')
    # ax.errorbar(2458717.7415-t0, 19.88-dm-mw_ext_r, 0.17, 
    #         mec='red', mfc='white', c='red')

    # the last g-band upper limit
    x = 2458694.9546-t0
    y = 21.12-mw_ext_g-dm
    ax.scatter(x, y, marker='o', c='k')
    ax.arrow(x, y, 0, 0.3, length_includes_head=True, 
            head_width=0.2, color='k', head_length=0.15)


#start, let's make a nice plot of two light curves: SN1998bw and SN2006aj
fig,axarr = plt.subplots(3,2,figsize=(7,6), sharex=True, sharey=True)
ax = axarr[0,0]
plot_19aaxfcpq(ax)
ax.set_xlim(-0.5,4)
ax.set_ylim(-14.5,-18.5)
plot_2006aj(ax)
plot_980425(ax)

ax = axarr[0,1]
ax.axis('off')
# plot_19abqwtfu(ax)
# ax.set_ylim(-15.5,-16.5)
 
ax = axarr[1,0]
plot_19abupned(ax)
 
ax = axarr[1,1]
plot_20aaiqiti(ax)

ax = axarr[2,0]
plot_19abqshry(ax)

ax = axarr[2,1]
plot_19ablesob(ax)

for ax in axarr.flatten():
    ax.tick_params(axis='both', labelsize=12)

fig.text(0.5, 0.04, '$\Delta$ t (days)', ha='center', fontsize=14)
fig.text(0.04, 0.5, 'Abs Mag', va='center', rotation='vertical', fontsize=14)
fig.subplots_adjust(
        left=0.16, bottom=0.13, right=None, top=None, wspace=0.3, hspace=0)
#plt.savefig("ztf_early_lc_collage.png", dpi=300, bbox_inches='tight')

plt.show()
