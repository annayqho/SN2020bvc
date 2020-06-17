""" Plot full light curves, one panel per band
"""


import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
import glob
import extinction
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/SN2020bvc/code")
from get_lc import *

zp = 2458883.17


def merge_lc():
    # get optical light curves 
    t, mag, emag, maglim, filt, inst = get_opt_lc()
    dt = t-zp
    choose = emag < 99
    dt = dt[choose]
    mag = mag[choose]
    emag = emag[choose]
    filt = filt[choose]

    # add the UV light curves
    add_dt, add_filt, fnu_mjy, efnu_mjy = get_uv_lc()
    # convert to AB mag
    add_mag = -2.5 * np.log10(fnu_mjy*1E-3) + 8.90
    add_emag = (efnu_mjy/fnu_mjy) # I think it's just the ratio
    choose = add_emag < 50
    dt = np.append(dt, add_dt[choose])
    filt = np.append(filt, add_filt[choose])
    mag = np.append(mag, add_mag[choose])
    emag = np.append(emag, add_emag[choose])

    return dt, filt, mag, emag


def plot_lc():
    dt, filt, mag, emag = merge_lc()
    det = np.logical_and(mag<99, ~np.isnan(mag))
    nondet = np.logical_or(mag==99, np.isnan(mag))

    fig,axarr = plt.subplots(
            4, 3, figsize=(8,8), sharex=True, sharey=True)

    bands = ['UVW2', 'UVM2', 'UVW1', 'U', 'u', 'B', 'g', 'V', 'r', 'i', 'z']
    for ii,use_f in enumerate(bands):
        ax = axarr.reshape(-1)[ii]
        choose = np.logical_and(det, filt == use_f)
        order = np.argsort(dt[choose])
        ax.errorbar(
                dt[choose][order], mag[choose][order]-ext[use_f], 
                emag[choose][order], c='black', fmt='o', ms=3,
                alpha=1.0, zorder=5)
        # for each panel, plot all of them as a grey background
        for f in bands:
            choose = np.logical_and(det, filt == f)
            order = np.argsort(dt[choose])
            ax.plot(
                    dt[choose][order], mag[choose][order]-ext[f], 
                    c='lightgrey', alpha=0.7, zorder=0)

        # for each panel, show the last non-detection from ATLAS (which was o-band)
        ax.arrow(
                2458883.17-zp, 19.4, 0, 1.0, length_includes_head=True,
                head_width=2, head_length=0.4, fc='k', ec='k', zorder=10)
        # for each panel, show the r-band peak with a cross
        ax.scatter(
                2458897-zp, 16.3, marker='+', c='r', zorder=10) 

        ax.yaxis.set_tick_params(labelsize=14)
        ax.xaxis.set_tick_params(labelsize=14)

        # for each panel, also show absolute mag
        # if ii % 3 == 2:
        #     ax2 = ax.twinx()
        #     y_f = lambda y_i: y_i-Planck15.distmod(z=0.03154).value
        #     ymin, ymax = ax.get_ylim()
        #     ax2.set_ylim((y_f(ymin), y_f(ymax)))
        #     ax2.plot([],[])
        #     ax2.tick_params(axis='both', labelsize=14)

        # label with the band
        if use_f == 'r':
            usecol = 'red'
            wid = 1.0
        else:
            usecol = 'k'
            wid = 0.5
        ax.text(
                0.04, 0.04, "$%s$" %use_f, 
                fontsize=14, transform=ax.transAxes,
                horizontalalignment='left',
                verticalalignment='bottom',
                bbox=dict(
                    boxstyle="round", fc='white', ec=usecol, 
                    lw=wid, alpha=1.0, pad=0.4))


    # Final reconfiguring
    axarr.reshape(-1)[-1].set_visible(False)
    plt.subplots_adjust(hspace=0, wspace=0)
    ax.set_xlim(-2, 32)
    ax.set_ylim(16, 22.5)
    ax.invert_yaxis()
    fig.text(0.5, 0.04, 
        r"Days since $t_0=$(UT 2020 Feb 03.67)", 
        ha='center', fontsize=16) 
    fig.text(0.04, 0.5, 'Apparent Mag (AB)', fontsize=16, rotation='vertical')

    plt.savefig("full_lc.png", dpi=500, bbox_inches='tight')
    #plt.show()


if __name__=="__main__":
    plot_lc()
