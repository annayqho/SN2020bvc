""" Plot a grid of radio light curves """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
import glob

dark = '#140b34'
purp = '#84206b'
orag = '#e55c30'
yell = '#f6d746'

ddir = "/Users/annaho/Dropbox/Projects/Research/IcBL/data/xray_compilations"

dcm = Planck15.luminosity_distance(z=0.02507).cgs.value

def plot_source(ax):
    dt = np.array([2.5, 14.3, 26.6])
    f = np.array([1.9E-13, 1.1E-14, 7.5E-15])
    l = 4 * np.pi * dcm**2 * f
    ax.errorbar(dt[0], l[0], yerr=(1.4E-13)*4*np.pi*dcm**2, c=dark, marker='s')
    ax.errorbar(dt[1], l[1], 
        yerr=[[l[1]-(0.8E-14)*4*np.pi*dcm**2, l[1]+(0.9E-14)*4*np.pi*dcm**2]], 
        c=dark, marker='s')
    ax.errorbar(dt[2], l[2], 
            yerr=[[l[2]-4.9E-15*4*np.pi*dcm**2, l[2]+7.2E-15*4*np.pi*dcm**2]], 
            c=dark, marker='s')
    ax.plot(dt,l,c='k',lw=2)


def plot_98bw(ax, background=False):
    """ Plot GRB 980425 X-ray light curve
    I scraped the data from Alessandra's iPTF17cw figure
    """
    dat = np.loadtxt(ddir + "/sn1998bw.txt", delimiter=',')

    col = 'lightgrey'
    if background is False:
        col = dark

    t = dat[:,0] 
    lum = dat[:,1]

    if background is False:
        col = '#e55c30'
    ax.plot(t, lum, c=col, label="_nolegend_")
    #ax.scatter(t, lum, c=col, marker='.', label="_nolegend_") 

    if background==False:
        ax.text(0.1, 0.1, "SN1998bw", fontsize=12, transform=ax.transAxes)


def plot_09bb(ax, background=False):
    """ Plot X-ray data from Corsi+ 2017
    """

    col = 'lightgrey'
    if background is False:
        col = dark

    ax.scatter(31, 4.5E39, c=col, marker='o')

    if background is False:
        col = '#e55c30'
    ax.scatter(31, 4.5E39, c=col, label="_nolegend_", marker='o')

    if background==False:
        ax.text(0.1, 0.1, "SN2009bb", fontsize=12, transform=ax.transAxes)


def plot_0316d(ax, background=False):
    """ Plot GRB 060218 X-ray light curve
    I scraped the data from Alessandra's iPTF17cw figure
    """
    z = 0.0593
    dcm = Planck15.luminosity_distance(z=z).cgs.value
    dat = np.loadtxt(ddir + "/sn2010bh.txt", delimiter=',')

    col = 'lightgrey'
    if background is False:
        col = dark

    t = dat[:,0] 
    lum = dat[:,1]

    if background is False:
        col = '#e55c30'
    ax.plot(t, lum, c=col, label="_nolegend_")
    #ax.scatter(t, lum, c=col, marker='.', label="_nolegend_") 

    if background==False:
        ax.text(0.1, 0.1, "SN2010bh", fontsize=12, transform=ax.transAxes)


def plot_06aj(ax, background=False):
    """ Plot GRB 060218 X-ray light curve
    This data is from Campana (2006)

    """
    z = 0.032
    dcm = Planck15.luminosity_distance(z=z).cgs.value
    dat = np.loadtxt(ddir + "/sn2006aj.txt", delimiter=',')

    col = 'lightgrey'
    if background is False:
        col = dark

    t = dat[:,0] / 86400
    lum = (dat[:,1]) * 4 * np.pi * dcm**2

    if background is False:
        col = '#e55c30'
    ax.plot(t, lum, c=col, label="_nolegend_")
    #ax.scatter(t, lum, c=col, marker='.', label="_nolegend_") 

    if background==False:
        ax.text(0.1, 0.1, "SN2006aj", fontsize=12, transform=ax.transAxes)



def plot_17cw(ax, background=False):
    """ Plot X-ray data
    from Corsi+ 2017
    """

    col = 'lightgrey'
    if background is False:
        col = dark

    ax.scatter(39, 1.2E41, c=col, marker='o')

    if background is False:
        col = '#e55c30'
    ax.scatter(39, 1.2E41, c=col, label="_nolegend_", marker='o')

    if background==False:
        ax.text(0.1, 0.1, "iPTF17cw", fontsize=12, transform=ax.transAxes)


def plot_12ap(ax, background=False):
    """ Plot SN2012ap data
    """
    # From Margutti+ 2014

    col = 'lightgrey'
    if background is False:
        col = dark

    ax.scatter(24, 2.4E39, c=col, marker='o')
    ax.arrow(24, 2.4E39, 0, -2E39, length_includes_head=True,
            head_length=5E38, head_width=2, color=col)

    if background is False:
        col = '#e55c30'
    ax.scatter(24, 2.4E39, c=col, label="_nolegend_", marker='o')
    ax.arrow(24, 2.4E39, 0, -2E39, length_includes_head=True,
            head_length=5E38, head_width=2, color=col)

    if background==False:
        ax.text(0.1, 0.1, "SN2012ap", fontsize=12, transform=ax.transAxes)


if __name__=="__main__":
    fig,axarr = plt.subplots(
            3, 2, figsize=(6,6), sharex=True, sharey=True)

    for ii,ax in enumerate(axarr.reshape(-1)):
        # LC of ZTF18aaqjovh
        plot_source(ax)
        plot_06aj(ax, background=True)
        plot_0316d(ax, background=True)
        plot_12ap(ax, background=True)
        plot_17cw(ax, background=True)
        plot_98bw(ax, background=True)
        plot_09bb(ax, background=True)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.xaxis.set_tick_params(labelsize=12)
        ax.set_yscale('log')
        ax.set_xlim(0, 40)
        ax.set_ylim(1E38, 1E45)

    plot_06aj(axarr[2,0], background=False)
    plot_0316d(axarr[1,0], background=False)
    plot_12ap(axarr[1,1], background=False)
    plot_17cw(axarr[2,1], background=False)
    plot_98bw(axarr[0,0], background=False)
    plot_09bb(axarr[0,1], background=False)
    fig.text(0.5, 0.04, r"$\Delta t$ (days)", ha='center', fontsize=12) 
    fig.text(
            0.04, 0.5, r'0.3--10 keV X-ray Luminosity (erg/s)', 
            fontsize=12, rotation='vertical', horizontalalignment='center',
            verticalalignment='center')

    fig.subplots_adjust(wspace=0.1, hspace=0.1)

    #plt.show()
    plt.savefig(
       "xray_lc.png", dpi=500, bbox_inches='tight', pad_inches=0.1)
