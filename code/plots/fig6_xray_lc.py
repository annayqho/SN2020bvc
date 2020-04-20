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
    dt = np.array([3, 13, 25])
    l = np.array([1.4E41, 1.5E40, 1.5E40])
    el = np.array([[0.5E41, 3.1E41],[0.7E40, 2.8E40], [0.7E40, 2.8E40]])

    for ii,dt_val in enumerate(dt):
        ax.errorbar(dt[ii], l[ii], yerr=[el[ii]], c=dark, marker='s')
    ax.plot(dt,l,c='k',lw=5)


def plot_98bw(ax, background=False):
    """ Plot GRB 980425 X-ray light curve
    I scraped the data from Alessandra's iPTF17cw figure
    """
    dat = np.loadtxt(ddir + "/sn1998bw.txt", delimiter=',')
    t = dat[:,0] 
    lum = dat[:,1]
    ax.plot(t, lum, c='#e55c30', label="98bw")


def plot_17iuk(ax, background=False):
    """ Plot SN 2017iuk X-ray light curve
    I scraped this from  the D'Elia 2018 paper, where values
    come in mJy at 1 keV"""
    dcm_17iuk = Planck15.luminosity_distance(z=0.0368).cgs.value
    dat = np.loadtxt(ddir + "/sn2017iuk.txt", delimiter=',')
    t = dat[:,0] / 86400
    lum = dat[:,1] * 1E-3 * 1E-23 * 4 * np.pi * dcm_17iuk**2 * 2.4E17
    order = np.argsort(t)
    ax.plot(t[order], lum[order], ls=':', c='#84206b', label="17iuk")


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
    t = dat[:,0] 
    lum = dat[:,1]
    ax.scatter(t, lum, edgecolor='grey', facecolor='white', marker='D',
            lw=0.5, label="10bh")


def plot_06aj(ax, background=False):
    """ Plot GRB 060218 X-ray light curve
    This data is from Campana (2006)

    """
    z = 0.032
    dcm = Planck15.luminosity_distance(z=z).cgs.value
    dat = np.loadtxt(ddir + "/sn2006aj.txt", delimiter=',')
    t = dat[:,0] / 86400
    lum = (dat[:,1]) * 4 * np.pi * dcm**2
    ax.plot(t, lum, c='k', label="06aj", linestyle='--')


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
    fig,ax = plt.subplots(
            1, 1, figsize=(6,4), sharex=True, sharey=True)

    plot_source(ax)
    plot_06aj(ax, background=False)
    plot_0316d(ax, background=False)
    #plot_12ap(ax, background=False)
    #plot_17cw(ax, background=False)
    plot_98bw(ax, background=False)
    #plot_09bb(ax, background=False)
    plot_17iuk(ax, background=False)

    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1E-1, 40)
    ax.set_ylim(5E39, 1E45)
    ax.set_xlabel(r"$\Delta t$ (days)", fontsize=16) 
    ax.set_ylabel(
            r'$L_{0.3\mathrm{-}10\,\mathrm{keV}}$ (erg/s)', fontsize=16)
    ax.legend(fontsize=14, loc='upper right')
    plt.tight_layout()
    #plt.show()
    plt.savefig(
       "xray_lc_log.png", dpi=500, bbox_inches='tight', pad_inches=0.1)
