""" 
Plot of luminosity over time
"""


import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("../../../AT2018cow/code")
sys.path.append("../../../../projects_proto/IcBL/old/data/radio_compilations/Zauderer2011")
from astropy.cosmology import Planck15
from astropy.io import ascii
from get_radio import *
from scale_fluxes import sma_lc
from read_table import *


def plot_limits(ax, x, y, ratiox, ratioy, col):
    """ Plot two arrows from the point """
    ax.annotate('', xy=(x*ratiox, y), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))
    ax.annotate('', xy=(x, y*ratioy), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))


def plot_line(ax, d, t, nufnu, name, label, col, legend=False, zorder=1):
    """ Plot a line
    If nu > 90 GHz, make it on the left axis
    If nu < 10 GHz, make it on the right axis
    
    Parameters
    ----------
    nufnu: Hz * mJy 
    name: name of the source
    label: label to use as the legend (which also determines the col)
    """
    lum = nufnu * 1e-23 * 1e-3 * 4 * np.pi * d**2
    fs = 11
    nsize = 10 # normal size for points
    if name=='AT2018cow':
        marker='*'
        fcol = col
        s=70
    else:
        if label=='SN':
            marker='o'
            s=nsize
            fcol = col # fill color
            label = 'CC SN'
        if label=='Ic-BL SN':
            marker='o'
            s=nsize
            fcol = col # fill color
            label = 'Ic-BL SN'
        elif label=='GRB':
            marker='o'
            fcol = 'white' # unfilled
            s=nsize
            label = 'GRB'
        elif label=='Rel. SN':
            marker='s'
            fcol = col 
            s=nsize
            label = 'GRB-SN'
        elif label=='TDE':
            marker='s'
            fcol = 'white' #unfilled
            s=nsize
            label='TDE'
    ax.scatter(
            t, lum, facecolor=fcol, edgecolor=col, 
            marker=marker, s=s, zorder=zorder, label="_none")
    if legend:
        ax.plot(t, lum, c=col, ls='-', label=label, zorder=zorder)
    else:
        ax.plot(t, lum, c=col, ls='-', label="__none", zorder=zorder)
    return lum


def plot_points(ax, d, nu, t, f, marker, name=None):
    """ Plot set of two points """
    lums = []
    for ii,nuval in enumerate(nu):
        if nuval > 90E9:
            lum = plot_point(ax, d, nuval, t[ii], f[ii], marker, name=name)
            lums.append(lum)
        else:
            lum = plot_point(ax, d, nuval, t[ii], f[ii], marker)
            lums.append(lum)
    ax.plot(
        t, lums, ls='--', c='k', zorder=0)
    return lums


def koala(ax, col, legend):
    # This is the X-band light curve
    dt = np.array([81,310,352,396])
    f = np.array([8.3E39,1.5E39,1.1E39,8.16E38])-2.65E38
    ax.errorbar(dt/(1.2714), f, 0.0006*f, c='black', fmt='*', ms=10)
    ax.plot(dt/(1.2714), f, c='black', lw=2)
    ax.text(
            60, 1.2E40, "ZTF18abvkwla", 
            fontsize=11, horizontalalignment='center',
            verticalalignment='bottom')


def at2018cow(ax, col, legend):
    """ 231.5 GHz light curve and 9 GHz light curve """
    d = Planck15.luminosity_distance(z=0.014).cgs.value

    # high frequency
    a, b, c = sma_lc()
    dt, f, ef = b
    ef_comb = np.sqrt(ef**2 + (0.15*f)**2)
    nu = 231.5E9

    # low frequency
    nu = 9E9
    data_dir = "../../../AT2018cow/data"
    dat = Table.read(
        "%s/radio_lc.dat" %data_dir, delimiter="&",
        format='ascii.no_header')
    tel = np.array(dat['col2'])
    choose = np.logical_or(tel == 'SMA', tel == 'ATCA')

    days = np.array(dat['col1'][choose])
    freq = np.array(dat['col3'][choose]).astype(float)
    flux_raw = np.array(dat['col4'][choose])
    flux = np.array(
            [float(val.split("pm")[0][1:]) for val in flux_raw])
    eflux_sys = np.array([0.1*f for f in flux])
    eflux_form = np.array(
            [float(val.split("pm")[1][0:-1]) for val in flux_raw])
    eflux = np.sqrt(eflux_sys**2 + eflux_form**2)
    choose = freq == 9

    # add the Margutti point
    x = np.hstack((days[choose], 83.51))
    y = np.hstack((flux[choose], 9.1)) * nu
    lum = plot_line(
            ax, d, x, y,
            'AT2018cow', None, col, legend, zorder=10)
    ax.text(x[-1], lum[-1]*1.2, 'AT2018cow', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='center')


def maxi(ax):
    """ data from Dillon """
    d = 4.62E26
    t = 4*365
    nu = 6E9
    f = 6.16
    lum = plot_point(ax, d, nu, t, f, '*')
    plt.text(t/1.2, lum, 'MAXI 140814A', fontsize=11,
            horizontalalignment='right', verticalalignment='center')


def tde(ax, col, legend):
    """  Plot the 225/230 GHz light curve from the SMA
    
    Plot the 4.9 GHz light curve from the VLA
    """
    z = 0.354
    d = Planck15.luminosity_distance(z=z).cgs.value

    # In the Eftekhari paper, it says that although the event was first
    # triggered by Swift/BAT on 2011 March 28.55 UT, subsequent
    # analysis of the BAT data revealed discernible emission as early as
    # 2011 March 25. All times should therefore be shifted relative to Mar 25.5

    # Need to add 3.04 to the Zauderer points
    nu, dt, f, ef, islim = zauderer()
    t = (dt+3.04)/(1+z)

    # Low frequency
    nu_plt = 4.9E9
    choose = np.logical_and(~islim, nu == nu_plt/1E9)
    dt_all = t[choose]
    nufnu_all = nu_plt*f[choose]

    # adding the set from Berger2012
    # and making the same correction as above
    # this is 4.9 GHz
    t = (np.array([3.87, 4.76, 5.00, 5.79, 6.78, 7.77, 9.79, 14.98, 22.78,
        35.86, 50.65, 67.61, 94.64, 111.62, 126.51, 143.62, 164.38, 174.47,
        197.41, 213.32])) / (1+z)
    f = np.array([0.25, 0.34, 0.34, 0.61, 0.82, 1.48, 1.47, 1.80, 2.10, 4.62,
        4.84, 5.86, 9.06, 9.10, 9.10, 11.71, 12.93, 12.83, 13.29, 12.43])

    # Berger 2012: use the 8.4 GHz light curve, since that's closest in freq
    #t = (np.array([14.97, 127.69, 159.77, 174.47, 177.50, 197.41, 213.32, 219.22]))/(1+z)
    #f = np.array([5.49, 19.03, 22.15, 23.19, 23.65, 22.42, 22.04, 21.52])
    #dt_all = np.append(dt_all, t)
    #nufnu_all = np.append(nufnu_all, f*nu_plt)

    # adding the set from Zauderer2013
    # they also say it's relative to March 25.5...
    # so I think I need to subtract 3.04 days from here too
    t = (np.array([245.23, 302.95, 383.92, 453.66, 582.31]))/(1+z)
    f = np.array([12.17, 12.05, 12.24, 11.12, 8.90])
    dt_all = np.append(dt_all, t)
    nufnu_all = np.append(nufnu_all, f*nu_plt)

    # adding the set from Eftekhari 2018
    t = np.array([645, 651.1, 787.6, 1032, 1105, 1373, 1894])
    f = np.array([8.24, 8.63, 6.23, 4.21, 3.52, 2.34, 1.47])
    dt_all = np.append(dt_all, t)
    nufnu_all = np.append(nufnu_all, f*nu_plt)

    order = np.argsort(dt_all)
    lum = plot_line(
            ax, d, dt_all[order], nufnu_all[order], 
            'SwiftJ1644+57', 'TDE', col, legend)
    ax.text(dt_all[order][10], lum[10]*1.1, 'Swift J1644+57', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='left')


def asassn14li(ax, col, legend):
    """ Alexander et al. 2016 """
    nu = 5.0E9
    d = Planck15.luminosity_distance(z=0.0206).cgs.value
    t = np.array([80, 141.38, 207.33, 246.25, 303.01, 375.94, 389.96])
    flux = np.array([2, 1.91, 1.74, 1.56, 1.26, 0.81, 0.89])
    lum = plot_line(
            ax, d, t, nu*flux, 'ASASSN14li', 'TDE', col, legend,
            zorder=10)
    ax.text(t[-1]/2, lum[-1]/1.5, 'ASASSN14li', fontsize=11,
            verticalalignment='top',
            horizontalalignment='left')


def sn2003L(ax, col, legend):
    """ Soderberg et al
    Values at 8.5 GHz """
    d = 2.8432575937224894e+26
    nu_plt = 8.5E9
    nu, dt, f, ef = read_2003L()
    choose = nu == nu_plt
    lum = plot_line(
            ax, d, dt[choose], 1E-3*f[choose]*nu_plt, 
            'SN2003L', 'SN', col, legend)
    ax.text(dt[choose][-1]/1.05, lum[-1], 'SN2003L', fontsize=11,
            verticalalignment='center',
            horizontalalignment='center')
    

def sn1979c(ax, col, legend):
    """ Weiler 1986 and Weiler 1991
    This is a IIL 
    
    Too old to have the raw table on arXiv
    """
    d = 5.341805643483106e+25
    nu = 1.4E9 # 20cm
    t = np.array(
            [437,594,631,663,679,684,727,747,786,822,839,876,882,
             914,937,973,995,
             1026,1071,1091,1127,1156,1168,1212,1243,1277,1314,1358,1390,
             1415,1435,1466,1513,1565,1600,1634,1659,1698,1714,1750,1771,
             1931,2027])
    flux = np.array(
            [0.2,2.1,2.5,2.7,2.8,2.8,4.4,4.8,6.0,7.1,7.1,7.6,8.6,
             9.8,6.5,8.6,9.5,
             10.2,10.8,10.3,10.4,12.2,10.1,10.2,11.5,11.2,13.0,11.3,10.2,
             9.6,11.2,13.2,11.1,9.1,8.5,9.1,8.8,10.1,9.7,9.1,8.9,
             7.0,7.7])
    lum = plot_line(ax, d, t, nu*flux, 'SN1979c', 'SN', col, legend)
    ax.text(t[0]/1.05, lum[0], 'SN1979C', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')
    

def sn1993J(ax, col, legend):
    """ SN 1993J from Weiler et al. 
    This is the peak of the 99.4 GHz light curve
    There is also a 110 GHz point, but only one,
    so I can't get the peak.
    1.4 GHz / 20 cm, 4.9 GHz / 6 cm
    values come from the best-fit model,
    but by eye they are clearly pretty close
    """
    d = 1.1e25
    freq = 5E9
    nu, dt, f, ef, islim = read_1993J_low_freq()
    choose = np.logical_and(~islim, nu==freq)
    lum = plot_line(
            ax, d, dt[choose], freq*f[choose], 
            'SN1993J', 'SN', col, legend)
    ax.text(dt[choose][0]/1.05, lum[0], 'SN1993J', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')

    freq = 99.4E9
    nu, dt, f, ef, islim = read_1993J_high_freq()
    choose = np.logical_and(~islim, nu==freq)


def sn2011dh(ax, col, legend):
    """ SN 2011dh
    Horesh et al. 2013
    Krauss et al. 2012
    M51: d = 8.03 Mpc; expl date May 31.58
    """
    d = 2.5E25

    # HIGH FREQUENCY
    # use two freq: 107E9 and 93E9
    dt, nu, f, ef, islim = read_2011dh()
    choose = np.logical_and(~islim, np.logical_or(nu==107E9, nu==93E9))

    # LOW FREQUENCY
    # from Horesh 2013 and Krauss 2012
    dt_all = []
    f_all = []
    nu_all = []

    freq = 8.5E9
    choose = np.logical_and(~islim, nu==freq)
    dt_all.extend(dt[choose])
    f_all.extend(f[choose])
    nu_all.extend([freq]*sum(choose))

    freq = 6.7E9
    dt_all.extend([16.4, 20.4, 25.4, 35.3, 45.3, 58.2, 92.9])
    f_all.extend([4.09, 4.8, 5.98, 7.222, 6.987, 6.11, 3.941])
    nu_all.extend([freq]*7)

    dt_all = np.array(dt_all)
    f_all = np.array(f_all)
    nu_all = np.array(nu_all)

    lum = plot_line(
            ax, d, dt_all, nu_all*f_all, 
            'SN2011dh', 'SN', col, legend)
    ax.text(dt[choose][0]/1.05, lum[0], 'SN2011dh', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


def grb030329(ax, col, legend):
    """ 
    Berger 2003
    Van der Horst et al. 2008
    
    Explosion day was obviously 03/29
    """
    z = 0.1686
    d = Planck15.luminosity_distance(z=z).cgs.value

    # LOW FREQUENCY

    # Berger: this is the best frequency to pick from this paper
    t = np.array(
            [0.58, 1.05, 2.65, 3.57, 4.76, 6.89, 7.68, 9.49, 11.90, 
                12.69, 14.87, 16.66, 18.72, 20.58, 25.70, 28.44, 31.51, 
                33.58, 36.52, 42.55, 44.55, 59.55, 66.53]) / (1+z)
    f = np.array(
            [3.50, 1.98, 8.50, 6.11, 9.68, 15.56, 12.55, 13.58, 17.70, 
                17.28, 19.15, 17.77, 15.92, 16.08, 15.34, 12.67, 13.55, 
                13.10, 10.64, 8.04, 8.68, 4.48, 4.92])
    nu = np.array([8.5E9]*len(f))

    # Van der Horst: best frequency is 2.3 GHz
    t = np.append(t, np.array([268.577, 306.753, 365.524, 420.168, 462.078, 
        583.683, 743.892, 984.163]) / (1+z))
    f = np.append(
            f, np.array([1613, 1389, 871, 933, 707, 543, 504, 318]) * 1E-3)
    nu = np.append(nu, np.array([2.3E9]*8))
    lum = plot_line(ax, d, t, nu*f, 'GRB030329', 'GRB', col, legend)
    ax.text(t[6]*1.05, lum[10]*1.05, 'GRB030329', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='left')
    

def grb130427A(ax, col, legend):
    """ Perley et al
    They have data from CARMA/PdBI at 90 GHz (3mm)
    But by the time they caught it, it was fading
    """
    z = 0.340
    d = Planck15.luminosity_distance(z=z).cgs.value

    freq = 5.10E9
    t = np.array([0.677, 2.04, 4.75, 9.71, 17.95, 63.78, 128.34]) / (1+z)
    f = np.array([1290, 1760, 648, 454, 263, 151, 86]) * 1E-3

    freq = 6.8E9
    t = np.array([0.677, 2.04, 4.75, 9.71, 9.95, 12.92, 27.67, 59.8, 128]) / (1+z)
    f = np.array([2570, 1820, 607, 374, 385, 332, 243, 109, 91]) * 1E-3

    lum = plot_line(ax, d, t, freq*f, 'GRB130427A', 'GRB', col, legend)
    ax.text(t[-4], lum[-4]/1.3, 'GRB130427A', fontsize=11,
            verticalalignment='top',
            horizontalalignment='center')


def sn2007bg(ax, col, legend):
    """ Salas et al. 2013
    Peak is resolved for 4.86, 8.46 GHz
    """
    nu = 8.46E9
    d = Planck15.luminosity_distance(z=0.0346).cgs.value
    t = np.array(
            [13.8, 19.2, 26.1, 30.9, 41.3, 55.9, 66.8, 81.8, 98.8, 124, 
                144, 159.8, 189.9, 214.9, 250.9, 286.8, 314.8, 368.8, 
                386.8, 419.9, 566.9, 623.8, 720.8, 775.8, 863.8])
    f = np.array(
            [480, 753, 804, 728, 1257, 1490, 1390, 1325, 1131, 957, 
                621, 316, 379, 404, 783, 1669, 2097, 2200, 
                2852, 3344, 3897, 3891, 3842, 3641, 3408]) * 1E-3
    lum = plot_line(ax, d, t, nu*f, 'SN2007bg', 'Ic-BL SN', col, legend)
    ax.text(t[0]/1.05, lum[0], 'SN2007bg', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='right', zorder=0)


def ptf11qcj(ax, col, legend):
    """ Corsi et al. 2014
    """
    nu = 5E9
    t0 = 55842
    dat = ascii.read("../../data/ptf11qcj.txt")
    t = dat['col1']-t0
    nu_all = dat['col3']
    choose = nu_all == 5
    f = dat['col4']
    d = Planck15.luminosity_distance(z=0.0287).cgs.value

    lum = plot_line(ax, d, t[choose], nu*f[choose], 'PTF11qcj', 'Ic-BL SN', col, legend)
    ax.text(t[choose][15], lum[15]*1.5, 'PTF11qcj', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='center', zorder=0)


def sn2003bg(ax, col, legend):
    """ Soderberg et al.
    Peak is resolved for 22.5, 15, 8.46, 4.86, 1.43
    Again, there are two peaks...
    Let's choose the first peak, 8.46
    """
    nu = 8.46E9
    d = 6.056450393620008e+25

    t = np.array(
                [10, 12, 23, 35, 48, 58, 63, 73, 85, 91, 115, 129,
                132, 142, 157, 161, 181, 201, 214, 227, 242, 255,
                266, 285, 300, 326, 337, 351, 368, 405, 410, 424,
                434, 435, 493, 533, 632, 702, 756, 820, 902, 978])
    f = np.array(
                [2.51, 3.86, 12.19, 24.72, 40.34, 51.72, 49.64, 46.20,
                38.638, 33.85, 45.74, 53.94, 54.27, 54.83, 48.43,
                47.43, 35.76, 31.35, 28.67, 27.38, 24.57, 22.30,
                21.67, 21.31, 20.88, 20.33, 19.85, 18.84, 17.14,
                14.61, 14.49, 14.16, 13.25, 13.08, 10.04, 8.92,
                6.23, 6.18, 4.62, 3.93, 4.69, 4.48])
    lum = plot_line(ax, d, t, nu*f, 'SN2003bg', 'SN', col, legend, zorder=0)
    ax.text(t[-1]/1.05, lum[-1]/1.05, 'SN2003bg', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


def sn2009bb(ax, col, legend):
    """ expl date Mar 19 
    Soderberg 2010 """
    nu = 8.46E9
    d = 1.237517263280789e+26
    t_apr = 11 + np.array([5.2, 8.2, 13.2, 15.1, 23.2, 29.1])
    t_may = 11 + 30 + np.array([3.1, 10.1, 13, 20.1, 27])
    t_jun = 11 + 30 + 31 + np.array([6, 17, 26])
    t_jul = 11 + 30 + 31 + 30 + np.array([18.9])
    t_aug = 11 + 30 + 31 + 30 + 31 + np.array([11.8])
    t = np.hstack((t_apr, t_may, t_jun, t_jul, t_aug))
    flux = np.array([24.681, 17.568, 16.349, 13.812, 8.881,
        7.714, 8.482, 6.824, 6.327, 3.294, 4.204, 3.203, 2.392,
        1.903, 1.032, 1.084])
    lum = plot_line(ax, d, t, nu*flux, 'SN2009bb', 'Ic-BL SN', col, legend)
    ax.text(t[0]/1.05, lum[0], '2009bb', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


def sn2020bvc(ax, col, legend):
    nu = 10E9
    d = Planck15.luminosity_distance(z=0.02507).cgs.value
    t = np.array([13, 35, 77, 100, 163])
    flux = np.array([0.066, 0.051, 0.186, 0.195, 0.048])
    lum = nu*flux*1E-3*1E-23*4*np.pi*d**2
    print(lum)
    ax.errorbar(t,lum,fmt='*',c=col,ms=10)
    ax.plot(t,lum,c=col,lw=2)
    #lum = plot_line(ax, d, t, nu*flux, 'SN20bvc', 'Rel. SN', col, legend)
    ax.text(t[-1]*1.1, lum[-1]/1.05, r'SN20bvc', fontsize=11, 
            verticalalignment='top', horizontalalignment='left')


def sn1998bw(ax, col, legend):
    """ SN 1998bw
    from Kulkarni et al. 1998
    """
    d = 1.17E26 # cm
    nu = 150E9
    t = np.array([12.4])
    f = np.array([39])
    nu = 2.3E9
    t = np.array([11.7, 14.6, 15.7, 16.5, 17.8, 19.7, 21.6, 23.6, 25.9, 26.8, 28.8, 30.0, 32.9, 34.7, 36.8, 38.8, 40.0, 45.7, 51.7, 57.7, 64.7, 67.7, 80.5])
    f = np.array([19.7, 22.3, 23.5, 23.9, 25.1, 25.3, 20.9, 22.9, 28.0, 28.7, 31.1, 31.3, 27.3, 33.5, 31.8, 31, 31.3, 26.8, 23.1, 18.5, 15.6, 15.6, 9.6])
    nu = 10E9
    t = np.array([3, 4, 9.9, 11.7, 14.6, 15.7, 16.5, 17.8, 19.7, 21.6, 23.6, 
        25.9, 26.8, 28.8, 30.0, 32.9, 34.7, 36.8, 38.8, 40.0, 45.7, 51.7, 57.7,
        64.7, 67.7, 80.5])
    f = np.array([13, 13, 48, 49.4, 37.6, 34.3, 31.4, 26.2, 21.6, 18.8, 17.6, 
        20.9, 21.7, 22, 22.1, 21.3, 20.2, 18.4, 16.1, 16.6, 13.2, 10.5, 8.1, 
        7.7, 7.2, 4.1])
    lum = plot_line(ax, d, t, nu*f, 'SN1998bw', 'Rel. SN', col, legend)
    ax.text(t[0]/1.05, lum[0], '1998bw', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


def sn2006aj(ax, col, legend):
    """ SN 2006aj
    Soderberg 2006
    """
    d = Planck15.luminosity_distance(z=0.033).cgs.value
    nu = 8.46E9
    t = np.array([1.87,3,3.83,4.85,6.97,7.94,9.95,
        12.88,16.74,19.86,21.96,24.91,30.71,34.81,41.74,50.70,104.52])
    f = np.array([453,381,269,280,164,30,39,15,
        75,48,87,20,32,15,22,25,17])*1E-3
    lum = plot_line(ax, d, t, nu*f, 'SN1998bw', 'Rel. SN', col, legend)
    ax.text(t[0], lum[0], '2006aj', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='center')


def sn2010bh(ax, col, legend):
    """ SN 2010bh
    Margutti 2014
    """
    d = Planck15.luminosity_distance(z=0.0593).cgs.value
    nu = 5.4E9
    t = np.array([18.93, 29.87, 36.97, 69.87])
    f = np.array([90, 128, 68, 50])*1E-3
    lum = plot_line(ax, d, t, nu*f, 'SN2010bh', 'Rel. SN', col, legend)
    ax.text(t[0], lum[0], '2010bh', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='right')


def ztf18aaqjovh(ax, col, legend):
    """ Ho 2020 """
    d = Planck15.luminosity_distance(z=0.054).cgs.value
    nu = 6E9
    t = np.array([16, 22, 34])
    f = np.array([32.5, 29.6, 26.6])*1E-3
    lum = plot_line(
            ax, d, t, nu*f, 'ZTF18jovh', 'Ic-BL SN', col, legend, zorder=5)
    ax.text(t[0], lum[0], 'ZTF18jovh', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='right')


def sn2017iuk(ax, col, legend):
    """ SN 2017iuk
    from the GCN:
    https://gcn.gsfc.nasa.gov/gcn/gcn3/22216.gcn3
    """
    d = Planck15.luminosity_distance(z=0.037).cgs.value
    nu = 6E9
    t = np.array([4.3, 9, 18, 33, 50])
    f = np.array([3, 7.36, 10.1, 7.80, 7.89])
    lum = plot_line(ax, d, t, nu*f, 'SN1998bw', 'Rel. SN', col, legend)
    #lum = f*1E-3*1E-23*4*np.pi*d**2*nu
    #print(lum)
    #ax.scatter(t, lum, c=col, marker='s', s=6)
    ax.text(t[0], lum[0], 'SN2017iuk', fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='right')


def ztf20abzoeiw(ax, col, legend):
    """ Ho 2020 """
    d = Planck15.luminosity_distance(z=0.01).cgs.value
    nu = np.array([10E9, 10E9, 6E9])
    t = np.array([20,35,110])
    f = np.array([286,564,95])*1E-3
    #lum = plot_line(
    #        ax, d, t, nu*f, 'SN1998bw', 'Ic-BL SN', col, legend, zorder=5)
    lum = nu*f*1E-3*1E-23*4*np.pi*d**2
    ax.plot(t, lum, c=col)
    ax.scatter(t, lum, c=col, s=10)
    #print(lum)
    ax.text(t[-1]/1.1, lum[-1], 'ZTF20eiw', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


if __name__=="__main__":
    fig, ax = plt.subplots(1, 1, figsize=(6,6), sharex=True, sharey=True)
    props = dict(boxstyle='round', facecolor='white')

    # viridis color palette
    # 440154 dark purple
    # 21918c dark green/turquoise
    # 3b528b dark blue
    # 5ec962 light green
    # fde725 yellow

    # inferno color palette
    # 000004 black
    # 420a68 dark purple
    # 932667 kind of a purple/red color
    # dd513a dark orange
    # fca50a light orange
    # fcffa4 really light yellow
    # 57106e purple
    # bc3754 pink/red
    # f98e09 orange

    # OK so I want to use four colors. I think it's nice to have two light ones
    # and two dark ones.

    #sn2003L(ax, 'lightblue', legend=True) # Ib
    sn1979c(ax, 'lightblue', None) # II
    sn1993J(ax, 'lightblue', None)
    sn2011dh(ax, 'lightblue', None)
    sn2007bg(ax, 'darkblue', legend=True)
    sn2003bg(ax, 'lightblue', None)
    ptf11qcj(ax, 'darkblue', None)
    sn2009bb(ax, 'darkblue', None)
    ztf18aaqjovh(ax, 'darkblue', None)
    ztf20abzoeiw(ax, 'darkblue', None)

    grb030329(ax, '#f98e09', legend=True)
    grb130427A(ax, '#f98e09', None)

    sn1998bw(ax, '#bc3754', legend=True)
    sn2006aj(ax, '#bc3754', None)
    sn2010bh(ax, '#bc3754', None)
    sn2020bvc(ax, 'darkblue', None)
    sn2017iuk(ax, '#bc3754', None)

    #at2018cow(ax, 'k', legend=True)
    #koala(ax, 'k', None)

    # new object...SN2021aabp
    x = 10
    dcm = Planck15.luminosity_distance(z=0.063).cgs.value
    y = (3.6*3)*1E-6*1E-23*4*np.pi*dcm**2*10E9
    ax.scatter(x, y, marker='o', c='darkblue', s=50)
    ax.arrow(x, y, 0, -y/1.8, length_includes_head=True,
            head_width=1.5, head_length=y/5, color='darkblue',
            zorder=100)
    ax.text(x, y/3, 'SN2021aabp', color='k', ha='center', va='top')

    moddir = "../../../../projects_proto/IcBL/old/data"

    # Plot model light curves
    dcm = Planck15.luminosity_distance(z=0.065).cgs.value
    # dat = np.loadtxt("../../data/offaxis_model_lc_1.57.txt", delimiter=',')
    # dt = dat[:,0] * 365.25
    # mjy = dat[:,1] 
    # lum = mjy * 1E-3 * 1E-23 * 4 * np.pi * dcm**2 * 8.5E9
    # ax.plot(dt,lum,c='k',ls='--', lw=0.5)

    dat = np.loadtxt(moddir + "/offaxis_model_lc_0.8.txt", delimiter=',')
    dt = dat[:,0] * 365.25
    mjy = dat[:,1] 
    lum = mjy * 1E-3 * 1E-23 * 4 * np.pi * dcm**2 * 5E9
    #ax.plot(dt,lum,c='k',ls='--', lw=0.5)

    dat = np.loadtxt(moddir + "/offaxis_model_lc_0.4.txt", delimiter=',')
    dt = dat[:,0] * 365.25
    mjy= dat[:,1]
    lum = mjy * 1E-3 * 1E-23 * 4 * np.pi * dcm**2 * 4.9E9
    #ax.plot(dt,lum,c='k',ls='--', lw=0.5)

    dat = np.loadtxt(moddir + "/offaxis_model_lc_0.txt", delimiter=',')
    dt = dat[:,0] * 365.25
    mjy = dat[:,1] 
    lum = mjy * 1E-3 * 1E-23 * 4 * np.pi * dcm**2 * 4.9E9
    #ax.plot(dt,lum,c='k',ls='--', lw=0.5)

    dat = np.loadtxt(moddir + "/model_cocoon_bottom.txt", delimiter=',')
    dt = dat[:,0]
    mjy = dat[:,1] # at 40 Mpc maybe?
    dcm = Planck15.luminosity_distance(z=0.009).cgs.value
    lum_bottom = mjy * 1E-3 * 1E-23 * 4 * np.pi * (dcm)**2 * 8.5E9

    dat = np.loadtxt(moddir + "/model_cocoon_top.txt", delimiter=',')
    mjy = np.interp(dt, dat[:,0], dat[:,1]) # at 40 Mpc maybe?
    lum_top = mjy * 1E-3 * 1E-23 * 4 * np.pi * (dcm)**2 * 8.5E9
    #ax.fill_between(dt, y1= lum_bottom, y2=lum_top, color='lightgrey', zorder=0)

    ax.set_ylabel(
            r"Luminosity $\nu L_{\nu}$ [erg s$^{-1}$]", 
            fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlim(1, 1000) 
    ax.set_ylim(1E34, 1E42)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r"Time [days; rest frame]", fontsize=16)

    #ax.scatter(
    #        0,0,c='k',marker='*',s=100,label="Fast-Lum. Opt. Transient")
    ax.legend(fontsize=12, loc='upper right', ncol=2, columnspacing=1)

    #ax.axhspan(1E34,1E37,edgecolor='k', fc='white', lw=3)
    #ax.axhline(y=1E37, c='k', ls='--')
    #ax.axvspan(1,2000,edgecolor='k', fc='white', lw=3)
    #ax.text(
    #        1.1,9E36,"Ordinary SNe", fontstyle='italic', fontsize=12,
    #        verticalalignment='top')


    plt.tight_layout()
    #plt.show()
    plt.savefig("lum_evolution_with_sn2021aabp.png", format='png', dpi=300)
    plt.close()
