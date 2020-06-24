""" Plot the light curve of this SN and compare it to SN 2006aj """

import numpy as np
from astropy.io import ascii
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.cosmology import Planck15
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/SN2020bvc/code")
from get_lc import get_opt_lc

rcol = 'Crimson'
gcol = 'Aquamarine'
icol = 'Goldenrod'
z = 0.025201
t0 = 2458883.17
msize = 6

""" a panel showing the full r, g, and i light curve from the P48 """
t,mag,emag,maglim,filt,inst = get_opt_lc()

fig,ax = plt.subplots(1,1, figsize=(10,5))
dt = (t-t0)/1.0252

# r-band light curve
choose = np.logical_and.reduce(
        (inst=='P48+ZTF', filt=='r', mag < 99, emag < 99))
ax.errorbar(dt[choose], mag[choose], emag[choose], ms=msize, fmt='o', 
        mfc=rcol, mec=rcol, label='P48 $r$', c=rcol, zorder=1)

# r-band upper limits
choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='r', emag >= 99, dt<5))
ax.scatter(
        dt[choose], maglim[choose], marker='o', 
        c=rcol, label=None, zorder=1)
for ii,val in enumerate(dt[choose]):
    ax.arrow(
            val, maglim[choose][ii], 0, 0.25, 
            color=rcol, head_length=0.2, head_width=0.3, zorder=1)

# g-band light curve
choose = np.logical_and.reduce(
        (inst=='P48+ZTF', filt=='g', mag < 99, emag < 99))
ax.errorbar(dt[choose], mag[choose], emag[choose], ms=msize, fmt='s', 
        mfc=gcol, mec=gcol, label='P48 $g$', c=gcol, zorder=2)

# add the ASAS-SN point
#ax.scatter(2458884.12-t0, 17, s=msize, marker='s',
#        facecolor='white', edgecolor=gcol, label='ASAS-SN $g$', zorder=2)

# g-band upper limits
choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='g', emag >= 99, dt<5))
ax.scatter(
        dt[choose], maglim[choose], marker='s', c=gcol, 
        label=None, zorder=2)
for ii,val in enumerate(dt[choose]):
    ax.arrow(
            val, maglim[choose][ii], 0, 0.25, 
            color=gcol, head_length=0.2, head_width=0.3, zorder=2)


# upper limits from ATLAS
val = 2458883.17-t0
ax.scatter(val, 19.4, marker='.', c='k')
ax.arrow(
        val, 19.4, 0, 0.25, 
        color='k', head_length=0.2, head_width=0.3, zorder=2)
ax.text(-0.2, 19.4, 'ATLAS $o$', fontsize=10,
        verticalalignment='center', horizontalalignment='right')

choose = np.logical_and.reduce(
        (inst=='P48+ZTF', filt=='i', mag < 99, emag < 99))
ax.errorbar(dt[choose], mag[choose], emag[choose], ms=msize, fmt='D', 
        mfc='white', mec=icol, label='P48 $i$', c=icol, zorder=0)
ax.text(val, 20, 'GRB 060218', rotation=270, fontsize=10)

# Spectral epochs
sps = [2458883.925137, 2458886.8544916, 58887.244648+2400000.5, 2458888.862627,
      Time("2020-02-12T12:16:26.713", format='isot').jd, 2458892.8246301,
      2458894.8325712, Time("2020-02-16T03:18:40.973", format='isot').jd,
      2458900.9285102, 2458908.9163794, 
      Time("2020-03-02T03:19:08.621", format='isot').jd]
for sp in sps:
    ax.text(sp-t0, 21.15, 'S', fontsize=10)


# Blackbody epochs
bb = np.array(
            [0.9, 1.36, 1.8, 2.8, 3.8, 4.74, 5.78, 6.27, 7.8, 9.1, 9.8, 10.75,
             11.09, 11.77, 12.47, 15.49, 20, 21.75, 23.77, 25.65, 26.5, 28.73,
             29.48])
for b in bb:
    ax.text(b, 15.9, 'B', fontsize=10)

# compare to 2006aj
offset = 0
dm = Planck15.distmod(z=0.0335).value-Planck15.distmod(z=0.025235).value
dat = ascii.read("../../data/lc_060218_full.txt")
t = dat['time']
mag = dat['magnitude']
emag = dat['e_magnitude']
band = dat['band']


# from the open SN catalog
dat = ascii.read("../../data/lc_060218_full.txt", delimiter=',')
b_plot = dat['band']
dt_plot = dat['time']-53784.149
t_plot = dat['time']
m_plot = dat['magnitude']
choose = np.logical_and(b_plot == 'B', dat['upperlimit']=='F')
# B-band extinction: 0.527
plt.scatter(
        dt_plot[choose][::1]/1.033, m_plot[choose][::1]-dm-0.527, 
        color=gcol, s=1, zorder=0, label=None)
plt.plot(
        dt_plot[choose][::1]/1.033, m_plot[choose][::1]-dm-0.527, 
        linestyle='-', lw=0.5, color=gcol, zorder=0, alpha=1, label="2006aj $B$")


choose = np.logical_and(b_plot == 'V', dat['upperlimit']=='F')
plt.scatter(
        dt_plot[choose][::1]/1.033, m_plot[choose][::1]-dm-0.399, 
        c=rcol, s=1, zorder=0, label=None)
plt.plot(
        dt_plot[choose][::1]/1.033, m_plot[choose][::1]-dm-0.399, 
        linestyle='--', lw=0.5, color=rcol, zorder=0, alpha=1, label="2006aj $V$")

# epoch of 060218
plt.axvline(x=0, c='k', lw=0.5, ls=':', label=None)

ax2 = ax.twinx()
ax2.set_ylabel(
        "Absolute Mag",
        fontsize=14, rotation=270, labelpad=15.0)
y_f = lambda y_i: y_i - Planck15.distmod(z=z).value
ymin, ymax = ax.get_ylim()
ax2.set_ylim((y_f(ymin), y_f(ymax)))
ax2.plot([],[])
ax2.tick_params(axis='both', labelsize=14)
ax2.invert_yaxis()

ax.legend(loc=(0.82, 0.06), fontsize=14)
ax.invert_yaxis()
ax.set_ylabel(r"Apparent Mag ($z=0.025201$)", fontsize=14)
ax.set_xlabel("Rest-Frame Days Since First Light", fontsize=16)
ax.tick_params(labelsize=14)
ax.tick_params(labelsize=14)
ax.set_xlim(-2.5,31)
ax.set_ylim(21.2,15.7)
plt.tight_layout()
plt.savefig("lc.eps", dpi=300, bbox_inches='tight')
#plt.show()
