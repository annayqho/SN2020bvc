""" Plot the light curve of this SN and compare it to SN 2006aj """

import numpy as np
from astropy.io import ascii
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.cosmology import Planck15

EXT_G = 0.041
EXT_R = 0.029
EXT_I = 0.021
rcol = '#e55c30'
gcol = '#140b34'
icol = 'k'
z = 0.02520
t0 = 2458882.0568
msize = 6

""" a panel showing the full r, g, and i light curve from the P48 """
fig,ax = plt.subplots(1,1, figsize=(10,4))
dat = ascii.read("../../data/marshal_lc.txt")
dt = dat['jdobs']-t0
filt = dat['filter']
mag = dat['magpsf']
emag = dat['sigmamagpsf']
inst = dat['instrument']
maglim = dat['limmag']

# r-band light curve
choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='r', mag < 99))
ax.errorbar(dt[choose], mag[choose]-EXT_R, emag[choose], ms=msize, fmt='o', 
        mfc=rcol, mec=rcol, label='P48 $r$', c=rcol, zorder=1)

# r-band upper limits
choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='r', mag == 99, dt<5))
ax.scatter(
        dt[choose], maglim[choose]-EXT_R, marker='o', 
        c=rcol, label=None, zorder=1)
for ii,val in enumerate(dt[choose]):
    ax.arrow(
            val, maglim[choose][ii]-EXT_R, 0, 0.25, 
            color=rcol, head_length=0.2, head_width=0.3, zorder=1)

# g-band light curve
choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='g', mag < 99))
ax.errorbar(dt[choose], mag[choose]-EXT_G, emag[choose], ms=msize, fmt='s', 
        mfc=gcol, mec=gcol, label='P48 $g$', c=gcol, zorder=2)

# add the ASAS-SN point
#ax.scatter(2458884.12-t0, 17-EXT_G, s=msize, marker='s',
#        facecolor='white', edgecolor=gcol, label='ASAS-SN $g$', zorder=2)

# g-band upper limits
choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='g', mag == 99, dt<5))
ax.scatter(
        dt[choose], maglim[choose]-EXT_G, marker='s', c=gcol, 
        label=None, zorder=2)
for ii,val in enumerate(dt[choose]):
    ax.arrow(
            val, maglim[choose][ii]-EXT_G, 0, 0.25, 
            color=gcol, head_length=0.2, head_width=0.3, zorder=2)


# upper limits from ATLAS
val = 2458883.17-t0
ax.scatter(val, 19.4, marker='.', c='k')
ax.arrow(
        val, 19.4, 0, 0.25, 
        color='k', head_length=0.2, head_width=0.3, zorder=2)
ax.text(val, 19.4, 'ATLAS $o$', fontsize=10,
        verticalalignment='bottom', horizontalalignment='center')

choose = np.logical_and.reduce((inst=='P48+ZTF', filt=='i', mag < 99))
ax.errorbar(dt[choose], mag[choose]-EXT_G, emag[choose], ms=msize, fmt='D', 
        mfc='white', mec=icol, label='P48 $i$', c=icol, zorder=0)

# Spectral epochs
sps = [2458883.925137, 2458886.8544916, 58887.244648+2400000.5, 2458888.8626271,
      Time("2020-02-12T12:16:26.713", format='isot').jd, 2458892.8246301,
      2458894.8325712, Time("2020-02-16T03:18:40.973", format='isot').jd,
      2458900.9285102, 2458908.9163794, 
      Time("2020-03-02T03:19:08.621", format='isot').jd]
for sp in sps:
    ax.text(sp-t0, 21, 'S', fontsize=12)


# compare to 2006aj
dat = ascii.read("lc_060218_swift.dat")
t = dat['t']
mag = dat['mag']
emag = dat['emag']
plt.plot(t, mag-dm, linestyle='-', lw=1, color='k', label="2006aj (UVOT/B)")

ax2 = ax.twinx()
ax2.set_ylabel(
        "Absolute Mag",
        fontsize=16, rotation=270, labelpad=15.0)
y_f = lambda y_i: y_i - Planck15.distmod(z=z).value
ymin, ymax = ax.get_ylim()
ax2.set_ylim((y_f(ymin), y_f(ymax)))
ax2.plot([],[])
ax2.tick_params(axis='both', labelsize=14)
ax2.invert_yaxis()

ax.legend(loc='lower right', fontsize=14)
ax.invert_yaxis()
ax.set_ylabel(r"Apparent Mag", fontsize=14)
ax.set_xlabel("Days since last ZTF non-detection", fontsize=16)
ax.tick_params(labelsize=14)
ax.tick_params(labelsize=14)
ax.set_xlim(-1.5,30)
ax.set_ylim(21.2,16)
plt.tight_layout()
#plt.savefig("lc.png", dpi=200)
plt.show()
