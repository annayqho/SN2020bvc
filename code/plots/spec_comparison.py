""" Comparison of SN2020bvc spectra with other Ic-BL SNe
Combining three panels into one """

import matplotlib.pyplot as plt
from day1_bbfit import *
from day4_spec_comparison import *
from peak_spec_comparison import *

fig,axarr = plt.subplots(3,1,figsize=(7,10), sharex=True,
        gridspec_kw = {'height_ratios': [2,1.5,2]})

ax = axarr[0]
panel(ax)
ax.legend(fontsize=12)
ax.set_yscale('log')
ax.tick_params(axis='both', labelsize=16)
ax.set_xticks([])
ax.set_ylim(8E-2, 3)
ax.text(0.5,0.9,'$\Delta t$=0.7d', fontsize=16, transform=ax.transAxes)

ax = axarr[1]
plot_day4(ax)
ax.set_yscale('log')
ax.tick_params(axis='both', labelsize=16)
ax.set_ylabel(
    r"Flux ($10^{-15}$ erg\,s$^{-1}$\,cm$^{-2}\,$\AA$^{-1}$)", fontsize=16)
ax.set_xticks([])
ax.set_ylim(1E-2, 2)
ax.text(0.5,0.9,'$\Delta t$=3.7d', fontsize=16, transform=ax.transAxes)

ax = axarr[2]
plot_peak(ax)
ax.set_xlabel("Rest Wavelength (\AA)", fontsize=16)
ax.set_xlim(3300, 9300)
ax.set_xticks([4000,5000,6000,7000,8000,9000])
ax.tick_params(axis='both', labelsize=16)
ax.text(0.5,0.9,'$\Delta t$=12.5d', fontsize=16, transform=ax.transAxes)

fig.subplots_adjust(hspace=0)

#plt.tight_layout()
#plt.show()
plt.savefig("full_spec_comparison.png", dpi=300)
