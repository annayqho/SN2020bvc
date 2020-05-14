import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import numpy as np

fig,ax = plt.subplots(1,1,figsize=(7,4))

ax.errorbar(
        [10, 6, 3], [63, 83, 111], [6, 6, 10], 
        fmt='s', c='k', label="Feb")
ax.errorbar(
        [15, 10, 3, 6], [33, 50, 106, 60], [4, 5, 10, 6], 
        fmt='D', mec='k', mfc='white', label="March", c='k')
ax.errorbar(
        [3, 6, 10], [282, 209, 186], [10, 6, 7], 
        fmt='o', mec='k', mfc='red', label="April", c='k')
ax.errorbar(
        [3, 6, 10], [356, 258, 195], [10, 6, 7], 
        fmt='s', mec='k', mfc='blue', label="May", c='k')

xplot = np.linspace(3,15)
yplot = 113*(xplot/3)**(-0.75)
ax.plot(xplot, yplot, ls='--', c='k')
ax.text(3, 70, r"$F_\nu \propto \nu^{-0.75}$", fontsize=14)

ax.set_xlabel("Frequency (GHz)", fontsize=16)
ax.set_ylabel("Flux Density ($\mu$Jy)", fontsize=16)

ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=12)
ax.tick_params(axis='both', labelsize=14)
plt.tight_layout()

plt.show()
#plt.savefig("radio_sed.png", dpi=300)
