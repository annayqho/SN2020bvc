""" Plot the bolometric light curve 
and show the best-fit explosion parameters """
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)

dt = np.array([1.36, 2.88, 3.81, 4.61, 6.27, 9.09, 10.86, 15.49, 26.51, 29.48])
lbol = np.array([4.97, 2.26, 1.98, 2.62, 2.95, 3.50, 3.67, 2.85, 2.19, 1.85])
elbol = np.array([0.15, 0.12, 0.07, 0.15, 0.15, 0.13, 0.09, 0.08, 0.06, 0.05])

plt.errorbar(dt/1.02507, lbol, yerr=elbol, c='k', fmt='o')
plt.xlabel("Rest-frame days since last non-detection", fontsize=14)
plt.ylabel("Bolometric luminosity ($10^{42}$ erg/s)", fontsize=14)
plt.tick_params(axis='both', labelsize=14)
plt.tight_layout()
plt.show()
