""" calculate the nickel mass """

import numpy as np

eni = 3.9E10
eco = 6.8E9
beta = 9/8
tni = 8.8
tco = 111.3
L = (3.85)*1E42

tterms = 0.83*(1-beta*tni)*np.exp(-beta*tni) + 26.56 * (1-(1+beta*tco)*np.exp(-beta*tco))
tterms_s = tterms * 86400
mni = (L * beta**2 * (tni*86400)**2 / (2*eni)) / tterms
print(mni/2E33)
