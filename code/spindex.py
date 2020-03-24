""" Use the radio data to estimate the spectral index """

import numpy as np

nu = np.array([6, 10])
fnu = np.array([84, 66])
efnu = np.array([5, 5])

nsim = 1000
alphas = np.zeros(nsim)
ysamples = np.zeros((nsim, len(nu)))

for jj,val in enumerate(fnu):
    ysamples[:,jj] = np.random.normal(loc=val,scale=efnu[jj],size=nsim)

for jj in np.arange(nsim):
    alphas[jj] = np.log10(
        ysamples[jj,0]/ysamples[jj,1])/np.log10(nu[0]/nu[1])

np.mean(alphas) # = -0.47
np.std(alphas) # 0.18

# Now do it again for epochs around 30 days

nu = np.array([10,15])
fnu = np.array([51,33])
efnu = np.array([5,4])

nsim = 1000
alphas = np.zeros(nsim)
ysamples = np.zeros((nsim, len(nu)))

for jj,val in enumerate(fnu):
    ysamples[:,jj] = np.random.normal(loc=val,scale=efnu[jj],size=nsim)

for jj in np.arange(nsim):
    alphas[jj] = np.log10(
        ysamples[jj,0]/ysamples[jj,1])/np.log10(nu[0]/nu[1])

np.mean(alphas) # = -1.05
np.std(alphas) # 0.39


