""" Print table showing blackbody fits """
import numpy as np

dat = np.loadtxt("../plots/bol_lc.txt")
dt = dat[:,0]
n = len(dt)

for ii in np.arange(n):
    dtstr = str(np.round(dt[ii],1))
    scale = 1E42
    root = 1
    lstr = str(np.round(dat[:,root][ii]/scale,1))
    ulstr = str(np.round((dat[:,root+1][ii]-dat[:,root][ii])/scale,1))
    llstr = str(np.round((dat[:,root][ii]-dat[:,root+2][ii])/scale,1))
    scale = 1E3
    root = 4
    tstr = str(np.round(dat[:,root][ii]/scale,1))
    utstr = str(np.round((dat[:,root+1][ii]-dat[:,root][ii])/scale,1))
    ltstr = str(np.round((dat[:,root][ii]-dat[:,root+2][ii])/scale,1))
    scale = 1E14
    root = 7
    rstr = str(np.round(dat[:,root][ii]/scale,1))
    urstr = str(np.round((dat[:,root+1][ii]-dat[:,root][ii])/scale,1))
    lrstr = str(np.round((dat[:,root][ii]-dat[:,root+2][ii])/scale,1))
    print("%s & $%s^{+%s}_{-%s}$ & $%s^{+%s}_{-%s}$ & $%s^{+%s}_{-%s}$ \\\\" %(
        dtstr,lstr,ulstr,llstr,tstr,utstr,ltstr,rstr,urstr,lrstr))
