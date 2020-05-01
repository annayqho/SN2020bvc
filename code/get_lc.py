""" Return the extinction-subtracted photometry for SN2020bvc """

import numpy as np
import extinction
from astropy.io import ascii

# a mapping from filter to central wavelength
bands = {}
bands['V'] = 5468
bands['B'] = 4392
bands['U'] = 3465
bands['UVW1'] = 2600
bands['UVM2'] = 2246
bands['UVW2'] = 1928
# p48
bands['g'] = 4722.7
bands['r'] = 6339.6
bands['i'] = 7886.1
# lt: assume everything is the same except
bands['u'] = 3513.7
bands['z'] = 8972.9

# extinction for each filter
ext = {}
for band in bands.keys():
    ext[band] = extinction.fitzpatrick99(
        np.array([bands[band]]), 0.034, 3.1)[0]


def get_opt_lc():
    dat = ascii.read("/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/marshal_lc.txt")
    t = dat['jdobs']
    filt = dat['filter']
    mag = dat['magpsf']
    emag = dat['sigmamagpsf']
    inst = dat['instrument']
    maglim = dat['limmag']

    mag_corr = np.zeros(len(mag))
    for ii,f in enumerate(np.unique(filt)):
        choose = filt == f
        mag_corr[choose] = mag[choose]-ext[f]

    return t,mag_corr,emag,maglim,filt,inst


def get_uv_lc():
    # How do I correct these for Galactic extinction?
    uvdat = ascii.read("../../data/UVOT_hostsub.ascii")
    uvt = uvdat['MJD']+2400000.5
    uvdt = uvt-t0
    uvfilt = uvdat['FILTER']
    uvflux = uvdat['AB_FNU_mJy']
    uveflux = uvdat['AB_FNU_mJy_ERRM']
