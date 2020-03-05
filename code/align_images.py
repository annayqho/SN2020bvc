""" Align images using SWARP

I ran this on gayatri,
but I wasn't able to get this to work :(
"""

import os

swarp_path = '/usr/local/optical/swarp/bin/swarp'

ddir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/host_sdss"

gfits = "../data/host_sdss/frame-g-003813-4-0393.fits"
rfits = "../data/host_sdss/frame-r-003813-4-0393.fits"
ifits = "../data/host_sdss/frame-i-003813-4-0393.fits"
zfits = "../data/host_sdss/frame-z-003813-4-0393.fits"
ufits = "../data/host_sdss/frame-u-003813-4-0393.fits"

# centering
ra = 218.487548
dec = 40.243758
size = 600

swarp_command = swarp_path \
        + " %s %s %s %s %s -c config.swarp -CENTER '" %(gfits,rfits,ifits,zfits,ufits)\
        + str(ra) + " " + str(dec) \
        + "' -SUBTRACT_BACK Y -RESAMPLE Y -COMBINE N -IMAGE_SIZE '" \
        + str(size) + "," + str(size) + "'"
os.system(swarp_command)
