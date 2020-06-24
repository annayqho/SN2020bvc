""" Plot an image of the host galaxy from SDSS 

This is Figure 1 in our paper.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import fits
import astropy.wcs
from astropy import coordinates as coords
from astropy.visualization import make_lupton_rgb

ra = 218.487548 
dec = 40.243758

ddir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data"

im = fits.open(ddir + "/cutout_218.4877_40.2435.fits")
g = im[0].data[0]  # 0 to 0.005 looks fine
r = im[0].data[1] # 0 to 0.01 looks fine
z = im[0].data[2] # 0 to 0.03

imsize = g.shape[0]

rgb = make_lupton_rgb(
        z/3, r/2, g, Q=2, stretch=0.1)

fig,ax = plt.subplots()

ax.text(0.95, 0.1, "Legacy Survey $grz$", fontsize=24, transform=ax.transAxes,
        horizontalalignment='right', color='white')

ax.imshow(rgb, origin='lower')

# mark the transient location
ax.plot([imsize/2, imsize/2], [imsize/2, imsize/2-10], c='white')
ax.plot([imsize/2, imsize/2+10], [imsize/2, imsize/2], c='white')

ax.plot((imsize-10,imsize-10), (imsize-10,imsize-20), color='white', lw=2)
ax.text(
        imsize-10, imsize-23, "S", color='white', fontsize=16,
        horizontalalignment='center', verticalalignment='top')
ax.plot((imsize-10,imsize-20), (imsize-10,imsize-10), color='white', lw=2)
ax.text(
        imsize-23, imsize-10, "E", color='white', fontsize=16,
        horizontalalignment='right', verticalalignment='center')
ax.axis('off')
# I think that the pixel scale is 0.3 arcsec
x = 10
y = 20
x2 = x + 5/0.3
ax.plot((x,x2), (y,y), c='white', lw=2)
ax.text((x2+x)/2, y/1.1, "5''", color='white', fontsize=16, 
        verticalalignment='top', horizontalalignment='center')

#plt.show()
plt.savefig(
"host.eps", format='eps', dpi=300, bbox_inches = 'tight')

