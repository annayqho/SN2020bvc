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

ddir = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/host_sdss"

# images have been resampled using SWARP,
# because they were previously out of alignment
gim = fits.open(ddir + "/frame-g-003813-4-0393.resamp.fits")
g = gim[0].data # 0 to 0.1 looks good
r = fits.open( # 0 to 0.2 looks good
        ddir + "/frame-r-003813-4-0393.resamp.fits")[0].data
z = fits.open( # 0 to 0.6 looks good
        ddir + "/frame-z-003813-4-0393.resamp.fits")[0].data
u = fits.open( 
        ddir + "/frame-u-003813-4-0393.resamp.fits")[0].data
i = fits.open( 
        ddir + "/frame-i-003813-4-0393.resamp.fits")[0].data

# Figure out pos from header
head = gim[0].header
wcs = astropy.wcs.WCS(head)
print(3600*astropy.wcs.utils.proj_plane_pixel_scales(wcs))
target_pix = wcs.all_world2pix([(np.array([ra,dec], np.float_))], 1)[0]
xpos = target_pix[0]
ypos = target_pix[1]
 
imsize = 50
# Plot cutout
gcut = g[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]
rcut = r[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]
zcut = z[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]
ucut = u[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]
icut = i[int(ypos-imsize):int(ypos+imsize),int(xpos-imsize):int(xpos+imsize)]

# icut: vmin = 0, vmax = 0.5
# rcut: vmin = 0, vmax = 0.3
# gcut: vmin = 0, vmax = 0.2

rgb = make_lupton_rgb(
        icut/2, rcut/1.5, gcut, Q=2, stretch=0.25)

fig,ax = plt.subplots()

ax.imshow(rgb, origin='lower')
ax.scatter(imsize, imsize, c='white', marker='x', s=50)
ax.plot((2*imsize-10,2*imsize-10), (2*imsize-10,2*imsize-20), color='white', lw=2)
ax.text(
        2*imsize-10, 2*imsize-23, "S", color='white', fontsize=16,
        horizontalalignment='center', verticalalignment='top')
ax.plot((2*imsize-10,2*imsize-20), (2*imsize-10,2*imsize-10), color='white', lw=2)
ax.text(
        2*imsize-23, 2*imsize-10, "E", color='white', fontsize=16,
        horizontalalignment='right', verticalalignment='center')
ax.axis('off')
# I think that the pixel scale is 0.3 arcsec
x = 10
y = 10
x2 = x + 5/0.3
ax.plot((x,x2), (y,y), c='white', lw=2)
ax.text((x2+x)/2, y/1.1, "5''", color='white', fontsize=16, 
        verticalalignment='top', horizontalalignment='center')

#plt.show()
plt.savefig(
"host.png", dpi=500, bbox_inches = 'tight', pad_inches = 0)

