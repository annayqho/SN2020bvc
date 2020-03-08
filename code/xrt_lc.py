""" Get the XRT LC info """

import numpy as np
from astropy.time import Time

f = "/Users/annaho/Dropbox/Projects/Research/SN2020bvc/data/xray/swift_lc.txt"
dat = np.loadtxt(f)
dt_swift = dat[:,0]
t0_swift = Time("2013-05-04T14:30:32.484", format='isot')
t_swift = t0_swift + dt_swift/86400
print(t_swift.mjd)

t0_sn = Time(2458882.0568, format='jd') # last non-detection
dt_sn = t_swift-t0_sn
print(dt_sn.value)
