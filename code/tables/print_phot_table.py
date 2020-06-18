""" Print table of photometry, ground-based and UVOT """

import numpy as np
from math import floor, log10
from astropy.time import Time
from astropy.cosmology import Planck15
from astropy.io import ascii
import sys
import extinction
sys.path.append("/Users/annaho/Dropbox/Projects/Research/SN2020bvc/code")
from get_lc import get_opt_lc,get_uv_lc

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

def round_sig(x, sig=2):
    print(x)
    if x < 0:
        return -round(-x, sig-int(floor(log10(-x)))-1)
    return round(x, sig-int(floor(log10(x)))-1)


def ndec(num):
    dec = str(num).split('.')[-1]
    return len(dec)


d = Planck15.luminosity_distance(z=0.025201).cgs.value

headings = np.array(
        ['Date', '$\Delta t$', 'Inst.', 'Filt.', 
         'Mag'])
units = np.array(
        ['(MJD)', '(d)', '', '', '(AB)'])
label = "uvot-phot"
caption = "Optical and ultraviolet photometry for SN\,2020bvc,\
corrected for Milky Way extinction."

# Print the table headers
ncol = len(headings)
colstr = ""
colstr += 'l'
for col in np.arange(ncol-1): colstr+="r"
print(colstr)

colheadstr = ""
for col in np.arange(ncol-1):
    colheadstr += "\colhead{%s} & " %headings[col]
colheadstr += "\colhead{%s}" %headings[-1]

unitstr = ""
for col in np.arange(ncol-1):
    unitstr += "\colhead{%s} & " %units[col]
unitstr += "\colhead{%s}" %units[-1]

rowstr = ""
for col in np.arange(ncol-1):
    rowstr += "%s & "
rowstr += "%s \\\ \n"

outputf = open("table_%s.txt" %label, "w")
outputf.write("\\startlongtable \n")
outputf.write("\\begin{deluxetable}{%s} \n" %colstr)
outputf.write("\\tablecaption{%s\label{tab:%s}} \n" %(caption,label))
outputf.write("\\tablewidth{0pt} \n")
outputf.write("\\tablehead{ %s \\\ %s} \n" %(colheadstr, unitstr))
#outputf.write("\\rotate \n")
outputf.write("\\tabletypesize{\scriptsize} \n")
outputf.write("\startdata \n")

# time of the last ATLAS non-detection
t0 = 2458883.17

# UV 
uvjd,uvdt,uvfilt,uvflux,uveflux= get_uv_lc()
uvmjd = uvjd-2400000.5
choose = uvdt > 0
uvmjd = uvmjd[choose]
uvdt = uvdt[choose]
uvfilt = uvfilt[choose]
uvflux = uvflux[choose]
uveflux = uveflux[choose]

# optical
t,mag,emag,maglim,filt,instr = get_opt_lc()
# only detections
choose = emag < 90
t = t[choose]
mag = mag[choose]
emag = emag[choose]
maglim = maglim[choose]
filt = filt[choose]
instr = instr[choose]
mjd = Time(t, format='jd').mjd
dt = t-t0

# Add the UVOT data
mjd = np.append(mjd, uvmjd)
dt = np.append(dt, uvdt)
filt = np.append(filt, uvfilt)
mag = np.append(mag, -2.5*np.log10(uvflux*1E-3/3631))
emag = np.append(emag, uveflux/uvflux)
instr = np.append(instr, np.array(['Swift+UVOT']*len(uvdt)))

# Now, reorder according to dt
order = np.argsort(dt)
mjd = mjd[order]
dt = dt[order]
filt = filt[order]
mag = mag[order]
emag = emag[order]
instr = instr[order]

for ii in np.arange(len(dt)):
    # Convert the flux into a fluxstr
    if mag[ii] < 99.0:
        # If not an upper limit, print row
        mjd_str = round_sig(mjd[ii], 11)
        dt_str = np.round(dt[ii], 2)
        mag_corr = mag[ii]-ext[filt[ii]]
        mag_str = str('{:.2f}'.format(round_sig(mag_corr, 4)))
        emag_str = "%.2f" %np.round(emag[ii],2)
        if emag_str == "0.00":
            emag_str = "0.01"
        print(emag_str)
        row = rowstr %(
                mjd_str, dt_str, instr[ii], "$%s$" %filt[ii], 
                "$%s \pm %s$"%(mag_str, emag_str))
        outputf.write(row)

outputf.write("\enddata \n")
outputf.write("\end{deluxetable} \n")
