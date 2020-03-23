""" Print table of photometry, ground-based and UVOT """

import numpy as np
from math import floor, log10
from astropy.time import Time
from astropy.cosmology import Planck15
import sys

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
        ['Date (JD)', '$\Delta t$', 'Instrument', 'Filter', 
         'AB Mag', 'Error in AB Mag'])
label = "uvot-phot"
caption = "Optical and ultraviolet photometry for SN\,2020bvc"

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

rowstr = ""
for col in np.arange(ncol-1):
    rowstr += "%s & "
rowstr += "%s \\\ \n"

outputf = open("table_%s.txt" %label, "w")
outputf.write("\\startlongtable \n")
outputf.write("\\begin{deluxetable}{%s} \n" %colstr)
outputf.write("\\tablecaption{%s\label{tab:%s}} \n" %(caption,label))
outputf.write("\\tablewidth{0pt} \n")
outputf.write("\\tablehead{ %s } \n" %colheadstr)
#outputf.write("\\rotate \n")
outputf.write("\\tabletypesize{\scriptsize} \n")
outputf.write("\startdata \n")

dat = ascii.read("../../data/marshal_lc.txt")
uvdat = ascii.read("../../data/UVOT_hostsub.ascii")
uvt = uvdat['MJD']+2400000.5
uvdt = uvt-t0
uvfilt = uvdat['FILTER']
uvflux = uvdat['AB_FNU_mJy']
uveflux = uvdat['AB_FNU_mJy_ERRM']

t = dat['jdobs']
dt = t-t0
mag = dat['magpsf']
emag = dat['sigmamagpsf']
instr = dat['instrument']
filt = dat['filter']

# mjd = np.append(mjd.value, Time(dt_uv + t0.value, format='mjd').value)
# dt = np.append(dt, dt_uv)
# filt = np.append(filt, filt_uv)
# mag = np.append(mag, -2.5*np.log10(fnu_mjy_uv*1E-3/3631))
# emag = np.append(emag, efnu_mjy_uv/fnu_mjy_uv)
# tel = np.append(tel, np.array(['UVOT']*len(dt_uv)))

for ii in np.arange(len(dt)):
    # Convert the flux into a fluxstr
    if mag[ii] < 99.0:
        # If not an upper limit, print row
        mjd_str = round_sig(mjd[ii], 11)
        dt_str = np.round(dt[ii], 2)
        mag_str = str('{:.2f}'.format(round_sig(mag[ii], 4)))
        emag_str = "%.2f" %np.round(emag[ii],2)
        if emag_str == "0.00":
            emag_str = "0.01"
        print(emag_str)
        row = rowstr %(
                mjd_str, dt_str, tel[ii], filt[ii], 
                mag_str, emag_str)
        outputf.write(row)

outputf.write("\enddata \n")
outputf.write("\end{deluxetable} \n")
