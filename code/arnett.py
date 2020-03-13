"""  Calculate the nickel mass of the explosion using the formulation in 
Valenti+ 2008 """
import numpy as np
from scipy.integrate import quad


lph = 3.67E42
eni = 3.90E10
eco = 6.78E9
tni = 8.8
tco = 113.6
tm = 8
y = tm/(2*tni)
s = (tm*(tco-tni)/(2*tco*tni))
mni = 0.11*2E33

def a(z):
    return 2*z*np.exp(-2*z*y+z**2)

def b(z):
    return 2*z*np.exp(-2*z*y+2*z*s+z**2)

def lph(t):
    return mni*np.exp(-(t/tm)**2) * \
            ((eni-eco)*quad(a,0,t/tm)[0]+eco*quad(b,0,t/tm)[0])
