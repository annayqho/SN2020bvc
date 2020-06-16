"""  Calculate the nickel mass of the explosion using the formulation in 
Valenti+ 2008 """
import numpy as np
from scipy.integrate import quad

# CONSTANTS
eni = 3.90E10
eco = 6.78E9
tni = 8.8
tco = 113.6
tm = 13
y = tm/(2*tni)
s = (tm*(tco-tni)/(2*tco*tni))
mni = 0.16*2E33
B = 9/8


def a(z):
    return 2*z*np.exp(-2*z*y + z**2)

def b(z):
    return 2*z*np.exp(-2*z*y+2*z*s+z**2)

def lph(t):
    return mni*np.exp(-(t/tm)**2) * \
            ((eni-eco)*quad(a,0,t/tm)[0]+eco*quad(b,0,t/tm)[0])

def c(z):
    return z*((eni-eco)*np.exp(-z/tni)+eco*np.exp(-z/tco))*np.exp((z/tm)**2)


def lph_khatami(t):
    return (2*mni/tm**2)*np.exp(-(t/tm)**2) * \
            quad(c,0,t)[0]


def mej():
    kappa = 0.1
    c = 3E10
    beta = 13.8
    Lam = 6/5
    v = 0.07*3E10
    mej = np.sqrt((3*(8*86400)**4 * v**2 / (10*Lam)) * (beta*c/kappa)**2) 
    print(mej/2E33)
    ek = (3/10)*mej*v**2
    print(ek)


if __name__=="__main__":
    mej()
