""" Get the LC from ZUDS """

import penquins
from penquins import Kowalski
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF_fast_transient_search/code")
import ForcePhotZTF
from ForcePhotZTF.phot_class import ZTFphot
from run_forced_phot import get_forced_phot


def logon():
    """ Log onto Kowalski """
    username = 'ah'
    password = 'TetraodonInsists'
    s = Kowalski(
        protocol='https', host='kowalski.caltech.edu', port=443,
        verbose=False, username=username, password=password)
    return s

# lc doesn't look very good...not using this
def from_danny():
    s = logon()
    q = {"query_type": "find",
         "query": {
            "catalog": "ZUDS_alerts_aux",
            "filter": {'_id': {'$eq': 'ZUDS20eekej'}},
            "projection": {"light_curve": 1}}} 
    query_result = s.query(query=q)
    out = query_result['result_data']['query_result'][0]['light_curve']

    mjd = np.array([point['mjd'] for point in out])
    filt = np.array([point['filter'] for point in out])
    flux = np.array([point['flux'] for point in out])
    fluxerr = np.array([point['fluxerr'] for point in out])
    zp = np.array([point['zp'] for point in out])
    mag = -2.5*np.log10(flux)+zp


# regular forced photometry from Yuhan
# I'm honestly not sure about this either!
zp,filt,jd,flux,eflux,mag,emag = get_forced_phot('ZTF20aalxlis', 
        218.487548, 40.243758, 2458882.0568, [-2, 30])

# Maybe I'll just stick with the regular IPAC photometry for now...
# the alerts, that is.
