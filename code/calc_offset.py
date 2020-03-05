""" Calc the offset from the HII region """

from penquins import kowalski


# from the legacysurvey distance
s = logon()
q = {"query_type": "find",
     "query": {
        "catalog": "ZUDS_alerts",
        "filter": {'objectId': {'$eq': 'ZUDS20eekej'}},
        "projection": {"candidate.lsdistnr1": 1}}}
query_result = s.query(query=q)
out = query_result['result_data']['query_result']
dist = np.array([val['candidate']['lsdistnr1'] for val in out])
np.mean(dist), np.std(dist)


# distnr in the reference image
q = {"query_type": "find",
     "query": {
        "catalog": "ZTF_alerts",
        "filter": {'objectId': {'$eq': 'ZTF20aalxlis'}},
        "projection": {"candidate.magpsf": 1, "candidate.distnr": 1}}}
query_result = s.query(query=q)
out = query_result['result_data']['query_result']
dist = np.array([val['candidate']['distnr'] for val in out])
mag = np.array([val['candidate']['magpsf'] for val in out])
np.mean(dist), np.std(dist)


