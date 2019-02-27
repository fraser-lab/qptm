from __future__ import division
from ptm_util import PTM_lookup

#I 15 G N2 m2G (N2-Methylguanosine) 20.5871583956 7.36784997944 0.715771437502

expected = {
  "I":{
    745:"m1G",
    # 746:"PSU",
    747:"m5U",
    # 955:"PSU",
    1618:"m6A",
    1835:"m2G",
    # 1911:"PSU",
    1915:"m3U", # actually 3-methylpseudouridine
    # 1917:"PSU",
    1939:"m5U",
    1962:"m5C",
    2030:"m6A",
    2069:"m7G",
    2251:"m(2'O)G",
    2445:"m2G",
    # 2449:"DHU", # 5,6-dihydrouridine
    # 2457:"PSU",
    2498:"m(2'O)C",
    2503:"m2A",
    # 2504:"PSU",
    2552:"m(2'O)U",
    # 2580:"PSU",
    # 2604:"PSU",
    # 2605:"PSU"
  },
  "J":{
    # 516:"PSU",
    527:"m7G",
    966:"m2m2G",
    967:"m5C",
    1207:"m2G",
    1402:"m4C m(2'O)C",
    1407:"m5C",
    1498:"m3U",
    1516:"m2G",
    1518:"m6m6A",
    1519:"m6m6A"
  }
}

if __name__ == "__main__":

  ptms_found = []

  with open("ptms.out") as ptmsfile:
    for line in ptmsfile:
      try:
        ch, resid, resn, loc, ptm, ptm_full, _, _, _ = line.split()
        resid = int(resid)
        if resid in expected[ch].keys() and ptm in expected[ch][resid].split():
          print "+++  true positive: ", " ".join([ch, str(resid), ptm])
          print line
          ptms_found.append(resid)
        # else:
        #   print "--- false positive: ", " ".join([ch, str(resid), ptm])
      except Exception:
        pass

  for resid, ptm in expected["I"].iteritems():
    if not resid in ptms_found:
      print "---  false negative: ", " ".join(["I", str(resid), ptm]), "\n"
