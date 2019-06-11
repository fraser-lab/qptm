from __future__ import division
from ptm_util import PTM_lookup

#I 15 G N2 m2G (N2-Methylguanosine) 20.5871583956 7.36784997944 0.715771437502

expected = {
  "I":{
    745:"m1G",
    ## 746:"PSU",
    747:"m5U",
    ## 955:"PSU",
    1618:"m6A",
    1835:"m2G",
    ## 1911:"PSU",
    1915:"m3U", # actually 3-methylpseudouridine
    ## 1917:"PSU",
    1939:"m5U",
    1962:"m5C",
    2030:"m6A",
    2069:"m7G",
    2251:"m(2'O)G",
    2445:"m2G",
    ## 2449:"DHU", # 5,6-dihydrouridine
    ## 2457:"PSU",
    2498:"m(2'O)C",
    2503:"m2A",
    ## 2504:"PSU",
    2552:"m(2'O)U",
    ## 2580:"PSU",
    ## 2604:"PSU",
    ## 2605:"PSU"
  },
  "J":{
    ## 516:"PSU",
    # 527:"m7G",
    # 966:"m2m2G",
    # 967:"m5C",
    # 1207:"m2G",
    # 1402:"m4C m(2'O)C",
    # 1407:"m5C",
    # 1498:"m3U",
    # 1516:"m2G",
    # 1518:"m6m6A",
    # 1519:"m6m6A"
  }
}

def compare(outfile):

  ptms_found = []
  true_positives = 0
  false_positives = 0
  false_negatives = 0

  print "format: chain, resid, ptm_code, ref_atom, ptm_abbr, ptm_name, d_ref, d_new_diff, score\n"

  with open(outfile) as ptmsfile:
    for line in ptmsfile:
      try:
        ch, resid, resn, loc, ptm, ptm_full, _, _, _ = line.split()
        resid = int(resid)
        if resid in expected[ch].keys() and ptm in expected[ch][resid].split():
          print "+++  true positive: ", " ".join([ch, str(resid), ptm])
          print line
          ptms_found.append(resid)
          true_positives += 1
        else:
          # print "--- false positive: ", " ".join([ch, str(resid), ptm])
          false_positives += 1
      except Exception:
        pass

  for ch in expected.keys():
    for resid, ptm in expected[ch].iteritems():
      if not resid in ptms_found:
        print "---  false negative: ", " ".join([ch, str(resid), ptm]), "\n"
        false_negatives += 1

  print "true_positives: %d\nfalse_negatives: %d\nfalse_positives: %d\n" %\
    (true_positives, false_negatives, false_positives)

if __name__ == "__main__":
  compare("ptms.out")
