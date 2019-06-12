from __future__ import division
from scitbx.array_family import flex
from scitbx.python_utils import robust_statistics as rs

def import_ptms(ptms_outfile):
  chain_id = flex.std_string()
  resid = flex.int()
  resname = flex.std_string()
  goto_atom = flex.std_string()
  short_name = flex.std_string()
  full_name = flex.std_string()
  cc = flex.double()
  d_ref = flex.double()
  d_mid = flex.double()
  d_new_in_ref = flex.double()
  d_new_diff = flex.double()
  d_far = flex.double()
  ratio = flex.double()
  score = flex.double()
  imported_ptms = (chain_id, resid, resname, goto_atom, short_name, full_name, cc, 
    d_ref, d_mid, d_new_in_ref, d_new_diff, d_far, ratio, score)
  with open(ptms_outfile, "rb") as out:
    for line in out.readlines():
      items = line.split()
      for array, coerce_lambda in zip(imported_ptms,
        (lambda x: x, lambda x: int(x), lambda x: x, lambda x: x, lambda x: x, lambda x: x,
          lambda x: float(x), lambda x: float(x), lambda x: float(x), lambda x: float(x),
          lambda x: float(x), lambda x: float(x), lambda x: float(x), lambda x: float(x))):
        array.append(coerce_lambda(items.pop(0)))
  return imported_ptms

def apply_filters(imported_ptms, cc_threshold=0.7, ref_frac=1, dif_frac=1,
  ratio_d_far_d_new_in_ref=2, ratio_d_far_d_mid=1, ratio_d_ref_d_new_in_ref=3,
  score_threshold=0):
  # use import_ptms above if reading from file
  chain_id, resid, resname, goto_atom, short_name, full_name, cc, d_ref, d_mid, \
    d_new_in_ref, d_new_diff, d_far, ratio, score = imported_ptms
  # apply cc threshold
  cc_threshold_sel = cc >= cc_threshold
  keep_selection = cc_threshold_sel
  # apply ratio filters
  ratio_d_far_d_new_in_ref_sel = d_far <= ratio_d_far_d_new_in_ref * d_new_in_ref
  keep_selection = keep_selection & ratio_d_far_d_new_in_ref_sel
  ratio_d_far_d_mid_sel = d_far <= ratio_d_far_d_mid * d_mid
  keep_selection = keep_selection & ratio_d_far_d_mid_sel
  ratio_d_ref_d_new_in_ref_sel = d_ref <= ratio_d_ref_d_new_in_ref * d_new_in_ref
  keep_selection = keep_selection & ratio_d_ref_d_new_in_ref_sel
  # filter on reference density
  if ref_frac < 1:
    keep_density = rs.percentile(d_ref, 1-ref_frac)
    ref_frac_sel = d_ref >= keep_density
    keep_selection = keep_selection & ref_frac_sel
  else:
    ref_frac_sel = flex.bool(len(keep_selection), True)
  # filter on difference density
  if dif_frac < 1:
    keep_density = rs.percentile(d_new_diff, 1-dif_frac)
    dif_frac_sel = d_new_diff >= keep_density
    keep_selection = keep_selection & dif_frac_sel
  else:
    dif_frac_sel = flex.bool(len(keep_selection), True)
  # filter on score
  score_threshold_sel = score >= score_threshold
  keep_selection = keep_selection & score_threshold_sel
  # apply filters and update logs
  log = flex.std_string(len(keep_selection))
  for i in xrange(len(keep_selection)):
    if not cc_threshold_sel[i]:
      log[i] += "rejected for cc < %4.2f. " % cc_threshold
    if not ratio_d_far_d_new_in_ref_sel[i]:
      log[i] += "rejected for d_far > %4.2f * d_new_in_ref. " % ratio_d_far_d_new_in_ref
    if not ratio_d_far_d_mid_sel[i]:
      log[i] += "rejected for d_far > %4.2f * d_mid. " % ratio_d_far_d_mid
    if not ratio_d_ref_d_new_in_ref_sel[i]:
      log[i] += "rejected for d_ref > %4.2f * d_new_in_ref. " % ratio_d_ref_d_new_in_ref
    if not ref_frac_sel[i]:
      log[i] += "rejected for reference density not in best %4.2f. " % ref_frac
    if not dif_frac_sel[i]:
      log[i] += "rejected for difference density not in best %4.2f. " % dif_frac
    if not score_threshold_sel[i]:
      log[i] += "rejected for score < %4.2f. " % score_threshold
  # prepare results
  accepted = (
    chain_id.select(keep_selection),
    resid.select(keep_selection),
    resname.select(keep_selection),
    goto_atom.select(keep_selection),
    short_name.select(keep_selection),
    full_name.select(keep_selection),
    cc.select(keep_selection),
    d_ref.select(keep_selection),
    d_mid.select(keep_selection),
    d_new_in_ref.select(keep_selection),
    d_new_diff.select(keep_selection),
    d_far.select(keep_selection),
    ratio.select(keep_selection),
    score.select(keep_selection),
    log.select(keep_selection))
  all_tested = (
    chain_id,
    resid,
    resname,
    goto_atom,
    short_name,
    full_name,
    cc,
    d_ref,
    d_mid,
    d_new_in_ref,
    d_new_diff,
    d_far,
    ratio,
    score,
    log)
  return (keep_selection, accepted, all_tested, log)

def write_ptms_from_flex_arrays(tuple_of_arrays, outfile):
  chain_id, resid, resname, goto_atom, short_name, full_name, cc, d_ref, d_mid, \
    d_new_in_ref, d_new_diff, d_far, ratio, score, log = tuple_of_arrays
  with open(outfile, "wb") as out:
    for i in xrange(len(chain_id)):
      out.write(" ".join([chain_id[i], str(resid[i]), resname[i], goto_atom[i],
        short_name[i], full_name[i], str(cc[i]), str(d_ref[i]), str(d_mid[i]),
        str(d_new_in_ref[i]), str(d_new_diff[i]), str(d_far[i]), str(ratio[i]),
        str(score[i]), log[i]]) + "\n")

