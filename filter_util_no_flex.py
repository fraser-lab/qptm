from __future__ import division

def import_ptms(ptms_outfile):
  chain_id = []
  resid = []
  resname = []
  goto_atom = []
  short_name = []
  full_name = []
  cc = []
  d_ref = []
  d_mid = []
  d_new_in_ref = []
  d_new_diff = []
  d_far = []
  ratio = []
  score = []
  model_id = []
  d_min = []
  b_factor = []
  imported_ptms = (chain_id, resid, resname, goto_atom, short_name, full_name, cc, 
    d_ref, d_mid, d_new_in_ref, d_new_diff, d_far, ratio, score, model_id, d_min, b_factor)
  with open(ptms_outfile, "r") as out:
    for line in out.readlines():
      items = line.split()
      for array, coerce_lambda in zip(imported_ptms,
        (lambda x: x, lambda x: int(x), lambda x: x, lambda x: x, lambda x: x, lambda x: x,
          lambda x: float(x), lambda x: float(x), lambda x: float(x), lambda x: float(x),
          lambda x: float(x), lambda x: float(x), lambda x: float(x), lambda x: float(x),
          lambda x: int(x), lambda x: float(x), lambda x: float(x))):
        array.append(coerce_lambda(items.pop(0)))
  return imported_ptms

def import_ptms_with_predictions(ptms_outfile):
  chain_id = []
  resid = []
  resname = []
  goto_atom = []
  short_name = []
  full_name = []
  cc = []
  d_ref = []
  d_mid = []
  d_new_in_ref = []
  d_new_diff = []
  d_far = []
  ratio = []
  score = []
  model_id = []
  d_min = []
  b_factor = []
  predictions = []
  probabilities = []
  imported_ptms = (chain_id, resid, resname, goto_atom, short_name, full_name, cc, 
    d_ref, d_mid, d_new_in_ref, d_new_diff, d_far, ratio, score, model_id, d_min, b_factor,
    predictions, probabilities)
  with open(ptms_outfile, "r") as out:
    for line in out.readlines():
      items = line.split()
      for array, coerce_lambda in zip(imported_ptms,
        (lambda x: x, lambda x: int(x), lambda x: x, lambda x: x, lambda x: x, lambda x: x,
          lambda x: float(x), lambda x: float(x), lambda x: float(x), lambda x: float(x),
          lambda x: float(x), lambda x: float(x), lambda x: float(x), lambda x: float(x),
          lambda x: int(x), lambda x: float(x), lambda x: float(x),
          lambda x: int(x), lambda x: float(x))):
        array.append(coerce_lambda(items.pop(0)))
  return imported_ptms

def import_synthetic_ptms(synthetic_ptms_outfile):
  chain_id = []
  resid = []
  resname = []
  goto_atom = []
  short_name = []
  full_name = []
  model_id = []
  d_min =  []
  b_factor = []
  imported_synthetic_ptms = (chain_id, resid, resname, goto_atom, short_name, full_name,
    model_id, d_min, b_factor)
  with open(synthetic_ptms_outfile, "r") as out:
    for line in out.readlines():
      items = line.split()
      for array, coerce_lambda in zip(imported_synthetic_ptms,
        (lambda x: x, lambda x: int(x), lambda x: x, lambda x: x, lambda x: x, lambda x: x,
          lambda x: int(x), lambda x: float(x), lambda x: None if x == 'None' else float(x))):
        array.append(coerce_lambda(items.pop(0)))
  return imported_synthetic_ptms
