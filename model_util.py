def avg_b(residue):
  average = lambda l: sum(l)/len(l)
  return average([a.b for a in residue.atom_groups()[0].atoms()])

        constant_str = " %s %s %s" % tuple(map(str, (model_id, d_min, b_factor)))