      if not None in (model_id, d_min, b_factor):
        constant_str = " %d %f %f" % (model_id, d_min, b_factor)
      else:
        constant_str = " %s %s %s" % map(str, (model_id, d_min, b_factor))