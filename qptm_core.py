from __future__ import division

# QTPM: Iterate over the model supplied and find locations of possible
# posttranslational modifications.

from cctbx import crystal
from scitbx.array_family import flex
from ptm_util import PTM_lookup, PTM_reverse_lookup, get_cc_of_residue_to_map, prune_confs, rename
from map_util import get_fcalc_map, get_diff_map, write_ccp4_map
from libtbx.utils import Sorry
import math

class LookForPTMs(object):
  """Step through the hierarchical molecular model looking for evidence of each
  known posttranslational modification at each appropriate site, and add this
  information to the identified_ptms list (as it is not persistent in the model).
  If a score threshold is provided, apply the best-scoring modifications meeting the
  threshold at each possible site, or if a list of modifications is supplied, make
  those modifications and write out the updated model."""
  def __init__(self, model_in, hierarchical_molec_model, emmap, model_id=None, diff_map=None, calc_map=None, modeled_ptms=None, params=None):
    self.hier = hierarchical_molec_model
    self.modeled_ptms = modeled_ptms
    def data_from_map(map_obj):
      data = map_obj.data*(len(map_obj.data)/flex.sum(flex.abs(map_obj.data)))
      map_obj.data = data
      return data.as_double()
    self.emmap = emmap
    self.mapdata = data_from_map(self.emmap)
    self.frac_matrix = self.emmap.grid_unit_cell().fractionalization_matrix()
    self.ucell_params = self.emmap.unit_cell_parameters
    self.symmetry = crystal.symmetry(
      space_group_symbol="P1",
      unit_cell=self.ucell_params)
    self.params = params
    if model_id is None:
      from random import randint
      self.id = randint(0,999999)
    else:
      self.id = model_id
    # keep a list of which PTMs are possible, let the user examine the results,
    # and be able to go back and make the changes when supplied this list
    self.identified_ptms = []
    self.all_tested_ptms = []
    self.ccs = []
    # if user specified params.selected_ptms, try to read that file
    if self.params.selected_ptms:
      self.read_selected_ptms()
    # if user specified params.synthetic_data, randomly modify a copy of the model
    # and generate fake map data matching those modifications.
    if self.params.synthetic_data:
      self.generate_modified_model()
      self.mapdata = get_fcalc_map(self.synthetic_model_object, self.symmetry, self.params.d_min,
        self.emmap, scatterer=self.params.experiment)
    # set up the fcalc map to use in scoring residue fits
    self.fcalc_map = data_from_map(calc_map) if calc_map else get_fcalc_map(
      model_in, self.symmetry, self.params.d_min, self.emmap, scatterer=self.params.experiment)
    self.diff_map = data_from_map(diff_map) if diff_map else get_diff_map(
      self.symmetry, self.fcalc_map, self.mapdata, self.params.d_min)
    self.test_count = 0
    self.walk()
    self.filter_ptms()
    # if user specified params.synthetic_data, write out the true modifications to judge against
    # and report on this for the user
    if self.params.synthetic_data:
      self.write_synthetic_ptms()
      self.report_accuracy(reference=self.synthetic_ptms)
    else:
      self.report_accuracy(reference=self.modeled_ptms)

  def walk(self):
    """walk through the molecular model, stepping by one nucleotide or residue,
    searching for PTM at each site"""
    # roughly: for chain in model, do for residue in chain, do step(residue)
    for chain in self.hier.chains():
      chain_id = chain.id.strip()
      struct_type = "protein" if chain.is_protein() else "na"
      i = 0 # index of the residue in the chain object
      # below, we use the residue_groups() instead of the residues() accessor so
      # we don't lose the connection to the parent model. It's a little less
      # convenient to use but being able to copy residues is worth it. Also note:
      # same accessor for protein and nucleic acids.
      for residue in chain.residue_groups():
        resid = residue.resid().strip()
        resname = residue.unique_resnames()[0].strip()
        ptms, cc = self.step(chain_id, chain, resid, resname, residue, i, type=struct_type)
        for ptm in ptms:
          self.all_tested_ptms.append((chain_id, resid, resname, ptm, cc))
        i += 1

  def step(self, chain_id, chain, resid, resname, residue, i, type="protein"):
    """step through atoms of the nucleic acid or protein residue, searching for the
    relevant PTM (if any) on each.
    # default: for ptm_mod in PTM[residue], check whether the ptm is present.
    # if self.selected_ptms: make the selected modifications to the model in place.
    # return the qualifying modifications to store
    """
    # if user has supplied specific ptms to apply, only step through those
    if hasattr(self, "selected_ptms"):
      try:
        sel = self.selected_ptms[chain_id][resid]
      except KeyError:
        sel = None
    else:
      sel = None
    ptms = []
    # get the map-model CC
    cc = get_cc_of_residue_to_map(
      residue, self.frac_matrix, self.ucell_params, self.mapdata, self.fcalc_map)
    # check whether this is a modified residue already. Find modifications based
    # on the unmodified residue. # FIXME: the atom names won't be the same for
    # some cases like pseudouridine. Do we need to add these to the main lookup
    # so that they have their own custom lambdas? Probably that would be the best
    # way to let the test atom positions match the modeled positions anyway.
    if resname in PTM_lookup.keys():
      ref_resname = resname
    else:
      try:
        ref_resname = PTM_reverse_lookup[resname]
      except KeyError:
        # print "skipping unrecognized residue", resname
        return ([], cc)
    self.test_count += len([k for k in PTM_lookup[ref_resname].keys() if k != "unmodified"])
    # use the lookup tables to iterate through all available modifications
    for (ptm_code, ptm_dict) in PTM_lookup[ref_resname].iteritems():
      if ptm_code == "unmodified": # for lsq-fitting purposes only
        continue
      elif sel is not None and ptm_dict["name"] not in sel:
        # user requested specific ptms and this is not among them
        continue
      # if a given posttranslational modification is evidenced, append the full
      # name to ptms to access later
      result = self.test_ptm(residue, ptm_dict, ptm_code)
      if result is not None:
        ptms.append(result)
    # depending on the mode of operation, curate the ptms accumulated and apply
    # to the model in place if requested
    def remove(residue):
      chain.remove_residue_group(residue)
    def place(fitted_modded):
      prune_confs(self.mapdata, self.frac_matrix, fitted_modded)
      # rename(fitted_modded, ptm_code[:3])
      chain.insert_residue_group(i, fitted_modded)
    if len(ptms) > 0:
      if self.params.apply_best:
        # determine the best scoring version of each modification if requested
        if len(ptms) == 1:
          keep_idx = 0
        else:
          mods = {name:[p for p in ptms if p[2] == name] for name in set(zip(*ptms)[2])}
          results = []
          for modname, modlist in mods.iteritems():
            scores = zip(*modlist)[9]
            keep_idx = scores.index(max(scores))
            results.append(modlist[keep_idx])
          ptms = results
      # remove the resi and replace it with the new one(s)
      remove(residue)
      if self.params.apply_best:
        place(ptms[keep_idx][0]) # keep only the best scoring one if requested
      else:
        for fitted_modded in ptms:
          place(fitted_modded[0]) # otherwise place all possible mods for visualization
    return (ptms, cc)

  def test_ptm(self, residue, ptm_dict, ptm_code):
    """test whether the given posttranslational modification is supported by the
    electron/charge density at the position of the residue. This method may access
    self.hier to check, for example, whether hydrogen bonding is possible between
    pseudouridine and a neighboring moiety, or may only use densities already stored
    for the atoms of the residue."""
    try:
      # if user has supplied specific ptms to apply, only apply those
      fitted_modded = ptm_dict["modify_lambda"](residue.detached_copy())
      (d_ref, d_new_diff, d_mid, d_new_in_ref, d_far, ratio) = ptm_dict["ratio_lambda"](
        self.hier, self.mapdata, self.diff_map, self.frac_matrix, fitted_modded)
      # discard any cases where the density shape doesn't match a single protrusion
      # --> this has been moved to the filter_ptms step
      # then filter by score
      # --> this filter has also been moved, but compute it here
      scaled_ratio = ptm_dict["scale_ratio_lambda"](self.hier, fitted_modded, d_ref, d_new_in_ref, ratio)
      score = scaled_ratio * math.log(d_new_diff)
      # if score >= self.params.score_threshold:
      #   # rename(fitted_modded, ptm_code[:3]) FIXME FIND THE RIGHT COOT SETTING TO ENABLE
      return (fitted_modded, ptm_dict["goto_atom"], ptm_dict["name"],
        d_ref, d_mid, d_new_in_ref, d_new_diff, d_far, scaled_ratio, score, "")
      # last emtpy string is the log of reasons we've rejected a modification, but we've
      # moved all these steps to filter_ptms so we don't have any to report yet
    except Sorry:
      return

  def filter_ptms(self):
    """Filter ptms based on the correlation coefficient threshold. Then,
    remove suggested ptms for which the reference atoms' density is below a
    given threshold. Absolute value is likely to be map-dependent so use a
    proportion instead. Then also remove suggested ptms for which the difference
    density for the PTM is below another user-supplied threshold."""
    from filter_util import apply_filters
    # make everything flex arrays and formatted as expected
    chain_id, resid, resname, ptm, cc = zip(*self.all_tested_ptms)
    flex_chain_id = flex.std_string(chain_id)
    flex_resid = flex.int(map(int, resid))
    flex_resname = flex.std_string(resname)
    flex_cc = flex.double(map(float, cc))
    length = len(flex_cc)
    # further unpack ptm
    fitted_modded, goto_atom, name, d_ref, d_mid, d_new_in_ref, d_new_diff, d_far, \
      scaled_ratio, score, log = zip(*ptm)
    # fitted_modded doesn't fit in a flex array -- we'll have to use list comprehensions
    # in a later step to apply the selectin we determine to the list containing this
    flex_goto_atom = flex.std_string(goto_atom)
    short_name, full_name = zip(*[n.split() for n in name])
    flex_short_name = flex.std_string(short_name)
    flex_full_name = flex.std_string(full_name)
    flex_d_ref = flex.double(map(float, d_ref))
    flex_d_mid = flex.double(map(float, d_mid))
    flex_d_new_in_ref = flex.double(map(float, d_new_in_ref))
    flex_d_new_diff = flex.double(map(float, d_new_diff))
    flex_d_far = flex.double(map(float, d_far))
    flex_scaled_ratio = flex.double(map(float, scaled_ratio))
    flex_score = flex.double(map(float, score))
    flex_model_id = flex.int(length, self.id)
    flex_d_min = flex.double(length, self.params.d_min)
    flex_b_factor = flex.double(length, self.params.set_b_factor) \
      if self.params.set_b_factor is not None else \
      flex.double([fm.atom_groups()[0].atoms()[0].b for fm in fitted_modded])
    # skip log, we'll overwrite it
    import_format_ptms = (flex_chain_id, flex_resid, flex_resname, flex_goto_atom,
      flex_short_name, flex_full_name, flex_cc, flex_d_ref, flex_d_mid, flex_d_new_in_ref,
      flex_d_new_diff, flex_d_far, flex_scaled_ratio, flex_score, flex_model_id,
      flex_d_min, flex_b_factor)
    keep_selection, accepted, all_tested, log = apply_filters(
      import_format_ptms,
      cc_threshold=self.params.cc_threshold,
      ref_frac=self.params.reference_densities_fraction,
      dif_frac=self.params.difference_densities_fraction,
      ratio_d_far_d_new_in_ref=self.params.ratio_d_far_d_new_in_ref,
      ratio_d_far_d_mid=self.params.ratio_d_far_d_mid,
      ratio_d_ref_d_new_in_ref=self.params.ratio_d_ref_d_new_in_ref,
      score_threshold=self.params.score_threshold)
    # save updated selections and logs
    self.identified_ptms = [self.all_tested_ptms[i] for i in
      xrange(len(self.all_tested_ptms)) if keep_selection[i]]
    updated_ptm = zip(fitted_modded, goto_atom, name, d_ref, d_mid, d_new_in_ref,
      d_new_diff, d_far, scaled_ratio, score, list(log))
    updated_all = zip(chain_id, resid, resname, updated_ptm, cc)
    self.all_tested_ptms = updated_all
    print "total possible ptms tested: %d" % self.test_count

  def generate_modified_model(self):
    """make a copy of the model and add random posttranslational modifications to
    10% of the available sites to generate synthetic training data. keep track of
    the modifications in self.synthetic_ptms (similar format to self.ptms)"""
    self.synthetic_ptms = []
    self.synthetic_model = self.hier.deep_copy()
    self.modify_model_at_random(self.synthetic_model)
    self.synthetic_model.write_pdb_file("synthetic.pdb")
    from iotbx import file_reader
    self.synthetic_model_object = file_reader.any_file("synthetic.pdb", force_type="pdb").file_object

  def modify_model_at_random(self, model_copy):
    """walk through the molecular model, modifying 10% of the recognized residues"""
    # several components borrowed from walk method
    from random import random
    # for testing purposes:
    # from random import seed
    # seed(32)
    for chain in self.synthetic_model.chains():
      chain_id = chain.id.strip()
      struct_type = "protein" if chain.is_protein() else "na"
      i = 0 # index of the residue in the chain object
      # below, we use the residue_groups() instead of the residues() accessor so
      # we don't lose the connection to the parent model. It's a little less
      # convenient to use but being able to copy residues is worth it. Also note:
      # same accessor for protein and nucleic acids.
      for residue in chain.residue_groups():
        resid = residue.resid().strip()
        resname = residue.unique_resnames()[0].strip()
        if self.params.set_b_factor is not None:
          atoms = residue.atom_groups()[0].atoms()
          atoms.set_b(flex.double(len(atoms), self.params.set_b_factor))
        if random() >= 0.9:
          try:
            self.modify_residue_at_random(
              chain_id, chain, resid, resname, residue, i, type=struct_type)
          except Sorry:
            continue
        i += 1

  def modify_residue_at_random(self, chain_id, chain, resid, resname, residue, i, type="protein"):
    """for a given residue, apply a random modification."""
    # figure out what we're looking at and how we know how to modify it
    if resname in PTM_lookup.keys():
      ref_resname = resname
    else:
      try:
        ref_resname = PTM_reverse_lookup[resname]
        # this residue has already been modified, but that's okay, we'll modify it some more!
      except KeyError:
        # for whatever reason, we don't recognize this residue, so leave it alone
        return
    # constructing a new dictionary is safer than referencing the dictionary and deleting
    # "unmodified", because that change would persist outside this method
    modifications_available = {key:value for key, value in PTM_lookup[ref_resname].iteritems()
      if key != "unmodified"}
    # select a random modification from those available
    from random import randint
    idx = randint(0,len(modifications_available)-1)
    key = modifications_available.keys()[idx]
    ptm_dict = modifications_available[key]
    # place that modification
    fitted_modded = ptm_dict["modify_lambda"](residue.detached_copy())
    rename(fitted_modded, key[:3])
    chain.remove_residue_group(residue)
    chain.insert_residue_group(i, fitted_modded)
    self.synthetic_ptms.append(
      (chain_id, resid, resname, ptm_dict["goto_atom"], ptm_dict["name"]))

  def report_accuracy(self, verbose=False, reference=None):
    """Look over the placed and identified modifications and report on our success rate."""
    if reference is None:
      raise Sorry, "report accuracy relative to what reference?"
    true_positives, false_positives, false_negatives = [],[],[]
    ptms_organized = {} # by chain, resi, goto_atom, and then mod_name
    for (chain_id, resid, resname, ptm, cc) in self.identified_ptms:
      goto_atom = ptm[1]
      mod_name = ptm[2]
      record = " ".join([chain_id, resid, goto_atom, mod_name])
      if not chain_id in ptms_organized:
        ptms_organized[chain_id] = {}
      if not resid in ptms_organized[chain_id]:
        ptms_organized[chain_id][resid] = {}
      if not goto_atom in ptms_organized[chain_id][resid]:
        ptms_organized[chain_id][resid][goto_atom] = []
      if mod_name in ptms_organized[chain_id][resid][goto_atom]:
        # this removes duplicates: multiple possible methylations at one position
        continue
      ptms_organized[chain_id][resid][goto_atom].append(mod_name)
      false_positives.append(record) # put everything in and then remove in next step
    for (chain_id, resid, resname, goto_atom, mod_name) in reference:
      record = " ".join([chain_id, resid, goto_atom, mod_name])
      if chain_id in ptms_organized and resid in ptms_organized[chain_id] and \
        goto_atom in ptms_organized[chain_id][resid] and \
        mod_name in ptms_organized[chain_id][resid][goto_atom]:
        # it was a true positive -- update
        del false_positives[false_positives.index(record)]
        true_positives.append(record)
      else:
        false_negatives.append(record)
    print "True positives count: %d\nFalse positives count: %d\nFalse negatives count: %d"%\
      (len(true_positives), len(false_positives), len(false_negatives))
    print "True negatives count: %d\n"%\
      (len(self.all_tested_ptms) - sum([len(true_positives), len(false_positives), len(false_negatives)]))
    if verbose:
      # print "True positives:"
      # print "\n".join(true_positives)
      print "\nFalse positives:"
      print "\n".join(false_positives)
      print "\nFalse negatives:"
      print "\n".join(false_negatives)
    print "\n"

  def write_synthetic_ptms(self):
    with open("synthetic_ptms.out", "wb") as out:
      for (chain_id, resid, resname, goto_atom, mod_name) in self.synthetic_ptms:
        out.write(" ".join([chain_id, resid, resname, goto_atom, mod_name,
          " ".join(map(str, [self.id, self.params.d_min, self.params.set_b_factor]))
          ]) + "\n")

  def write_identified_ptms(self):
    with open("ptms.out", "wb") as out:
      for (chain_id, resid, resname, ptm, cc) in self.identified_ptms:
        out.write(" ".join([chain_id, resid, resname, ptm[1], ptm[2], str(cc),
            " ".join(map(str, [ptm[i] for i in xrange(3,10)])),
            " ".join(map(str, [self.id, self.params.d_min, self.params.set_b_factor]))
            ]) + "\n")
        # ptm is (fitted_modded, ptm_dict["goto_atom"], ptm_dict["name"],
        # d_ref, d_mid, d_new_in_ref, d_new_diff, d_far, scaled_ratio, score)
        # we will write out (chain_id, resid, resname, goto_atom, short_name, full_name,
        # cc, d_ref, d_mid, d_new_in_ref, d_new_diff, d_far, scaled_ratio, score)
    print """\nMap densities for posttranslational modifications suggested for this model
written to file ptms.out. Columns are chain_id, resid, resname, atom, ptm
abbreviation, modified resname, density1, density2, score. Please curate
this file and rerun this program with selected_ptms=ptms.out to apply these
modifications.
"""

  def write_all_tested_ptms(self):
    with open("all_tested_ptms.out", "wb") as out:
      for (chain_id, resid, resname, ptm, cc) in self.all_tested_ptms:
        out.write(" ".join([chain_id, resid, resname, ptm[1], ptm[2], str(cc),
              " ".join(map(str, [ptm[i] for i in xrange(3,10)])),
              " ".join(map(str, [self.id, self.params.d_min, self.params.set_b_factor])),
              ptm[10]]) + "\n")

  def write_ccs(self):
    with open("ccs.out", "wb") as out:
      out.write("\n".join(self.ccs))

  def read_selected_ptms(self):
    """read self.params.selected_ptms into memory"""
    self.selected_ptms = {}
    with open(self.params.selected_ptms, "rb") as ref:
      for line in ref.readlines():
        chain_id, resid, resname, goto_atom, ptm_abbr, ptm_longname = line.split(" ")[:6]
        ptm = " ".join([ptm_abbr, ptm_longname])
        if chain_id not in self.selected_ptms.keys():
          self.selected_ptms[chain_id] = {}
        if resid not in self.selected_ptms[chain_id].keys():
          self.selected_ptms[chain_id][resid] = []
        self.selected_ptms[chain_id][resid].append(ptm)
        # if len(self.selected_ptms[chain_id][resid]) > 1:
        #   raise NotImplementedError, "Only one modification per residue is supported."

  def write_modified_model(self, filename="modified.pdb"):
    """write a modified model back out if supplied with selected ptms to apply"""
    self.hier.write_pdb_file(filename)
    print "\nSelected modifications in the provided file written to %s." % filename

  def write_difference_map(self, filename="difference_map.ccp4"):
    """write the calculated difference map out"""
    write_ccp4_map(self.diff_map, self.symmetry, filename)

  def write_synthetic_map(self, filename="synthetic_map.ccp4"):
    """write the synthetic map out"""
    write_ccp4_map(self.mapdata, self.symmetry, filename)

  def write_calculated_map(self, filename="calculated_map.ccp4"):
    """write the fcalc map out"""
    write_ccp4_map(self.fcalc_map, self.symmetry, filename)





