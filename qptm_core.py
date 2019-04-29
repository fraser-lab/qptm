from __future__ import division

# QTPMD: Iterate over the model supplied and find locations of possible
# posttranslational modification. 

from cctbx import crystal
from scitbx.array_family import flex
from ptm_util import PTM_lookup, PTM_reverse_lookup, get_cc_of_residue_to_map, prune_confs, rename
from map_util import get_fcalc_map, get_diff_map, write_ccp4_map

class LookForPTMs(object):
  """Step through the hierarchical molecular model looking for evidence of each
  known posttranslational modification at each appropriate site, and add this
  information to the identified_ptms list (as it is not persistent in the model).
  If a score threshold is provided, apply the best-scoring modifications meeting the
  threshold at each possible site, or if a list of modifications is supplied, make
  those modifications and write out the updated model."""
  def __init__(self, model_in, hierarchical_molec_model, emmap, params=None):
    self.hier = hierarchical_molec_model
    self.emmap = emmap
    self.emmap.data = self.emmap.data*(len(self.emmap.data)/flex.sum(flex.abs(self.emmap.data)))
    self.mapdata = self.emmap.data.as_double()
    self.frac_matrix = self.emmap.grid_unit_cell().fractionalization_matrix()
    self.ucell_params = self.emmap.unit_cell_parameters
    self.symmetry = crystal.symmetry(
      space_group_symbol="P1",
      unit_cell=self.ucell_params)
    self.params = params
    # keep a list of which PTMs are possible, let the user examine the results,
    # and be able to go back and make the changes when supplied this list
    self.identified_ptms = []
    self.ccs = []
    # if user specified params.selected_ptms, try to read that file
    if self.params.selected_ptms:
      self.read_selected_ptms()
    # set up the fcalc map to use in scoring residue fits
    self.fcalc_map = get_fcalc_map(
      model_in, self.symmetry, self.params.d_min, self.emmap, scatterer="electron")
      # model_in, self.symmetry, self.params.d_min, self.emmap, scatterer="xray") # FIXME PUT IT BACK TO ELECTRON
    self.diff_map = get_diff_map(
      self.symmetry, self.fcalc_map, self.mapdata, self.params.d_min)
    self.test_count = 0
    self.walk()

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
        ptms = self.step(chain_id, chain, resid, resname, residue, i, type=struct_type)
        for ptm in ptms:
          self.identified_ptms.append((chain_id, resid, resname, ptm))
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
        return []
    self.test_count += len([k for k in PTM_lookup[ref_resname].keys() if k != "unmodified"])
    # get the map-model CC and only test for PTMs if this exceeds a threshold
    cc = get_cc_of_residue_to_map(
      residue, self.frac_matrix, self.ucell_params, self.mapdata, self.fcalc_map)
    self.ccs.append("%s %s %f"%(residue.id_str().strip(), residue.unique_resnames()[0].strip(), cc))
    if not cc >= self.params.cc_threshold:
      return []
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
    if sel is not None and len(ptms) > 1:
      # we are applying the best scoring version of each modification if the user
      # requests this modification (for visualization, we just keep all versions)
      mods = {name:[p for p in ptms if p[1] == name] for name in set(zip(*ptms)[1])}
      results = []
      for modname, modlist in mods.iteritems():
        scores = zip(*modlist)[-1]
        keep_idx = scores.index(max(scores))
        results.append(modlist[keep_idx])
      ptms = results
    if len(ptms) > 0:
      remove(residue)
      for fitted_modded in ptms:
        place(fitted_modded[2]) # yes, this places multiple overlapping copies if asked!
    return ptms

  def test_ptm(self, residue, ptm_dict, ptm_code):
    """test whether the given posttranslational modification is supported by the
    electron/charge density at the position of the residue. This method may access
    self.hier to check, for example, whether hydrogen bonding is possible between
    pseudouridine and a neighboring moiety, or may only use densities already stored
    for the atoms of the residue."""
    # if user has supplied specific ptms to apply, only apply those
    fitted_modded = ptm_dict["modify_lambda"](residue.detached_copy())
    (d_ref, d_new_diff, d_mid, d_new_in_ref, d_far, ratio) = ptm_dict["ratio_lambda"](
      self.hier, self.mapdata, self.diff_map, self.frac_matrix, fitted_modded)
    # discard any cases where the density shape doesn't match a single protrusion
    # TODO this should be controlled by one or two tunable thresholds
    if d_far > d_new_in_ref: return
    if d_far > d_mid: return
    if d_new_in_ref < 0.3 * d_ref: return
    # then filter by score
    score = ptm_dict["score_lambda"](self.hier, fitted_modded, d_ref, d_new_in_ref, ratio)
    if score >= self.params.score_threshold:
      # rename(fitted_modded, ptm_code[:3]) FIXME FIND THE RIGHT COOT SETTING TO ENABLE
      return (ptm_dict["goto_atom"], ptm_dict["name"],
        fitted_modded, d_ref, d_new_diff, score)

  def filter_ptms(self):
    """remove suggested ptms for which the reference atoms' density is below a
    given threshold. Absolute value is likely to be map-dependent so use a
    proportion instead. Then also remove suggested ptms for which the difference
    density for the PTM is below another user-supplied threshold."""
    # calculate keep_density based on keep_fraction
    if self.params.reference_densities_fraction < 1:
      reference_densities = flex.double([tup[3][3] for tup in self.identified_ptms])
      from scitbx.python_utils import robust_statistics as rs
      keep_density = rs.percentile(reference_densities,
        1 - self.params.reference_densities_fraction)
      keep_selection = reference_densities >= keep_density
      self.identified_ptms = [self.identified_ptms[i] for i in 
        xrange(len(self.identified_ptms)) if keep_selection[i]]
    if self.params.difference_densities_fraction < 1:
      difference_densities = flex.double([tup[3][4] for tup in self.identified_ptms])
      from scitbx.python_utils import robust_statistics as rs
      keep_density = rs.percentile(difference_densities,
        1 - self.params.difference_densities_fraction)
      keep_selection = difference_densities >= keep_density
      self.identified_ptms = [self.identified_ptms[i] for i in 
        xrange(len(self.identified_ptms)) if keep_selection[i]]
    print "total possible ptms tested: %d" % self.test_count

  def write_identified_ptms(self):
    with open("ptms.out", "wb") as out:
      for (chain_id, resid, resname, ptm) in self.identified_ptms:
        out.write(" ".join([chain_id, resid, resname, ptm[0], ptm[1],
              " ".join(map(str, [ptm[3], ptm[4], ptm[5]]))]) + "\n")
        # ptm is (ptm_dict["goto_atom"], ptm_dict["name"], d_ref, d_new_diff, score)
    print """\nMap densities for posttranslational modifications suggested for this model
written to file ptms.out. Columns are chain_id, resid, resname, atom, ptm
abbreviation, modified resname, density1, density2, score. Please curate
this file and rerun this program with selected_ptms=ptms.out to apply these
modifications.
"""

  def write_ccs(self):
    with open("ccs.out", "wb") as out:
      out.write("\n".join(self.ccs))

  def read_selected_ptms(self):
    """read self.params.selected_ptms into memory"""
    self.selected_ptms = {}
    with open(self.params.selected_ptms, "rb") as ref:
      for line in ref.readlines():
        chain_id, resid, resname, goto_atom, ptm_abbr, ptm_longname, _, _, _ = line.split(" ")
        ptm = " ".join([ptm_abbr, ptm_longname])
        if chain_id not in self.selected_ptms.keys():
          self.selected_ptms[chain_id] = {}
        if resid not in self.selected_ptms[chain_id].keys():
          self.selected_ptms[chain_id][resid] = []
        self.selected_ptms[chain_id][resid].append(ptm)
        if len(self.selected_ptms[chain_id][resid]) > 1:
          raise NotImplementedError, "Only one modification per residue is supported."

  def write_modified_model(self, filename="modified.pdb"):
    """write a modified model back out if supplied with selected ptms to apply"""
    self.hier.write_pdb_file(filename)
    print "\nSelected modifications in the provided file written to %s." % filename

  def write_difference_map(self, filename="difference_map.ccp4"):
    """write the calculated difference map out"""
    write_ccp4_map(self.diff_map, self.symmetry, filename)





