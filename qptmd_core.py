from __future__ import division

# QTPMD: Iterate over the model supplied and find locations of possible
# posttranslational modification. 

from cctbx import crystal
from scitbx.array_family import flex
from ptm_util import PTM_lookup, PTM_reverse_lookup, get_cc_of_residue_to_map
from map_util import get_fcalc_map

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
    self.mapdata = self.emmap.data.as_double()
    self.frac_matrix = self.emmap.grid_unit_cell().fractionalization_matrix()
    self.ucell_params = emmap.unit_cell_parameters
    self.symmetry = crystal.symmetry(
      space_group_symbol="P1",
      unit_cell=emmap.unit_cell_parameters)
    self.params = params
    # keep a list of which PTMs are possible, let the user examine the results,
    # and be able to go back and make the changes when supplied this list
    self.identified_ptms = []
    self.ccs = []
    # if user specified params.selected_ptms, try to read that file
    if self.params.selected_ptms:
      self.read_selected_ptms()
    # if a score threshold is supplied, store that (default None)
    self.score_threshold = self.params.score_threshold
    # set up the fcalc map to use in scoring residue fits
    self.fcalc_map = get_fcalc_map(
      model_in, self.symmetry, self.params.d_min, self.emmap, scatterer="electron")
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
          self.identified_ptms.append(
            " ".join([chain_id, resid, resname, ptm[0], " ".join(map(str, [ptm[2], ptm[3], ptm[4]]))]))
        i += 1

  def step(self, chain_id, chain, resid, resname, residue, i, type="protein"):
    """step through atoms of the nucleic acid or protein residue, searching for the
    relevant PTM (if any) on each.
    # default: for ptm_mod in PTM[residue], check whether the ptm is present.
    # if self.selected_ptms: make the selected modifications to the model in place.
    # if self.score_threshold: make the best-scoring modification meeting this
    # threshold. Again, apply to the model in-place.
    # in any event: return the qualifying modifications to store
    """
    # if user has supplied specific ptms to apply, only step through those
    if hasattr(self, "selected_ptms"):
      sel = self.selected_ptms[chain_id][resid]
      if not sel:
        return
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
    # get the map-model CC and only test for PTMs if this exceeds a threshold
    cc = get_cc_of_residue_to_map(
      residue, self.frac_matrix, self.ucell_params, self.mapdata, self.fcalc_map)
    self.ccs.append(cc)
    if not cc >= 0.7:
      return []
    # use the lookup tables to iterate through all available modifications
    for (ptm_abbr, ptm_dict) in PTM_lookup[ref_resname].iteritems():
      if ptm_abbr == "unmodified": # for lsq-fitting purposes only
        continue
      elif sel is not None and ptm_dict["name"] not in sel:
        # user requested specific ptms and this is not among them
        continue
      # if a given posttranslational modification is evidenced, append the full
      # name to ptms to access later
      result = self.test_ptm(residue, ptm_dict)
      if result is not None:
        ptms.append(result)
    # depending on the mode of operation, curate the ptms accumulated and apply
    # to the model in place if requested
    def replace(fitted_modded):
      chain.remove_residue_group(residue)
      chain.insert_residue_group(i, fitted_modded)
    if self.score_threshold is not None and len(ptms) > 0:
      # select and apply the ptm with the highest score.
      selected = ptms[0]
      if len(ptms) > 1:
        for current in ptms[1:]:
          if ptms[-1] > selected[-1]:
            selected = current
      fitted_modded = selected[1]
      replace(fitted_modded)
      return [selected]
    elif sel is not None and len(ptms) > 0:
      # user requested specific ptms including the one(s) just identified.
      # right now only implemented for a single ptm per residue.
      if len(ptms) > 1:
        raise NotImplemented, "Only one modification per residue is supported."
      fitted_modded = ptms[0]
      replace(fitted_modded)
    return ptms

  def test_ptm(self, residue, ptm_dict):
    """test whether the given posttranslational modification is supported by the
    electron/charge density at the position of the residue. This method may access
    self.hier to check, for example, whether hydrogen bonding is possible between
    pseudouridine and a neighboring moiety, or may only use densities already stored
    for the atoms of the residue."""
    # if user has supplied specific ptms to apply, only apply those
    fitted_modded = ptm_dict["modify_lambda"](residue.detached_copy())
    (d_ref, d_ptm, ratio) = ptm_dict["ratio_lambda"](
      self.hier, self.mapdata, self.frac_matrix, fitted_modded)
    score = ptm_dict["score_lambda"](self.hier, fitted_modded, d_ref, d_ptm, ratio)
    if self.score_threshold is None or score >= self.score_threshold:
      return (ptm_dict["name"], fitted_modded, d_ref, d_ptm, score)

  def write_identified_ptms(self):
    with open("ptms.out", "wb") as out:
      for line in self.identified_ptms:
        out.write(line + "\n")
    print "Map densities for posttranslational modifications possible in this " +\
      "model\nwritten to file ptms.out. Columns are chain_id, resid, resname, " +\
      "density1,\ndensity2, score. Please curate this file and rerun this " +\
      "program with\nselected_ptms=ptms.out to apply these modifications."

  def write_ccs(self):
    with open("ccs.out", "wb") as out:
      for cc in self.ccs:
        out.write("%f\n" % cc)

  def read_selected_ptms(self):
    """read self.params.selected_ptms into memory"""
    self.selected_ptms = {}
    with open(self.params.selected_ptms, "rb") as ref:
      for line in ref.readlines():
        chain_id, resid, resname, ptm, d1, d2, score = line.split()
        if chain_id not in self.selected_ptms.keys():
          self.selected_ptms[chain_id] = {}
        if resid not in self.selected_ptms[chain_id].keys():
          self.selected_ptms[chain_id][resid] = []
        self.selected_ptms[chain_id][resid].append(ptm)
        if len(self.selected_ptms[chain_id][resid]) > 1:
          raise NotImplemented, "Only one modification per residue is supported."




