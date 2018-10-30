from __future__ import division

# QTPMD: Iterate over the model supplied and find locations of possible
# posttranslational modification. 

from ptm_lookup import NA_PTM, prot_PTM
from map_util import get_density_at_position

class LookForPTMs(object):
  """Step through the hierarchical molecular model looking for evidence of each
  known posttranslational modification at each appropriate site, and add this
  information to the model. Also supply methods for accessing this afterward."""
  def __init__(self, hierarchical_molec_model, emmap, params=None):
    self.hier = hierarchical_molec_model
    self.emmap = emmap
    self.params = params
    self.score_threshold = 0.2
    self.identified_ptms = []
    self.walk()

  def walk(self):
    """walk through the molecular model, stepping by one nucleotide or residue,
    searching for PTM at each site"""
    # roughly: for chain in model, do for residuein chain, do step(residue)
    for chain in self.hier.chains():
      struct_type = "protein" if chain.is_protein() else "na"
      for resi in chain.residues(): # same accessor for protein and nucleic acids
        self.step(resi, type-struct_type)

  def step(self, residue, type="protein"):
    """step through atoms of the nucleic acid or protein residue, searching for the
    relevant PTM (if any) on each"""
    # roughly: for atom in residue, do calculate density at atom and store.
    # for ptm_mod in NA_PTM[residue], check whether the modification is present.
    # later plan to have the option to expand beyond RNA PTMs to protein and DNA.
    resname = residue.resname.strip()
    residue.ptms = []
    for atom in residue.atoms():
      atom.density = get_density_at_position(atom.xyz, self.map)
    # skip the calc per atom -- let the lambda do it only when necessary
    PTM_lookup = prot_PTM if type == "protein" else NA_PTM
    for (ptm_name, ptm_dict) in PTM_lookup[resname].iteritems():
      if ptm_name == "unmodified": # for overlay purposes
        continue
      # if a given posttranslational modification is evidenced, append the full
      # name to residue.ptms to access later
      self.test_ptm(residue, ptm_name, ptm_dict)

  def test_ptm(self, residue, ptm_dict):
    """test whether the given posttranslational modification is supported by the
    electron/charge density at the position of the residue. This method may access
    self.hier to check, for example, whether hydrogen bonding is possible between
    pseudouridine and a neighboring moiety, or may only use densities already stored
    for the atoms of the residue."""
    fitted_modded = ptm_dict["position_lambda"](residue)
    (d_ref, d_ptm, ratio) = ptm_dict["ratio_lambda"](
      self.hier, self.emmap, fitted_modded)
    if ptm_dict["score_lambda"](
      self.hier, self.emmap, fitted_modded, ratio) >= self.score_threshold:
      residue.ptms.append((ptm_dict["name"], ))

  def report_ptms(self):
    """print all posttranslational modifications identified by the walk method"""
    for chain in self.hier.chains():
      for resi in chain.residues():
        if resi.ptms: # one or more ptms identified
          self.identified_ptms.append(chain.id + resi.resseq + resi.resname + \
            ": " + "\n     ".join(resi.ptms))
    print "Density ratios for posttranslational modifications possible in this model:\n"
    print "\n".join(self.identified_ptms)

