from __future__ import division
from scitbx.array_family import flex
from ptm_util import PTM_lookup, PTM_reverse_lookup, unmodified_residues

class ModifiedModel(object):
  """Class for models that may have already been modified, for examination
  independently of any map."""
  def __init__(self, hier_model):
    self.hier = hier_model
    self.modeled_ptms = []
    from random import randint
    self.id = randint(0,999999)
  def check_modeled_ptms(self, model_id=None, d_min=None, b_factor=None):
    """walk through the molecular model checking what PTMs are already modeled,
    and write the modeled ones out to a new file modeled_ptms.out"""
    # return [] # TMP FIXME
    with open("modeled_ptms.out", "wb") as outfile:
      assert not None in (model_id, d_min, b_factor)
      constant_str = " %d %f %f" % (model_id, d_min, b_factor)
      for chain in self.hier.chains():
        chain_id = chain.id.strip()
        struct_type = "protein" if chain.is_protein() else "na"
        # below, we use the residue_groups() instead of the residues() accessor so
        # we don't lose the connection to the parent model. It's a little less
        # convenient to use but being able to copy residues is worth it. Also note:
        # same accessor for protein and nucleic acids.
        for residue in chain.residue_groups():
          resid = residue.resid().strip()
          resname = residue.unique_resnames()[0].strip()
          ptm_tuple = self.check_residue(chain_id, resid, resname)
          if ptm_tuple:
            self.modeled_ptms.append(ptm_tuple)
            outfile.write(" ".join(ptm_tuple) + constant_str + "\n")
    return self.modeled_ptms
  def check_residue(self, chain_id, resid, resname):
    """determine whether the resname matches a known ptm, and if so, return a
    tuple of the same values that we'd write to ptms.out"""
    if resname in PTM_lookup.keys():
      return
    else:
      try:
        ref_resname = PTM_reverse_lookup[resname]
      except KeyError:
        # print "skipping unrecognized residue", resname
        return
    ptm_dict = PTM_lookup[ref_resname][resname]
    return (chain_id, resid, ref_resname, ptm_dict["goto_atom"], ptm_dict["name"])
  def clean_up(self, prune=True, b_factor=None, filename="clean.pdb"):
    """walk through the molecular model and remove any identified posttranslational
    modifications using the prune_lambda function associated with the modification
    in the PTM_lookup dictionary. Reset B factors to specified value if requested."""
    skipped = ["HOH", "WAT"]
    for chain in self.hier.chains():
      for residue in chain.residue_groups():
        if b_factor is not None:
          atoms = residue.atom_groups()[0].atoms()
          atoms.set_b(flex.double(len(atoms), b_factor))
        resname = residue.unique_resnames()[0].strip()
        if prune:
          if resname in unmodified_residues:
            continue
          elif resname in PTM_reverse_lookup.keys():
            pruned_resname = PTM_reverse_lookup[resname]
            PTM_lookup[pruned_resname][resname]["prune_lambda"](residue)
            for ag in residue.atom_groups():
              ag.resname = pruned_resname
          else:
            if resname not in skipped:
              print "Warning: skipping unrecognized residue, ligand or ion %s" % resname
              skipped.append(resname)
    self.hier.write_pdb_file(filename)

