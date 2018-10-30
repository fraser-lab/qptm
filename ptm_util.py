from __future__ import division

from scitbx.array_family import flex
from scitbx import matrix
from map_util import get_density_at_position

"""Lookups for distances: expect the given modifications at these
distances. Placing atoms this way is much faster than LSQ fitting
entire residues.
"""
methyl_distances = {"C":10,"O":10,"N":10} # FIXME: obviously update this

def place_modified(modified, unmodified, ref_residue):
  """find the transformation that places the umodified residue at the
  position of the reference residue in the model, and apply that
  transformation to the modified residue"""
  # use mmtbx/math/superpose.py. It also reports the transformation.
  pass

def resi_from_selected_atoms(ref_residue, atom_names):
  """copy ref_residue and produce a copy containing only the atoms
  specified. This can be used for a superposition of the part of a
  residue between the backbone and a PTM, excluding a flexible part."""
  # deep copy ref_residue and trim every atom not in atom_names
  # return the new residue and use this in position_lambdas
  pass

def average_position(*atoms):
  """arguments in atoms should be coordinates (x,y,z). This function
  returns the average position as a scitbx matrix.col."""
  listx, listy, listz = zip(*atoms)
  avg = []
  for l in listx, listy, listz:
    avg.append(sum(l)/len(l))
  return matrix.col(avg)

def place_methyl(methylated_atom, *reference_atoms, on="C"):
  """return the position where a methyl group should lie, at the
  standard methyl bond distance from the methylated_atom, in the
  direction of the vector from the average of the positions of
  the reference_atoms to the position of the methylated atom."""
  methyl_distance = methyl_distances[on]
  reference_position = average_position(*reference_atoms)
  direction = (methylated_atom - reference_position).normalize()
  methyl_position = methyl_distance * direction + methylated_atom
  return methyl_position

def densities_and_ratio(map, atom1, atom2):
  """return (density_at_atom1, density_at_atom2, ratio(2:1))"""
  d1 = get_density_at_position(atom1)
  d2 = get_density_at_position(atom2)
  r = d2/d1
  return (d1, d2, r)

"""NA_PTM: dict of posttranslational modifications known in the
nucleic acids. (FIXME: Will we need to differentiate between DNA and RNA?)

top level: nucleic acids
second level: modifications at the selected atom
third level: the full posttranslational modification name, molecular
structure (if it's a solitary fragment), structure of the complete
modified residue (for some we will substitute this for the model),
and methods for placing it in the map and computing its fit score
"""
NA_PTM = {
  "A":{},
  "U":{
    "unmodified":(,), # phenix-style structure of the unmodified residue
    "pseudouridine":{ # example of a modification to existing atoms
      "name":"pseudouridine",
      "model":None,
      "model_on_RNA":(,), # model of the complete pseudouridine
      "position_lambda":lambda residue:fitted_modded,
      # some lambda function that will align the unmodified uridine
      # to the one in the model and then apply that transformation
      # to this modified one to score. This is a lambda so that we can
      # align to the relevant (usually, modified) part of the residue.
      # *where possible, this should just be an extension along a vector!*
      "ratio_lambda":lambda model, map, fitted_modded:(0, 0, 0)
      # there is no density threshold since no new atoms are present in
      # this particular case. Ignore this for now.
      "score_lambda":lambda model, map, fitted_modded, ratio:score
      # for this one we want to look for a hydrogen bonding partner in the model
    }
  },
  "C":{},
  "G":{
    "unmodified":(,), # phenix-style structure of the unmodified residue
    "m^7":{ # example of an addition of new atoms only (ignoring H)
      "name":"N7-Methylguanosine",
      "model":None, # phenix-style hierarchical model of the modification,
      "model_on_RNA":(,), # model including the full nucleotide
      "position_lambda":lambda residue:fitted_modded,
        use place_methyl(methylated_atom, *reference_atoms)
      # some lambda function that will align the unmodified uridine
      # to the one in the model and then apply that transformation
      # to this modified one to score. This is a lambda so that we can
      # align to the relevant (usually, modified) part of the residue.
      "ratio_lambda":lambda model, map, fitted_modded:\
        densities_and_ratio(map, atom1, atom2)
        # returns tuple (density_at_atom1, density_at_atom2, ratio2:1)
      # *store the densities -- let the user analyze a distribution*
      "score_lambda":lambda model, map, fitted_modded, ratio:score
      # some lambda function to decide whether the PTM at this site looks
      # reasonably well evidenced by the map, as reported by this score,
      # hopefully usually just based on the ratio in most cases
    }
  }
}

prot_PTM = {} # not implemented


