from __future__ import division
from scitbx.array_family import flex
from scitbx import matrix
from map_util import get_density_at_position
from libtbx.utils import Sorry

"""Lookups for distances: expect the given modifications at these
distances. Placing atoms this way is much faster than LSQ fitting
entire residues."""
methyl_distances = {"C":1.4,"O":1.4,"N":1.4} # FIXME: more precise?

"""Amino acids and nucleotides: standard and rare (unmodified)"""
std_amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
rare_amino_acids = ['FME', 'MSE', 'SEC', 'PYL', 'XAA', 'UNK'] # includes unknown or other unspecified
std_nucleotides = ['A', 'U', 'T', 'C', 'G']
unmodified_residues = std_amino_acids + rare_amino_acids + std_nucleotides

def locate_atom_by_name(residue, name):
  """find the first copy of an atom by name on a hierarchical
  residue object"""
  atoms = residue.atom_groups()[0].atoms()
  for i in xrange(len(atoms)):
    if atoms[i].name.strip() == name.strip():
      return atoms[i]
  raise Sorry("couldn't locate requested atom")

def locate_atoms_by_name(residue, name):
  """find all copies of an atom by name on a hierarchical residue
  object"""
  matches = []
  atoms = residue.atom_groups()[0].atoms()
  for i in xrange(len(atoms)):
    if atoms[i].name.strip() == name.strip():
      matches.append(atoms[i])
  return matches

def average_position(*atoms):
  """arguments in atoms should be coordinates (x,y,z). This function
  returns the average position as a scitbx matrix.col."""
  try:
    listx, listy, listz = zip(*(a.xyz for a in atoms))
  except AttributeError:
    listx, listy, listz = zip(*atoms)
  avg = []
  for l in listx, listy, listz:
    avg.append(sum(l)/len(l))
  return matrix.col(avg)

def sp3_rot_methylate(residue, methylated_atom, reference_atoms_tup, on="C", name=" C  "):
  """place a new methyl group on the given residue at the
  standard methyl bond distance from the methylated_atom.
  If reference_atoms_tup is three atoms, interpret as an
  sp2 group adjacent to the methylated_atom and return
  the best rotamer for an sp3 methylated_atom.
  If reference_atoms_tup is four atoms, interpret as an
  sp3 group adjacent to the methylated_atom and return
  the best rotamer for an sp3 methylated_atom."""
  methyl_distance = methyl_distances[on]
  methylated_atom_obj = locate_atom_by_name(residue, methylated_atom)
  methylated_atom_position = matrix.col(methylated_atom_obj.xyz)
  reference_atoms = [
    locate_atom_by_name(residue, atom) for atom in reference_atoms_tup]
  reference_atom_positions = [matrix.col(atom.xyz) for atom in reference_atoms]
  # first reference atom should be adjacent to methylated_atom
  # and others should be adjacent to it
  axis = methylated_atom_position - reference_atom_positions[0]
  all_directions = []
  if len(reference_atoms) == 2: # ref atom 1 sp3, pick 3 rotamers
    # (use this option when not all ref atom positions available)
    direction = reference_atom_positions[0] - reference_atom_positions[1]
    all_directions.append(direction)
    all_directions.append(direction.rotate_around_origin(
        axis=axis, angle=120, deg=True))
    all_directions.append(direction.rotate_around_origin(
        axis=axis, angle=240, deg=True))
  elif len(reference_atoms) == 3: # ref atom 1 sp2, pick 6 rotamers
    init_directions = [
      reference_atom_positions[0] - reference_atom_positions[1],
      reference_atom_positions[0] - reference_atom_positions[2]]
    for i in init_directions:
      all_directions.append(i)
      all_directions.append(i.rotate_around_origin(
        axis=axis, angle=120, deg=True))
      all_directions.append(i.rotate_around_origin(
        axis=axis, angle=240, deg=True))
  elif len(reference_atoms) == 4: # ref atom 1 sp3, pick 3 rotamers
    # (prefer this option when all three ref atom positions known)
    all_directions = [
      reference_atom_positions[0] - reference_atom_positions[1],
      reference_atom_positions[0] - reference_atom_positions[2],
      reference_atom_positions[0] - reference_atom_positions[3]]
  methyl_positions = [
    methyl_distance * d + methylated_atom_position for d in all_directions]
  new_atoms = [
    methylated_atom_obj.detached_copy() for i in xrange(len(methyl_positions))]
  for i, atom in enumerate(new_atoms):
    atom.set_name(name)
    atom.set_element("C")
    atom.set_xyz(tuple(methyl_positions[i]))
    residue.atom_groups()[0].append_atom(atom)
  return residue

def sp2_methylate(residue, methylated_atom, reference_atoms_tup, on="C", name=" C  "):
  """place a new methyl group on the given residue at the
  standard methyl bond distance from the methylated_atom.
  If reference_atoms_tup is three atoms, the direction of
  methylation will be (methylated - avg(reference)).
  If reference_atoms_tup is two atoms, the direction of
  methylation will be (ref[0] - ref[1])."""
  methyl_distance = methyl_distances[on]
  methylated_atom_obj = locate_atom_by_name(residue, methylated_atom)
  methylated_atom_position = matrix.col(methylated_atom_obj.xyz)
  reference_atoms = [
    locate_atom_by_name(residue, atom) for atom in reference_atoms_tup]
  if len(reference_atoms) == 3:
    reference_position = average_position(*reference_atoms)
    all_directions = [(methylated_atom_position - reference_position).normalize()]
  elif len(reference_atoms) == 2:
    reference_positions = [matrix.col(atom.xyz) for atom in reference_atoms]
    all_directions = [(reference_positions[0] - reference_positions[1]).normalize()]
  elif len(reference_atoms) == 4: # same as above but two options
    reference_positions = [
      [matrix.col(atom.xyz) for atom in reference_atoms[0:2]],
      [matrix.col(atom.xyz) for atom in reference_atoms[2:4]]]
    all_directions = [
      (reference_positions[i][0] - reference_positions[i][1]).normalize()
      for i in (0, 1)]
  else:
    raise Sorry("Unrecognized number of reference atoms supplied to sp2_methylate")
  methyl_positions = [
    methyl_distance * d + methylated_atom_position for d in all_directions]
  new_atoms = [
    methylated_atom_obj.detached_copy() for i in xrange(len(methyl_positions))]
  for i, atom in enumerate(new_atoms):
    atom.set_name(name)
    atom.set_element("C")
    atom.set_xyz(tuple(methyl_positions[i]))
    residue.atom_groups()[0].append_atom(atom)
  return residue

def sp2_dimethylate(residue, methylated_atom, reference_atoms_tup, on="N", names=[]):
  """place two new methyl groups on the given residue at the
  standard methyl bond distance from the methylated_atom, placed
  along the vectors between the 2nd and 1st reference atoms and
  the 3rd and 1st reference atoms extending from the methylated
  atom."""
  methyl_distance = methyl_distances[on]
  methylated_atom_obj = locate_atom_by_name(residue, methylated_atom)
  methylated_atom_position = matrix.col(methylated_atom_obj.xyz)
  reference_atoms = [
    locate_atom_by_name(residue, atom) for atom in reference_atoms_tup]
  reference_positions = [matrix.col(atom.xyz) for atom in reference_atoms]
  dir1 = (reference_positions[0] - reference_positions[1]).normalize()
  dir2 = (reference_positions[0] - reference_positions[2]).normalize()
  pos1 = methyl_distance * dir1 + methylated_atom_position
  pos2 = methyl_distance * dir2 + methylated_atom_position
  new1 = methylated_atom_obj.detached_copy()
  new1.set_name(names[0])
  new1.set_element("C")
  new1.set_xyz(tuple(pos1))
  residue.atom_groups()[0].append_atom(new1)
  new2 = methylated_atom_obj.detached_copy()
  new2.set_name(names[1])
  new2.set_element("C")
  new2.set_xyz(tuple(pos2))
  residue.atom_groups()[0].append_atom(new2)
  return residue

def place_modified(modified, unmodified, ref_residue):
  """find the transformation that places the umodified residue at the
  position of the reference residue in the model, and apply that
  transformation to the modified residue"""
  # use mmtbx/math/superpose.py. It also reports the transformation.
  raise NotImplemented

def resi_from_selected_atoms(ref_residue, atom_names):
  """copy ref_residue and produce a copy containing only the atoms
  specified. This can be used for a superposition of the part of a
  residue between the backbone and a PTM, excluding a flexible part."""
  # deep copy ref_residue and trim every atom not in atom_names
  # return the new residue and use this in modify_lambdas
  raise NotImplemented

def densities_and_ratio(mapdata_ref, mapdata_new, frac_matrix, residue,
  atoms_ref=[], atoms_new=[], mid_atoms_pairs=[], far_atoms_pairs=[]):
  """mapdata1 is expected to be an EM map and mapdata2 is probably a difference
  map in which to look for evidence of new atom positions. atoms_ref should be the
  reference atom objects and atoms_new should be the prospective newly placed atoms,
  organized as follows:
  atoms_example = (atom1, atom2, ...) <-- e.g. 3 for an acetyl group
  residue is the whole object we'll use to look up atom positions
  mid_atoms_pairs = [(atom_ref, atom_new), ...] <-- reference and new atom names
    so we can get the midpoint along a bond that should exist
  far_atoms_pairs = [(atom_near, atom_far), ...] <-- near and far atoms along
    a trajectory we'll go half that distance further along
  NOTE: once we look up positions, each of the atoms will produce a tuple of
    positions, and we'll have to test densities at each and take the best one.
    This does add another layer of grouping!

  To sum up, given atom positions represented below:
  (atoms_ref): possibly multiple atoms
  (atoms_new): possibly multiple atoms
  [(mid1, mid2)...]: pairs of atoms from which to take the mid position
  [(far1, far2)...]: pairs of atoms from which to extrapolate a new position

  Return the densities at the following positions:
  (1): average density at atoms_ref positions
  (mid): average density between each pair of positions in mid_atoms_pairs.
    we hope to find nonnegative density at this position!
  (2): average density at atoms_new positions
  (far): average density extrapolated from each pair of positions in
    far_atoms_pairs as follows:
    take far2 - far1 vector, divide by half, add to far2 position
    we hope *not* to find a lot of density at this position!
  and also return the density ratio (2)/(1).
  """
  # this is going to take a handful of helper functions
  def get_positions(list_of_atom_names):
    return [map(lambda a: a.xyz, locate_atoms_by_name(residue, atom_name))
      for atom_name in list_of_atom_names]
  def get_midpoint(position_pair):
    return average_position(*position_pair)
  def get_extrapolated_far_point(position_pair):
    p0, p1 = [matrix.col(position_pair[i]) for i in (0,1)]
    return 0.5*(p1-p0) + p1
  def get_nonnegative_densities(positions, mapdata):
    # to avoid ratios involving negative densities looking favorable
    return [max(0.00000001, max(map(lambda xyz: get_density_at_position(
      xyz, frac_matrix, mapdata), list(xyzs)))) for xyzs in positions]
  def get_paired_nonnegative_densities(list_of_pairs_of_atom_tuples, mapdata, find="mid"):
    """
    if position == "mid":
      for pair in list_of_pairs_of_atom_tuples:
        calculate all the densities halfway between the pair of atoms (at all conformers)
        keep the best one
      calculate the average density representing where we expect a covalent bond

    elif position == "far":
      for pair in list_of_pairs_of_atom_tuples:
        calculate new positions halfway further along the trajectory between the atoms
        calculate densiites at these new positions where we expect density to have
          dropped off relative to the positions of proposed atoms
        keep the best one
      calculate the average density representing beyond the edge of the PTM
    """
    pairs_densities = []
    # TODO: rework this section as a list comprehension instead of a loop
    for pair in list_of_pairs_of_atom_tuples:
      positions = get_positions(pair)
      # e.g. [[a1, a2], [b1, b2, b3]] <-- a's and b's are positions at conformers
      paired_positions = [(a, b) for a in positions[0] for b in positions[1]]
      # e.g. [(a1, b1), (a1, b2), (a1, b3), (a2, b1), (a2, b2), (a2, b3)]
      if find == "mid":
        new_positions = [[get_midpoint(pair) for pair in paired_positions]]
      elif find == "far":
        new_positions = [[get_extrapolated_far_point(pair) for pair in paired_positions]]
      else:
        raise Sorry, ("Unrecognized type of position calculation %s" % find)
      pairs_densities.extend(get_nonnegative_densities(new_positions, mapdata))
    return sum(pairs_densities)/len(pairs_densities)

  # TODO: error handling each time we try to get the nonnegative densiites
  # from a list that might be empty (right now get a ValueError)

  # for the reference and placed atom positions:
  atoms_ref_xyz = get_positions(atoms_ref)
  atoms_new_xyz = get_positions(atoms_new)
  # densities ref: density at the reference positions in the original map
  densities_ref = get_nonnegative_densities(atoms_ref_xyz, mapdata_ref)
  # densities_atoms_new_in_map_ref: density at the new positions in the original map
  densities_atoms_new_in_map_ref = get_nonnegative_densities(atoms_new_xyz, mapdata_new)
  # densities_new: **difference density** at the new positions
  densities_atoms_new_in_diff_map = get_nonnegative_densities(atoms_new_xyz, mapdata_new)
  d_ref = sum(densities_ref)/len(densities_ref)
  d_new_in_ref = sum(densities_atoms_new_in_map_ref)/len(densities_atoms_new_in_map_ref)
  d_new_diff = sum(densities_atoms_new_in_diff_map)/len(densities_atoms_new_in_diff_map)
  # for the mid and far positions:
  d_mid = get_paired_nonnegative_densities(mid_atoms_pairs, mapdata_ref, find="mid")
  d_far = get_paired_nonnegative_densities(far_atoms_pairs, mapdata_ref, find="far")
  # ratio of the placed to reference atom densities:
  # FIXME TODO do we actually want to take the ratio with the difference density or no????
  ratio = d_new_diff/d_ref
  return (d_ref, d_new_diff, d_mid, d_new_in_ref, d_far, ratio)

def fractionalize(coords, ucell_params):
  """return ucell_params in fractional coordiantes"""
  return tuple([coords[i]/ucell_params[i] for i in xrange(3)])

def get_cc_of_residue_to_map(residue, frac_matrix, ucell_params, mapdata, fcalc_map):
  """generic function for computing a correlation coefficient between a set of model
     atom positions and an EM map"""
  positions = flex.vec3_double([a.xyz for a in residue.atom_groups()[0].atoms()])
  pos_frac = flex.vec3_double([fractionalize(p, ucell_params) for p in positions])
  # get densities from the EM map at the model atom positions
  values_em = flex.double([get_density_at_position(p, frac_matrix, mapdata)
    for p in positions])
  # get matching densities in the F_calc map
  values_fc = flex.double([fcalc_map.eight_point_interpolation(p) for p in pos_frac])
  cc = flex.linear_correlation(x=values_em, y=values_fc).coefficient()
  return cc

def prune_confs(mapdata, frac_matrix, residue):
  """remove duplicate atoms from a hierarchical residue model"""
  atom_groups = residue.atom_groups()
  for ag in atom_groups:
    atoms = ag.atoms()
    atoms_grouped = {name:[a for a in atoms if a.name.strip() == name]
      for name in [a.name.strip() for a in atoms]}
    atoms_duplicates = {name:values for name, values in atoms_grouped.iteritems()
      if len(values) > 1}
    for name, values in atoms_duplicates.iteritems():
      best_density = -100
      best_atom = None
      for a in values:
        test_density = get_density_at_position(a.xyz, frac_matrix, mapdata)
        if test_density > best_density:
          if best_atom is not None:
            ag.remove_atom(best_atom)
          best_density = test_density
          best_atom = a
        else:
          ag.remove_atom(a)

def prune_atoms(residue, keep_atom_list):
  """remove atoms not in keep_atom_list from residue"""
  atom_groups = residue.atom_groups()
  for ag in atom_groups:
    atoms = ag.atoms()
    atom_names = [a.name.strip() for a in atoms]
    for i, aname in enumerate(atom_names):
      if aname not in keep_atom_list:
        ag.remove_atom(atoms[i])

def rename(residue, resname):
  for ag in residue.atom_groups():
    ag.resname = resname

"""PTM categories grouped by position of the modification(s) on the residue for
adjustment of the score to match how strong we expect the density to be out there"""
PTM_score_scaling_by_position = {
  "nucleotides": {
    "sugar_ring":1, # modifications directly on a sugar ring atom
    "sugar_ter":1.5, # modifications directly on a sugar terminal atom, e.g. OMC
    "1-3_ring":1.2, # modifications 1-3 atoms from the sugar, directly on a ring, e.g. G7M
    "1-3_ter":2, # modifications 1-3 atoms from the sugar, terminal atom
    "4-6_ring":1.5, # modifications 4-6 atoms from the sugar, directly on a ring, e.g. 1MG
    "4-6_ter":3, # modifications 4-6 atoms from the sugar, terminal atom, e.g. 2MG
  }
}

"""PTM: dict of known posttranslational modifications

top level: nucleic acids or protein residues
second level: modifications at the selected atom
third level: the full posttranslational modification name, molecular
structure (if it's a solitary fragment), structure of the complete
modified residue (for some we will substitute this for the model),
and methods for placing it in the map and computing its fit score.
Note, lambdas may modify the residues passed to them (but not the
model or map) -- they should always be passed copies using the
residue_group.detached_copy() method.
"""
# See G7M for a fully annotated example
# TODO: get a document together with diagrams of all of these, with atom labels and
# annotation of which atoms are the reference, new, mid and far atom collections
PTM_lookup = {
  "A":{
    "unmodified":{
      "atoms":set(["N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9",
        "C1'", "C2'", "O2'", "C3'", "O3'", "C4'", "O4'", "C5'", "O5'", "P", "OP1", "OP2", "OP3"])
    },
    "MA6":{ # https://www.rcsb.org/ligand/MA6
      "name":"m6m6A (N6-Dimethyladenosine)",
      "goto_atom":"N6",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp2_dimethylate(residue, "N6", ("C6", "N1", "C5"), on="N", names=[" C9 ", " C10"]),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("N6",),
          atoms_new=("C9", "C10"),
          mid_atoms_pairs=[("N6", "C9"), ("N6", "C10")],
          far_atoms_pairs=[("N6", "C9"), ("N6", "C10")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["4-6_ter"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["A"]["unmodified"]["atoms"])
    },
    "6MZ":{ # https://www.rcsb.org/ligand/6MZ
      "name":"m6A (N6-Methyladenosine)",
      "goto_atom":"N6",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp3_rot_methylate(residue, "N6", ("C6", "N1", "C5"), on="C", name=" CZ "),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("N6",),
          atoms_new=("CZ",),
          mid_atoms_pairs=[("N6", "CZ")],
          far_atoms_pairs=[("N6", "CZ")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["4-6_ter"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["A"]["unmodified"]["atoms"])
    },
    "2MA":{ # https://www.rcsb.org/ligand/2MA
      "name":"m2A (2-Methyladenosine)",
      "goto_atom":"C2",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp2_methylate(residue, "C2", ("C2", "N3", "N1"), on="C", name=" CM2"),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("C2",),
          atoms_new=("CM2",),
          mid_atoms_pairs=[("C2", "CM2")],
          far_atoms_pairs=[("C2", "CM2")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["4-6_ring"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["A"]["unmodified"]["atoms"])
    },
    "A2M":{ # https://www.rcsb.org/ligand/A2M
      "name":"m(2'O)A (2'O-Methyladenosine)",
      "goto_atom":"O2'",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp3_rot_methylate(residue, "O2'", ("C2'", "C1'"), on="O", name=" CM'"),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("O2'",),
          atoms_new=("CM'",),
          mid_atoms_pairs=[("O2'", "CM'")],
          far_atoms_pairs=[("O2'", "CM'")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["sugar_ter"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["A"]["unmodified"]["atoms"])
    }
  },
  "U":{
    "unmodified":{
      "atoms":set(["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6",
        "C1'", "C2'", "O2'", "C3'", "O3'", "C4'", "O4'", "C5'", "O5'", "P", "OP1", "OP2", "OP3"])
    },
    # "PSU":{ # example of a modification to existing atoms
    #   "name":"PSU (pseudouridine)",
    #   "model":None,
    #   "model_on_struct":(,), # model of the complete pseudouridine
    #   "modify_lambda":lambda residue:residue, # not implemented (FIXME)
    #   # some lambda function that will align the unmodified uridine
    #   # to the one in the model and then apply that transformation
    #   # to this modified one to score. This is a lambda so that we can
    #   # align to the relevant (usually, modified) part of the residue.
    #   # *where possible, this should just be an extension along a vector!*
    #   "ratio_lambda":lambda model, mapdata, frac_matrix, fitted_modded:(0, 0, 0)
    #   # there is no density threshold since no new atoms are present in
    #   # this particular case. Ignore this for now.
    #   "scale_ratio_lambda":lambda model, fitted_modded, ratio:score
    #   # for this one we want to look for a hydrogen bonding partner in the model
    # },
    "UR3":{ # https://www.rcsb.org/ligand/UR3
      "name":"m3U (N3-Methyluridine)",
      "goto_atom":"N3",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp2_methylate(residue, "N3", ("C4", "N3", "C2"), on="N", name=" C3U"),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("N3",),
          atoms_new=("C3U",),
          mid_atoms_pairs=[("N3", "C3U")],
          far_atoms_pairs=[("N3", "C3U")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["1-3_ring"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["U"]["unmodified"]["atoms"])
    },
    "5MU":{ # https://www.rcsb.org/ligand/5MU
      "name":"m5U (5-Methyluridine)",
      "goto_atom":"C5",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp2_methylate(residue, "C5", ("C4", "C5", "C6"), on="C", name=" C5M"),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("C5",),
          atoms_new=("C5M",),
          mid_atoms_pairs=[("C5", "C5M")],
          far_atoms_pairs=[("C5", "C5M")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["1-3_ring"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["U"]["unmodified"]["atoms"])
    },
    "OMU":{ # https://www.rcsb.org/ligand/OMU
      "name":"m(2'O)U (2'O-Methyluridine)",
      "goto_atom":"O2'",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp3_rot_methylate(residue, "O2'", ("C2'", "C1'"), on="O", name=" CM2"),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("O2'",),
          atoms_new=("CM2",),
          mid_atoms_pairs=[("O2'", "CM2")],
          far_atoms_pairs=[("O2'", "CM2")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["sugar_ter"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["U"]["unmodified"]["atoms"])
    }
  },
  "C":{
    "unmodified":{
      "atoms":set(["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6",
        "C1'", "C2'", "O2'", "C3'", "O3'", "C4'", "O4'", "C5'", "O5'", "P", "OP1", "OP2", "OP3"])
    },
    "5MC":{ # https://www.rcsb.org/ligand/5MC
      "name":"m5C (5-Methylcytidine)",
      "goto_atom":"C5",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp2_methylate(residue, "C5", ("C4", "C5", "C6"), on="C", name=" CM5"),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("C5",),
          atoms_new=("CM5",),
          mid_atoms_pairs=[("C5", "CM5")],
          far_atoms_pairs=[("C5", "CM5")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["1-3_ring"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["C"]["unmodified"]["atoms"])
    },
    "OMC":{ # https://www.rcsb.org/ligand/OMC
      "name":"m(2'O)C (2'O-Methylcytidine)",
      "goto_atom":"O2'",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp3_rot_methylate(residue, "O2'", ("C2'", "C1'"), on="O", name=" CM2"),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("O2'",),
          atoms_new=("CM2",),
          mid_atoms_pairs=[("O2'", "CM2")],
          far_atoms_pairs=[("O2'", "CM2")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["sugar_ter"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["C"]["unmodified"]["atoms"])
    },
    "4OC":{ # https://www.rcsb.org/ligand/4OC
      "name":"m4C (N4,O2'-Dimethylcytidine)",
      "goto_atom":"N4",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp3_rot_methylate(
          sp3_rot_methylate(residue, "N4", ("C4", "N3", "C5"), on="N", name=" CM4"),
          "O2'", ("C2'", "C1'",), on="O", name=" CM2"),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("N4",),
          atoms_new=("CM4", "CM2"),
          mid_atoms_pairs=[("O2'", "CM2"), ("N4", "CM4")],
          far_atoms_pairs=[("O2'", "CM2"), ("N4", "CM4")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["1-3_ter"], # dual position, compromise
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["C"]["unmodified"]["atoms"])
    }
  },
  "G":{
    "unmodified":{
      "atoms":set(["N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", "N9",
        "C1'", "C2'", "O2'", "C3'", "O3'", "C4'", "O4'", "C5'", "O5'", "P", "OP1", "OP2", "OP3"])
    }, # phenix-style structure of the unmodified residue
    "G7M":{ # https://www.rcsb.org/ligand/G7M, example of an addition of new atoms only (ignoring H)
      "name":"m7G (N7-Methylguanosine)",
      "goto_atom":"N7",
      "model":None, # phenix-style hierarchical model of the modification,
      "model_on_struct":None, # model including the full nucleotide
      "modify_lambda":lambda residue:\
        sp2_methylate(residue, "N7", ("C5", "N7", "C8"), on="N", name=" CM7"),
      # some lambda function that will align the unmodified uridine
      # to the one in the model and then apply that transformation
      # to this modified one to score. This is a lambda so that we can
      # align to the relevant (usually, modified) part of the residue.
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("N7",),
          atoms_new=("CM7",),
          mid_atoms_pairs=[("N7", "CM7")],
          far_atoms_pairs=[("N7", "CM7")]),
        # returns tuple (density_at_atom1, density_at_atom2, ratio2:1)
      # *store the densities -- let the user analyze a distribution*
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["1-3_ring"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["G"]["unmodified"]["atoms"])
      # some lambda function to decide whether the PTM at this site looks
      # reasonably well evidenced by the map, as reported by this score,
      # hopefully usually just based on the ratio in most cases
    },
    "2MG":{ # https://www.rcsb.org/ligand/2MG
      "name":"m2G (N2-Methylguanosine)",
      "goto_atom":"N2",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp2_methylate(residue, "N2", ("C2", "N1", "C2", "N3"), on="N", name=" CM2"),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("N2",),
          atoms_new=("CM2",),
          mid_atoms_pairs=[("N2", "CM2")],
          far_atoms_pairs=[("N2", "CM2")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["4-6_ter"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["G"]["unmodified"]["atoms"])
    },
    "M2G":{ # https://www.rcsb.org/ligand/M2G
      "name":"m2m2G (N2-Dimethylguanosine)",
      "goto_atom":"N2",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp2_dimethylate(residue, "N2", ("C2", "N1", "N3"), on="N", names=[" CM1", " CM2"]),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("N2",),
          atoms_new=("CM1", "CM2"),
          mid_atoms_pairs=[("N2", "CM1"), ("N2", "CM2")],
          far_atoms_pairs=[("N2", "CM1"), ("N2", "CM2")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["4-6_ter"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["G"]["unmodified"]["atoms"])
    },
    "1MG":{ # https://www.rcsb.org/ligand/1MG
      "name":"m1G (N1-Methylguanosine)",
      "goto_atom":"N1",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp2_methylate(residue, "N1", ("C2", "N1", "C6"), on="N", name=" CM1"),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("N1",),
          atoms_new=("CM1",),
          mid_atoms_pairs=[("N1", "CM1")],
          far_atoms_pairs=[("N1", "CM1")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["4-6_ring"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["G"]["unmodified"]["atoms"])
    },
    "OMG":{ # https://www.rcsb.org/ligand/OMG
      "name":"m(2'O)G (2'O-Methylguanosine)",
      "goto_atom":"O2'",
      "model":None,
      "model_on_struct":None,
      "modify_lambda":lambda residue:\
        sp3_rot_methylate(residue, "O2'", ("C2'", "C1'",), on="O", name=" CM2"),
      "ratio_lambda":lambda model, mapdata, diffmapdata, frac_matrix, fitted_modded:\
        densities_and_ratio(
          mapdata, diffmapdata, frac_matrix, fitted_modded,
          atoms_ref=("O2'",),
          atoms_new=("CM2",),
          mid_atoms_pairs=[("O2'", "CM2")],
          far_atoms_pairs=[("O2'", "CM2")]),
      "scale_ratio_lambda":lambda model, fitted_modded, d1, d2, ratio:\
        ratio*PTM_score_scaling_by_position["nucleotides"]["sugar_ter"],
      "prune_lambda":lambda model:\
        prune_atoms(model, PTM_lookup["G"]["unmodified"]["atoms"])
    },
  }
}

PTM_reverse_lookup = {
  modname:resname for (resname,resdict) in PTM_lookup.iteritems() for modname in resdict
}
