from __future__ import division
from scitbx.array_family import flex
from cctbx import maptbx
from cctbx.maptbx import non_crystallographic_eight_point_interpolation

def get_density_at_position(position, frac_matrix, mapdata):
  """generic function for computing the density at a single point"""
  # raise Sorry if the position doesn't appear to be within the bounds of the map.
  # this may mean the user supplied a model that was not centered in the map.
  density = non_crystallographic_eight_point_interpolation(
    map=mapdata,
    gridding_matrix=frac_matrix,
    site_cart=position)
  return density

def get_fcalc_map(pdb, symmetry, d_min, emmap, scatterer="electron"):
  """Fcalc map generation borrowed from EMRinger"""
  xrs = pdb.input.xray_structure_simple(crystal_symmetry=symmetry)
  xrs.scattering_type_registry(
    d_min=d_min,
    table=scatterer)
  fc = xrs.structure_factors(d_min=d_min).f_calc()
  cg = maptbx.crystal_gridding(
    unit_cell=symmetry.unit_cell(),
    space_group_info=symmetry.space_group_info(),
    pre_determined_n_real=emmap.data.all())
  fc_map = fc.fft_map(
    crystal_gridding=cg).apply_sigma_scaling().real_map_unpadded()
  return fc_map
