from __future__ import division
from cctbx.maptbx import non_crystallographic_eight_point_interpolation

def get_density_at_position(position, frac_matrix, emmap):
  """generic function for computing the density at a single point"""
  # raise Sorry if the position doesn't appear to be within the bounds of the map.
  # this may mean the user supplied a model that was not centered in the map.
  density = non_crystallographic_eight_point_interpolation(
    map=emmap,
    gridding_matrix=frac_matrix,
    site_cart=position)
  return density
