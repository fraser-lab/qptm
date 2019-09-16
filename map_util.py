from __future__ import division
from scitbx.array_family import flex
from cctbx import maptbx, uctbx, crystal
from cctbx.maptbx import non_crystallographic_eight_point_interpolation
from mmtbx.command_line.real_space_diff_map import scale_k1, scale_two_real_maps_in_fourier_space
from mmtbx.maps.fem import ccp4_map
from libtbx.utils import Sorry

def get_density_at_position(position, frac_matrix, mapdata):
  """generic function for computing the density at a single point"""
  # raise Sorry if the position doesn't appear to be within the bounds of the map.
  # this may mean the user supplied a model that was not centered in the map.
  try:
    density = non_crystallographic_eight_point_interpolation(
      map=mapdata,
      gridding_matrix=frac_matrix,
      site_cart=position)
  except RuntimeError as e:
    raise Sorry("Could not compute map density at the requested position. "+\
      "Please check that model is entirely within map bounds.")
  return density

def neutralize_scatterers(xrs):
  scatterers = xrs.scatterers()
  for sc in scatterers:
    neutralized_scatterer = filter(lambda x: x.isalpha(), sc.scattering_type)
    sc.scattering_type = neutralized_scatterer

def get_fcalc_map(model, symmetry, d_min, emmap, scatterer="electron"):
  """Fcalc map generation borrowed from EMRinger"""
  xrs = model.input.xray_structure_simple(crystal_symmetry=symmetry)
  if scatterer == "electron":
    # until ions are added to the electron scattering tables
    neutralize_scatterers(xrs)
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
  scale = len(fc_map)/flex.sum(flex.abs(fc_map))
  fc_map_scaled = fc_map * scale
  return fc_map_scaled

def get_diff_map(symmetry, fcalc_map, expt_map, d_min):
  """Compute a difference map between the calculated and experimental maps.
  Modified from mmtbx.command_line.real_space_diff_map class compdiff."""
  scale = scale_k1(x=expt_map, y=fcalc_map)
  scaled_fcalc_map = fcalc_map * scale # should already be scaled but okay
  origin = expt_map.origin()
  shifted_expt_map = expt_map.shift_origin()
  shifted_diff_map = flex.double(flex.grid(expt_map.all()))
  diff_map = flex.double(flex.grid(expt_map.origin(), expt_map.focus()))
  ucell_params = symmetry.unit_cell().parameters()
  boxes = maptbx.boxes_by_dimension(
    n_real = expt_map.all(),
    dim    = 30,
    abc    = ucell_params[:3])
  i_box = 0
  test_map_box_obs_all = None
  for start, end in zip(boxes.starts, boxes.ends):
    i_box += 1
    map_box_obs = maptbx.copy(shifted_expt_map, start, end)
    map_box_calc = maptbx.copy(scaled_fcalc_map, start, end)
    map_box_obs.reshape(flex.grid(map_box_obs.all()))
    map_box_calc.reshape(flex.grid(map_box_calc.all()))
    # abc = [ucell_params[i]/expt_map.all()[i] for i in range(3)]
    # ucb = uctbx.unit_cell(parameters=(
      # abc[0],abc[1],abc[2],ucell_params[3],ucell_params[4],ucell_params[5]))
    # cs = crystal.symmetry(unit_cell=ucb, space_group="P1")
    cs = symmetry
    diff_map_part = scale_two_real_maps_in_fourier_space(
        m1         = map_box_obs,
        m2         = map_box_calc,
        cs         = cs,
        d_min      = d_min,
        vector_map = True) # this means yes, use the phases from fcalc
    maptbx.set_box(
        map_data_from = diff_map_part,
        map_data_to   = shifted_diff_map,
        start         = start,
        end           = end)
  # skip this additional scaling step -- it explodes everything
  # sd = diff_map.sample_standard_deviation()
  # if(sd!=0):
    # diff_map = diff_map/sd
  # somehow put the contents of shifted_diff_map into diff_map but keep origin?
  for i in xrange(len(shifted_diff_map)):
    diff_map[i] = shifted_diff_map[i]
    # FIXME SUPER SLOW but it should work for now
  return diff_map

def write_ccp4_map(mapdata, symmetry, filename):
  cg = maptbx.crystal_gridding(
    unit_cell=symmetry.unit_cell(),
    space_group_info=symmetry.space_group_info(),
    pre_determined_n_real=mapdata.all())
  return ccp4_map(cg, filename, map_data=mapdata)

