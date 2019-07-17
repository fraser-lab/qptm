from __future__ import division

# qPTxM: Quantifying Posttranscriptional Modifications
#
# This program accepts as input an EM map (.map, .mrc or .ccp4) and a model (.pdb or
# .cif) and produces the locations and identities of posttranscriptional modifications
# that may be present in the model based on density present in the map. It is based on
# cctbx and Phenix tools.
#
# citation: Stojkovic V, Myasnikov AG, Young ID, Frost A, Fraser JS, Fujimori DG.
# High-resolution cryo-electron microscopy structure of the Escherichia coli 50S
# subunit and validation of nucleotide modifications.
#
# Preprint available on bioRxiv (https://www.biorxiv.org/content/10.1101/695429v1)
# and pending peer review.

import libtbx.phil
from libtbx.utils import Sorry, Usage
from qptm_core import LookForPTMs
from ptm_util import prune_ptms

master_phil_str = """
model_file = None
  .type = path
  .help = "Path to the model to be analyzed, in pdb or mmcif format"
map_file = None
  .type = path
  .help = "Path to the map into which the model was built, in mrc or ccp4 format"
difference_map_file = None
  .type = path
  .help = "Path to a difference map to use directly (instead of calculating one)."
calculated_map_file = None
  .type = path
  .help = "Path to a calculated map to use directly (instead of calculating one)."
experiment = *electron xray
  .type = choice
  .help = "Type of experiment (dictating the scattering table to use for the"
  .help = "calculated map)"
d_min = 3
  .type = float
  .help = "Resolution to use for calculation of CC of residues with the EM map"
cc_threshold = 0.6
  .type = float
  .help = "Threshold to determine the minimum correlation coefficient between the map"
  .help = "and residue to search for modifications on that residue."
ratio_d_far_d_new_in_ref = 1.5
  .type = float
  .help = "Ratio guiding how strong the density just past the proposed atom positions"
  .help = "is permitted to be relative to the density at the proposed atom positions."
ratio_d_far_d_mid = 1
  .type = float
  .help = "Ratio guiding how strong the density just past the proposed atom positions"
  .help = "is permitted to be relative to the density midway between the proposed"
  .help = "atoms and the reference atoms."
ratio_d_ref_d_new_in_ref = 3
  .type = float
  .help = "Ratio guiding how strong the experimental map density at the reference"
  .help = "positions is permitted to be relative to the difference map density at the"
  .help = "proposed atom positions."
reference_densities_fraction = 0.5
  .type = float
  .help = "Fraction of modifications to keep, selecting those with the strongest"
  .help = "densities at the reference atom positions."
difference_densities_fraction = 0.1
  .type = float
  .help = "Fraction of modifications to keep, selecting those with the strongest"
  .help = "difference densities (unmodeled densities)."
score_threshold = 0.5
  .type = float
  .help = "Minimum score for a modification to be included in the suggested PTM list."
selected_ptms = None
  .type = path
  .help = "If provided, interpret this file as output from this program that has been"
  .help = "curated by the user, and apply any modifications listed."
synthetic_data = False
  .type = bool
  .help = "Randomly modify 10%% of the recognized residues to generate synthetic data."
adjust_filters_only = False
  .type = bool
  .help = "Don't rerun, just look for existing ptms.out and rejected_ptms.out and redo"
  .help = "the step to try out different filters (cc_threshold, score_threshold,"
  .help = "ratio_d_far_d_new_in_ref, ratio_d_far_d_mid, ratio_d_ref_d_new_in_ref,"
  .help = "reference_densities_fraction, difference_densities_fraction)"
plot = True
  .type = bool
  .help = "Show plots upon completion."
"""

master_phil = libtbx.phil.parse(master_phil_str)

helpstr = """
Hello! qPTxM is a tool for quantifying posttranscriptional modifications. To
run, please supply a model (pdb or mmcif) and map (mrc or ccp4) and specify
an estimated median resolution:

qptxm.py model.pdb map.mrc d_min=3

Any modifications already present (if recognized) will be stripped off and
searched for anew, and you will be shown some plots describing which possible
modifications look promising and which were rejected. You can also load a copy
of your model (ptms.pdb) with all proposed modifications modeled, and step
through these positions with a custom script for Coot:

coot goto_ptms.py ptms.pdb

You may also like to load pruned.pdb to compare against the model with no
modifications modeled. Under the menu option Custom, you will find the
option to call up a panel of buttons to take you directly to each modeled
modification.

Once you have decided on which proposed modifications to keep, remove any
you want to reject from the file ptms.out and rerun qptxm.py with the
additional argument selected_ptms=ptms.out (or whatever you may want to
rename this). It will produce another model with just those modifications.

You may also find it useful to run with adjust_filters_only=True in order
to test how many modifications are suggested if you adjust any of the
optional parameters, like the minimum correlation coefficient of the model
to the map (see below). This way is much faster but does not produce a new
model and goto_ptms.py script. It also won't be able to adjust results based
on an updated resolution estimate. You can also turn off plotting with
plot=False or test choices of parameters on synthetic data (generated from
your model but a calculated, noise-free map of the same dimensions) and see
how accurately qptxm finds a randomly-modified 10%% of positions in that
map by setting synthetic_data=True.

Available parameters:
""" + master_phil_str

def validate_params(params):
  import os
  if not params.model_file:
    raise Sorry("A molecular model (.pdb or .cif) is required.")
  if params.experiment == "electron" and not params.map_file:
    raise Sorry("A map file (.map, .mrc or .ccp4) is required for EM experiments.")
  if params.experiment == "xray" and not (params.difference_map_file and params.calculated_map_file):
    raise Sorry("Calculated and difference map files (.map, .mrc or .ccp4) are required for X-ray experiments.")
  if params.selected_ptms is not None and not os.path.exists(params.selected_ptms):
    raise Sorry("Could not locate the file provided: %s" % params.selected_ptms)
  if params.reference_densities_fraction < 0 or params.reference_densities_fraction > 1:
    raise Sorry("Please select a fraction between 0 and 1 for reference_densities_fraction.")
  if params.difference_densities_fraction < 0 or params.difference_densities_fraction > 1:
    raise Sorry("Please select a fraction between 0 and 1 for difference_densities_fraction.")

def run(args):
  if not args or "-h" in args or "--help" in args:
    raise Usage(helpstr)
  from iotbx import file_reader
  from goto_ptms_gen import gen_from_ptms
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    )
  if cmdline.unused_args:
    print "\nEncountered unrecognized parameters:", str(cmdline.unused_args), "\n"
    return
  params = cmdline.work.extract()
  if params.adjust_filters_only:
    # only redo the filtering step, rewriting all_tested_ptms.out, ptms.out,
    # params.out and the plots
    from filter_util import import_ptms, import_synthetic_ptms, apply_filters, write_ptms_from_flex_arrays
    imported_ptms = import_ptms("all_tested_ptms.out")
    synthetic_ptms = import_synthetic_ptms("synthetic_ptms.out") if params.synthetic_data else None
    keep_selection, accepted, all_tested_ptms, log = apply_filters(
      imported_ptms,
      imported_synthetic_ptms=synthetic_ptms,
      cc_threshold=params.cc_threshold,
      ref_frac=params.reference_densities_fraction,
      dif_frac=params.difference_densities_fraction,
      ratio_d_far_d_new_in_ref=params.ratio_d_far_d_new_in_ref,
      ratio_d_far_d_mid=params.ratio_d_far_d_mid,
      ratio_d_ref_d_new_in_ref=params.ratio_d_ref_d_new_in_ref,
      score_threshold=params.score_threshold)
    write_ptms_from_flex_arrays(accepted, "ptms.out")
    write_ptms_from_flex_arrays(all_tested_ptms, "all_tested_ptms.out")
    if params.plot:
      from plot_util import plot_densities_from_flatfile
      plot_densities_from_flatfile("all_tested_ptms.out", "ptms.out", params.cc_threshold)
    with open("params.out", "wb") as outf:
      outf.write(cmdline.work.as_str())
    return
  validate_params(params)
  with open("params.out", "wb") as outf:
    outf.write(cmdline.work.as_str())
  # process the model
  model_in = file_reader.any_file(params.model_file, force_type="pdb")
  model_in.check_file_type("pdb")
  hier_model = model_in.file_object.construct_hierarchy()
  hier_model.remove_alt_confs(True) # only accommodates the major conformer
  hier_model.remove_hd()
  prune_ptms(hier_model, filename="pruned.pdb")
  pruned_pdb_in = cmdline.get_file("pruned.pdb").file_object
  # process the map(s)
  def get_map_file_object(map_path):
    if not map_path: return
    map_in = file_reader.any_file(map_path, force_type="ccp4_map")
    map_in.check_file_type("ccp4_map")
    return map_in.file_object
  look_for_ptms = LookForPTMs(pruned_pdb_in, hier_model,
    emmap=get_map_file_object(params.map_file),
    diff_map=get_map_file_object(params.difference_map_file),
    calc_map=get_map_file_object(params.calculated_map_file),
    params=params)
  # TODO: some way to check whether the model is actually on the map?
  # maybe just be prepared to throw exceptions if we ask for density
  # at a position and can't get any
  if params.selected_ptms is None:
    look_for_ptms.write_identified_ptms()
    look_for_ptms.write_all_tested_ptms()
  look_for_ptms.write_ccs()
  look_for_ptms.write_difference_map(filename="difference_map.ccp4")
  if params.selected_ptms is not None:
    look_for_ptms.write_modified_model(filename="modified.pdb")
    gen_from_ptms(ptms_flatfile=params.selected_ptms)
  else:
    look_for_ptms.write_modified_model(filename="ptms.pdb")
    gen_from_ptms(ptms_flatfile="ptms.out")
    if params.plot:
      from plot_util import plot_densities_from_flatfile
      plot_densities_from_flatfile("all_tested_ptms.out", "ptms.out", params.cc_threshold)

if __name__=="__main__":
  import sys
  run(sys.argv[1:])


# WISHLIST: accept multiple models, e.g. if a ribosome has a large and small subunit as separate
# pdbs. Should be able to loop through the hierarchical models. 

# WISHLIST: a method to inject new PTMs for DNA, RNA or protein into the dictionaries before
# launching the analysis. Could set a default lambda that aligns the whole residue.




