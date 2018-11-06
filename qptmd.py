from __future__ import division

# QPTMD: Quantitative Posttranslational Modification Detection
#
# This program accepts as input an EM map (.map, .mrc or .ccp4) and a model (.pdb or
# .cif) and produces the locations and identities of posttranslational modifications
# that may be present in the model based on density present in the map. It is based on
# cctbx and Phenix tools.

import libtbx.phil
# import numpy
# from libtbx import easy_pickle
# from libtbx import easy_mp
# from libtbx.str_utils import make_header
from libtbx.utils import Sorry, Usage
# from libtbx import adopt_init_args, Auto
# from cStringIO import StringIO
# import time
# import os
# import sys
# import sqlite3
from qptmd_core import LookForPTMs

master_phil_str = """
model_file = None
  .type = path
  .help = "Path to the model to be analyzed, in pdb or mmcif format"
map_file = None
  .type = path
  .help = "Path to the map into which the model was built, in mrc or ccp4 format"
d_min = 3
  .type = float
  .help = "Resolution to use for calcualtion of CC of residues with the EM map"
score_threshold = None
  .type = float
  .help = "If provided, use this threshold to determine which modifications to apply"
  .help = "to the model (use the best-scoring modification if multiple qualify)."
selected_ptms = None
  .type = path
  .help = "If provided, interpret this file as output from this program that has been"
  .help = "curated by the user, and apply any modifications listed."
"""

master_phil = libtbx.phil.parse(master_phil_str)

helpstr = """
No help for you. (FIXME)

Available parameters:
""" + master_phil_str

def validate_params(params):
  if not params.model_file:
    raise Sorry("A molecular model (.pdb or .cif) is required.")
  if not params.map_file:
    raise Sorry("A map file (.map, .mrc or .ccp4) is required.")
  if params.selected_ptms is not None and not os.path.exists(params.selected_ptms):
    raise Sorry("Could not locate the file provided: %s" % params.selected_ptms)
  # if params.selected_ptms and params.score_threshold:
  #   raise Sorry("Cannot apply the selected modifications and also a score "+\
  #     "threshold. Please select only one of these options to apply.")

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
  params = cmdline.work.extract()
  validate_params(params)
  model_in = file_reader.any_file(params.model_file, force_type="pdb")
  model_in.check_file_type("pdb")
  hier_model = model_in.file_object.construct_hierarchy()
  hier_model.remove_alt_confs(True) # only accommodates the major conformer
  map_in = file_reader.any_file(params.map_file, force_type="ccp4_map")
  map_in.check_file_type("ccp4_map")
  emmap = map_in.file_object
  # TODO: some way to check whether the model is actually on the map?
  # maybe just be prepared to throw exceptions if we ask for density
  # at a position and can't get any
  pdb_in = cmdline.get_file(params.model_file).file_object
  look_for_ptms = LookForPTMs(pdb_in, hier_model, emmap, params=params)
  look_for_ptms.write_identified_ptms()
  look_for_ptms.write_ccs()
  gen_from_ptms(ptms_flatfile="ptms.out")
  from plot_util import plot_densities_from_flatfile
  plot_densities_from_flatfile("ptms.out")

if __name__=="__main__":
  import sys
  run(sys.argv[1:])


# WISHLIST: accept multiple models, e.g. if a ribosome has a large and small subunit as separate
# pdbs. Should be able to loop through the hierarchical models. 

# WISHLIST: a method to inject new PTMs for DNA, RNA or protein into the dictionaries before
# launching the analysis. Could set a default lambda that aligns the whole residue.




