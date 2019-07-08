# qPTxM
qPTxM stands for quantification of posttranscriptional modifications. It is a program in python built upon the [computational crystallography toolbox (cctbx)](https://github.com/cctbx/cctbx_project) and [Phenix](https://www.phenix-online.org/).

**To use this software you will need a recent Phenix installation.**

To run qPTxM, please supply a model (pdb or mmcif) and map (mrc or ccp4) and specify an estimated median resolution:

```qptxm.py model.pdb map.mrc d_min=3```

Any modifications already present (if recognized) will be stripped off and searched for anew, and you will be shown some plots describing which possible modifications look promising and which were rejected. You can also load a copy of your model (ptms.pdb) with all proposed modifications modeled, and step through these positions with a custom script for Coot:

```coot goto_ptms.py ptms.pdb```

You may also like to load pruned.pdb to compare against the model with no modifications modeled. Under the menu option Custom, you will find the option to call up a panel of buttons to take you directly to each modeled modification.

Once you have decided on which proposed modifications to keep, remove any you want to reject from the file ptms.out and rerun qptxm.py with the additional argument selected_ptms=ptms.out (or whatever you may want to rename this). It will produce another model with just those modifications.

You may also find it useful to run with adjust_filters_only=True in order to test how many modifications are suggested if you adjust any of the optional parameters, like the minimum correlation coefficient of the model to the map (see below). This way is much faster but does not produce a new model and goto_ptms.py script. It also won't be able to adjust results based on an updated resolution estimate. You can also turn off plotting with plot=False or test choices of parameters on synthetic data (generated from your model but a calculated, noise-free map of the same dimensions) and see how accurately qptxm finds a randomly-modified 10%% of positions in that map by setting synthetic_data=True.

## Current release:
[![DOI](https://zenodo.org/badge/195718850.svg)](https://zenodo.org/badge/latestdoi/195718850)


## Citation:
*Submitted to bioRxiv and pending peer review.*
