from __future__ import division
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from scitbx.array_family import flex
from filter_util import import_ptms

def plot_densities_from_flatfile(flatfile_tested, flatfile_accepted, cc_threshold):
  tested_ptms = import_ptms(flatfile_tested)
  accepted_ptms = import_ptms(flatfile_accepted)
  cmap = plt.get_cmap('plasma')
  test_color = cmap(0.05)
  accept_color = cmap(0.9)
  colors = (test_color, accept_color)
  alphas = (0.8, 0.5)
  legend_tags = ["tested (%d)" % len(tested_ptms[0]), "accepted (%d)" % len(accepted_ptms[0])]
  # chain_id, resid, resname, goto_atom, short_name, full_name, cc, d_ref, d_mid, \
  #   d_new_in_ref, d_new_diff, d_far, ratio, score = ptms
  for (xaxis, figname, array_index) in zip(
    ("Correlation coefficients of all residues to the map",
      "Reference position densities (rmsd)",
      "Proposed modification position densities (rmsd, experimental map)",
      "Proposed modification position densities (rmsd, difference map)",
      "Scores for modifications based on densities"),
    ("CCs.pdf", "Ref_densities.pdf", "New_densities.pdf", "Scores.pdf"),
    (6, 7, 9, 10, 13)):
    for (legend, array, color, alpha) in zip(
      legend_tags,
      (tested_ptms[array_index], accepted_ptms[array_index]),
      colors, alphas):
      not_too_small_part = array.select(array > 0.01)
      n, bins, patches = plt.hist(not_too_small_part, max(10, len(not_too_small_part)//30),
        facecolor=color, alpha=alpha)
      if array_index == 6 and legend == "tested": # CC
        annotation_height = 0.75*max([p._height for p in patches])
        plt.text(cc_threshold-0.02, annotation_height,
          "Threshold for testing a given\nresidue for modifications",
          horizontalalignment='right')
        plt.axvline(x=cc_threshold, color='k')
    handles = (Rectangle((0,0),1,1,color=c) for c in colors)
    labels = legend_tags
    plt.legend(handles, labels)
    plt.xlabel(xaxis)
    plt.ylabel("Frequency")
    plt.savefig(figname, dpi=300)
    plt.show()
  return

if __name__ == "__main__":
  import sys
  plot_densities_from_flatfile(sys.argv[1], sys.argv[2], float(sys.argv[3]))
