from __future__ import division
from matplotlib import pyplot as plt
from scitbx.array_family import flex

def plot_densities_from_flatfile(flatfile, cc_threshold):
  import csv
  densities1 = flex.double()
  densities2 = flex.double()
  scores = flex.double()
  with open(flatfile, 'rb') as flat:
    reader = csv.reader(flat, delimiter=' ')
    for row in reader:
      if len(row) < 9 : continue
      densities1.append(float(row[6]))
      densities2.append(float(row[7]))
      scores.append(float(row[8]))
  ratios = densities2/densities1
  nbins = max(10, len(densities1)//30)
  for (array, xaxis, figname) in (
      (densities1, "Reference position densities (rmsd)", "Ref_densities.pdf"),
      (densities2, "Proposed modification position densities (rmsd)", "New_densities.pdf"),
      (ratios, "Ratios of densities at proposed and reference positions", "Ratios.pdf")):
    n, bins, patches = plt.hist(array, nbins, facecolor='b', alpha=0.5)
    plt.xlabel(xaxis)
    plt.ylabel("Frequency")
    plt.savefig(figname, dpi=300)
    plt.show()
  with open("ccs.out", 'rb') as ccs:
    reader = csv.reader(ccs, delimiter=' ')
    values = flex.double()
    for row in reader:
      values.append(float(row[0]))
    n, bins, patches = plt.hist(values, max(10, len(values)//30), facecolor='b', alpha=0.5)
    plt.xlabel("Correlation coefficients of all residues to the map")
    annotation_height = 0.85*max([p._height for p in patches])
    plt.text(cc_threshold-0.02, annotation_height,
      "Threshold for testing a given\nresidue for modifications",
      horizontalalignment='right')
    plt.ylabel("Frequency")
    plt.axvline(x=cc_threshold, color='k')
    plt.savefig("CCs.pdf", dpi=300)
    plt.show()

if __name__ == "__main__":
  import sys
  plot_densities_from_flatfile(sys.argv[1], float(sys.argv[2]))
