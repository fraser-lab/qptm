from __future__ import division
from matplotlib import pyplot as plt
from scitbx.array_family import flex

def plot_densities_from_flatfile(flatfile):
  import csv
  densities1 = flex.double()
  densities2 = flex.double()
  scores = flex.double()
  with open(flatfile, 'rb') as flat:
    reader = csv.reader(flat, delimiter=' ')
    for row in reader:
      if len(row) < 9: continue
      densities1.append(float(row[6]))
      densities2.append(float(row[7]))
      scores.append(float(row[8]))
  mod_densities1 = max(densities1, flex.double(len(densities1), 0.0000001))
  ratios = densities2/mod_densities1
  for (array, xaxis) in ((densities1, "densities1"), (densities2, "densities2"),
    (ratios, "ratios")):
    n, bins, patches = plt.hist(array, 20, facecolor='b', alpha=0.5)
    plt.xlabel(xaxis)
    plt.ylabel("frequency")
    plt.show()
  with open("ccs.out", 'rb') as ccs:
    reader = csv.reader(ccs, delimiter=' ')
    values = flex.double()
    for row in reader:
      values.append(float(row[0]))
    n, bins, patches = plt.hist(values, 20, facecolor='b', alpha=0.5)
    plt.xlabel("correlation coefficients")
    plt.ylabel("frequency")
    plt.show()
