import pickle
import numpy as np

class ModificationClassifier(object):
  def __init__(self, rf):
    self.rf = rf
  def predict(self, features_as_lists):
    assert len(features_as_lists) == 17, "wrong number of features passed for prediction"
    features_hashes = [list(map(hash, features_as_lists[2]))] + \
                      [list(map(hash, features_as_lists[4]))] + \
                      [list(map(float, l)) for l in features_as_lists[6:15]] + \
                      [list(map(float, l)) for l in features_as_lists[16:18]]
    features_np = np.asarray(list(zip(*features_hashes)))
    predicted = self.rf.predict(features_np)
    probabilities = self.rf.predict_proba(features_np)
    features_with_predictions = \
      list(features_as_lists) + [predicted, probabilities[:,1]]
    return features_with_predictions
  def write_predicted(self, features_with_predictions, threshold=None,
    predicted_file="predicted_ptms.out"):
    if threshold is None:
      keep = lambda feature: feature[-2] == 1
    else:
      keep = lambda feature: feature[-1] >= threshold
    with open(predicted_file, "w") as out:
      for feature in zip(*features_with_predictions):
        if keep(feature):
          out.write(" ".join(list(map(str, feature[:-2]))) + "\n")
  def write_features_with_predictions(self, features_with_predictions,
    predicted_file="all_tested_ptms_with_predictions.out"):
    with open(predicted_file, "w") as out:
      for feature in zip(*features_with_predictions):
        out.write(" ".join(list(map(str, feature))) + "\n")

def write_predictions_for_all_tested_ptms(rf, ptms_file):
  features_as_lists = import_ptms(ptms_file)
  this_classifier = ModificationClassifier(rf)
  features_with_predictions = this_classifier.predict(features_as_lists)
  this_classifier.write_predicted(features_with_predictions)
  this_classifier.write_features_with_predictions(features_with_predictions)

if __name__ == "__main__":
  # this needs python3
  import os, sys
  from filter_util_no_flex import import_ptms
  from goto_ptms_gen import gen_from_ptms

  rf_basedir, ptms_file = sys.argv[1:3]

  # this next bit in particular **will not work in python 2.7**
  with open(os.path.join(
    rf_basedir, "synthetic_expanded_severaldmin_B10_combined_rf_py3.pickle"),
    "rb") as picklefile_read:
    this_rf = pickle.load(picklefile_read)

  write_predictions_for_all_tested_ptms(this_rf, ptms_file)
  gen_from_ptms(ptms_flatfile="predicted_ptms.out")
  print("Please rerun qptxm.py with the additional argument " +\
    "selected_ptms=predicted_ptms.out to generate an updated ptms.pdb file.")
