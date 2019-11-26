# run with python3

# General data manipulation packages
import numpy as np
import pandas as pd
import sqlite3 as sql

# Useful function to split data sets
from sklearn.model_selection import train_test_split

# Random Forest class:
from sklearn.ensemble import RandomForestClassifier

# For pickling results
import pickle

# Function to calculate the prediction accuracy from the classifier
from sklearn.metrics import accuracy_score

# Function to calculate a confusion matrix
from sklearn.metrics import confusion_matrix

# Packages used for plotting
import matplotlib.pyplot as plt
import seaborn as sns

# Function to read in the features in a pandas dataframe
from train_random_forest_sql import reformat_as_dataframe

class RFPredictions(object):
  """Keep track of a random forest's predictions on a test data set
  for analysis and visualization"""
  def __init__(self, rf, data_test):
    self.rf = rf
    self.data = data_test
    try:
      self.labels = data_test['labels']
    except Exception:
      self.labels = None
    other_columns = [c for c in data_test.columns if c != 'labels']
    self.data = data_test[other_columns]
  def predict(self):
    self.predictions = self.rf.predict(self.data)
    self.probabilities = self.rf.predict_proba(self.data)
    if self.labels is not None:
      self.accuracy = accuracy_score(self.labels, self.predictions)
  def print_score(self):
    print('Out-of-bag score estimate: %.3f' % self.rf.oob_score_)
    if self.labels is not None:
      print('Mean accuracy score: %.3f' % self.accuracy)
  def plot_confusion_matrix(self, title='confusion matrix', figname='cm.pdf', show=False):
    if self.labels is None:
      print('Cannot plot confusion matrix without data labels (true positives).')
      return
    cm_labels = ['unmodified', 'modified']
    cm = pd.DataFrame(confusion_matrix(self.labels, self.predictions))
    print(cm[::-1])

    plt.figure(figsize = (10,8))
    sns.set(font_scale=1.4)#for label size
    ax = plt.subplot()
    sns.heatmap(cm, annot=True, annot_kws={"size": 12}, ax=ax, fmt='g')# font size, axes, number format

    # labels, title and ticks
    ax.set_xlabel('Predicted')
    ax.set_ylabel('True')
    ax.set_title(title);
    ax.xaxis.set_ticklabels(['unmodified', 'modified'])
    ax.yaxis.set_ticklabels(['unmodified', 'modified'])
    ax.set(ylim=(0,2))
    plt.savefig(figname, dpi=300)
    print('Confusion matrix figure saved to %s' % figname)
    if show:
      plt.show()

def filter_all_tested_to_match_preds(predictions, all_tested_filename, new_filename):
  with open(all_tested_filename, "r") as old_file:
    all_tested = old_file.readlines()
  with open(new_filename, "w") as new_file:
    filter_print = lambda line: new_file.write(" ".join(line.split()[:17]) + "\n")
    for i in range(len(predictions)):
      if predictions[i]:
        filter_print(all_tested[i])
  print("wrote %s" % new_filename)


def plot_predictions(rf_filename, db_filename, title=None, show=False):
  with open(rf_filename, "rb") as rf_file:
    rf = pickle.load(rf_file)
  features_df = reformat_as_dataframe(db_filename, labels=True)
  predictions = RFPredictions(rf, features_df)
  predictions.predict()
  filter_all_tested_to_match_preds(predictions.predictions, "all_tested_ptms.out", "preds_by_rf.out")
  predictions.print_score()
  predictions.plot_confusion_matrix(title=title, show=show)

if __name__ == "__main__":
  import sys
  rf_filename, db_filename, plot_title = sys.argv[1:4]
  plot_predictions(rf_filename, db_filename, title=plot_title)
