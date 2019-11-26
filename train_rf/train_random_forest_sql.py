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

def reformat_as_dataframe(db_filename, labels=True, table="features"):
  """Concatenate features tables from many results.db files."""
  con = None
  con = sql.connect(db_filename)
  cur = con.cursor()

  # Gather combined features
  cur.execute("SELECT * FROM %s" % table)
  rows = cur.fetchall()

  # Save in a different format
  np_features = np.array([[hash(r[1]), hash(r[2])] + list(r[3:14]) for r in rows])
  df_categories = pd.DataFrame(
    np_features[:,0:2], columns=['resname_hash', 'short_name_hash'], dtype=int)
  df_floats = pd.DataFrame(
    np_features[:,2:12], columns=['cc', 'd_ref', 'd_mid', 'd_new_in_ref',
    'd_new_diff', 'd_far', 'ratio', 'score', 'd_min', 'b_factor'], dtype=float)
  if labels:
    df_labels = pd.DataFrame(
      np_features[:,12], columns=['labels'], dtype=int)
    df = pd.concat([df_categories, df_floats, df_labels], axis=1)
  else:
    df = pd.concat([df_categories, df_floats], axis=1)
  return df

def train_rf_on_dataframe(df):
  feature_names = ['resname_hash', 'short_name_hash', 'cc', 'd_ref', 'd_mid',
  'd_new_in_ref', 'd_new_diff', 'd_far', 'ratio', 'score', 'd_min', 'b_factor']
  X_train, X_test, y_train, y_test = train_test_split(
    df[feature_names],
    df['labels'],
    test_size=0.1,
    stratify=df['labels'],
    random_state=1000)
  rf = RandomForestClassifier(n_estimators=1000, oob_score=True)
  rf.fit(X_train, y_train)
  print("features and importances:")
  print(feature_names)
  print(rf.feature_importances_)
  picklefile_write = open("trained_rf.pkl", "wb")
  pickle.dump(rf, picklefile_write, protocol=4)
  print("wrote trained random forest classifier to trained_rf.pkl")
  picklefile_write.close()
  return rf

if __name__ == "__main__":
  import sys
  df_train = reformat_as_dataframe(sys.argv[1], table="combined_features")
  rf = train_rf_on_dataframe(df_train)
