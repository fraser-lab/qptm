# run with python3

import sqlite3 as sql

def combine_features_tables(*paths_to_dbs):
  """Concatenate features tables from many results.db files."""
  con = None
  con = sql.connect('combined_results.db')
  cur = con.cursor()

  # Create a table matching the schema of the 'features' tables
  cur.execute("CREATE TABLE combined_features(id INT PRIMARY KEY, \
    resname TEXT, short_name TEXT, cc REAL, d_ref REAL, d_mid REAL, \
    d_new_in_ref REAL, d_new_diff REAL, d_far REAL, \
    ratio REAL, score REAL, d_min REAL, b_factor REAL, true INT)")

  # Attach the other databases one by one and add them in
  for path in paths_to_dbs:
    cur.execute("ATTACH DATABASE '%s' AS this_results" % path)
    cur.execute("INSERT INTO combined_features(\
      resname, short_name, cc, d_ref, d_mid, \
      d_new_in_ref, d_new_diff, d_far, ratio, score, \
      d_min, b_factor, true) \
      SELECT resname, short_name, cc, d_ref, d_mid, \
      d_new_in_ref, d_new_diff, d_far, ratio, score, \
      d_min, b_factor, true \
      FROM this_results.features")
    con.commit()
    cur.execute("DETACH this_results")

  # save results to disk and close the connection
  con.commit()
  con.close()

if __name__ == "__main__":
  import sys
  combine_features_tables(*sys.argv[1:])
