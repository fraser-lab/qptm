# run with python3

import sqlite3 as sql

def generate_qptxm_features_db(
  all_tested_filename="all_tested_ptms.out",
  synthetic_filename="synthetic_ptms.out"):
  """Read in features from all_tested_filename and look up true positives from
  synthetic_filename to generate the matching labels. Write out a SQLite
  database results.db to be concatenated with many others from other models.
  """
  con = None
  con = sql.connect('results.db')
  cur = con.cursor()

  # Create tables to hold the information contained in each file written
  # by qPTxM. We will construct a unique identifier for each tested
  # modification (incorporating model_id so it is unique across runs).
  cur.execute("CREATE TABLE all_tested(identifier TEXT PRIMARY KEY, \
    resname TEXT, short_name TEXT, cc REAL, d_ref REAL, d_mid REAL, \
    d_new_in_ref REAL, d_new_diff REAL, d_far REAL, ratio REAL, score REAL, \
    d_min REAL, b_factor REAL)")
  cur.execute("CREATE TABLE synthetic(identifier TEXT PRIMARY KEY, true INT)")

  # Read in all_tested_ptms.out for almost all necessary information
  with open(all_tested_filename, "r") as all_tested_file:
    for line in all_tested_file.readlines():
      parts = line.replace("'","''").split() # SQL standard is to escape ' as ''
      model_id = parts[14]
      chain_id, resid, resname, goto_atom, short_name = parts[0:5]
      identifier = "%s_%s_%s_%s_%s_%s" % \
        (model_id, chain_id, resid, resname, goto_atom, short_name)
      cc, d_ref, d_mid, d_new_in_ref, d_new_diff, d_far, ratio, score = \
        parts[6:14]
      d_min, b_factor = parts[15:17]
      query = "INSERT INTO all_tested VALUES('%s', '%s', '%s', %s, %s, %s, \
      %s, %s, %s, %s, %s, %s, %s)" % \
        (identifier, resname, short_name, cc, d_ref, d_mid, d_new_in_ref, \
          d_new_diff, d_far, ratio, score, d_min, b_factor)
      cur.execute(query)

  # Read in synthetic_ptms.out to look up true positives (all
  # modifications identified in this file are true positives).
  with open(synthetic_filename, "r") as synthetic_file:
    for i, line in enumerate(synthetic_file.readlines()):
      parts = line.replace("'","''").split()
      model_id = parts[6]
      chain_id, resid, resname, goto_atom, short_name = parts[0:5]
      identifier = "%s_%s_%s_%s_%s_%s" % \
        (model_id, chain_id, resid, resname, goto_atom, short_name)
      query = "INSERT INTO synthetic VALUES('%s', 1)" % identifier
      cur.execute(query)

  # Create a new table with just the features we care about,
  # keeping track of the true positives with the 'true' column.
  cur.execute("CREATE TABLE features AS SELECT all_tested.identifier, \
    resname, short_name, cc, d_ref, d_mid, d_new_in_ref, d_new_diff, d_far, \
    ratio, score, d_min, b_factor, true \
    FROM all_tested LEFT JOIN synthetic \
    ON all_tested.identifier = synthetic.identifier")

  # COALESCE returns the first non-NULL value. We'll use this
  # to mark the 'true' column in all other entries as 0.
  cur.execute("UPDATE features SET true = COALESCE(true, 0)")

  # # writing it back out to test:
  # cur.execute("SELECT * FROM features")
  # rows = cur.fetchall()
  # for row in rows[:10]:
  #   print(f"{row[1]} {row[2]} ... {row[13]}")

  # save results to disk and close the connection
  con.commit()
  con.close()

if __name__ == "__main__":
  import sys
  generate_qptxm_features_db(
    all_tested_filename=sys.argv[1],
    synthetic_filename=sys.argv[2])
