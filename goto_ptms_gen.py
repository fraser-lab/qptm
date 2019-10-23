"""
Jiffy script to write a script to step through PTMs suggested by qptmd.py in Coot.
"""

template_start = """\"\"\"
Jiffy script to step through PTMs suggested by qptmd.py in Coot.
\"\"\"

def read_ptms(flatfile):
  import csv
  with open(flatfile, 'rb') as flat:
    reader = csv.reader(flat, delimiter=' ')
    chain = []
    resid = []
    atom = []
    abbr = []
    score = []
    for row in reader:
      if len(row) < 8: continue
      chain.append(row[0])
      resid.append(int(row[1]))
      atom.append(row[3])
      abbr.append(row[4])
      score.append(float(row[8]))
    name = ["%s %d %s %s %4.2f" % (ch, ri, at, ab, sc) \\
      for (ch, ri, at, ab, sc) in zip(chain, resid, atom, abbr, score)]
  sorted_score = sorted(score)[::-1]
  sorted_sites = []
  for sidx in xrange(len(sorted_score)):
    sc = sorted_score[sidx]
    uidx = score.index(sc)
    sorted_sites.append([
      chain[uidx],
      resid[uidx],
      atom[uidx],
      name[uidx],
      score[uidx],
      ])
    score[uidx] = None # prevent duplicate indices upon duplicate scores
  # sites = zip(chain, resid, atom, name, score)
  # sites_dict = {s[-1]:s for s in sites}
  return sorted_sites

# Make a button list for stepping through identified PTMS
def goto_ptms():"""
template_middle = """  sites = read_ptms("%s")
  mods_sites_buttons = []
  # super super hacky but all these lambdas have call time evaluation >:("""
template_end = """  generic_button_dialog("Possible modified sites:",mods_sites_buttons)

# Add the tool to the Display menu
menu=coot_menubar_menu("Custom")
submenu_display=gtk.Menu()
menuitem_2=gtk.MenuItem("Modifications...")
menuitem_2.set_submenu(submenu_display)
menu.append(menuitem_2)
menuitem_2.show()
add_simple_coot_menu_menuitem(submenu_display,
"Step through suggested posttranslational modifications", 
lambda func: goto_ptms())
"""

insert_unit = """  mods_sites_buttons.append([sites[%d][3], lambda func:
    set_go_to_atom_chain_residue_atom_name(sites[%d][0],sites[%d][1],sites[%d][2])])"""

def gen_from_ptms(ptms_flatfile="ptms.out"):
  nlines = 0
  with open(ptms_flatfile, 'r') as flat:
    for line in flat.readlines():
      nlines += 1
  sections = [template_start]
  sections.append(template_middle % ptms_flatfile)
  for i in range(nlines):
    sections.append(insert_unit % (i, i, i, i))
  sections.append(template_end)
  combined = "\n".join(sections)
  with open("goto_ptms.py", 'w') as pyfile:
    pyfile.write(combined)
  print ("\nWrote file goto_ptms.py. Pass this file to Coot to run.")

if __name__ == "__main__":
  import sys
  gen_from_ptms(ptms_flatfile=sys.argv[1])
