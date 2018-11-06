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
    for row in reader:
      if len(row) < 8: continue
      chain.append(row[0])
      resid.append(int(row[1]))
      atom.append(row[3])
      abbr.append(row[4])
    name = ["%s %d %s %s" % (ch, ri, at, ab) \
      for (ch, ri, at, ab) in zip(chain, resid, atom, abbr)]
  sites = zip(chain, resid, atom, name)
  # sites_dict = {s[-1]:s for s in sites}
  return sites

#Make a button list for stepping through identified PTMS
def goto_ptms():
  sites = read_ptms(os.path.join(os.getcwd(), "ptms.out"))
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

if __name__ == "__main__":
  nlines = 0
  with open("ptms.out", 'rb') as flat:
    for line in flat.readlines():
      nlines += 1
  sections = [template_start]
  for i in xrange(nlines):
    sections.append(insert_unit % (i, i, i, i))
  sections.append(template_end)
  combined = "\n".join(sections)
  with open("goto_ptms.py", 'wb') as pyfile:
    pyfile.write(combined)
  print "wrote file goto_ptms.py. Pass this file to Coot to run."
