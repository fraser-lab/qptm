"""
Jiffy script to step through PTMs suggested by qptmd.py in Coot.
"""

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
  # super super hacky but all these lambdas have call time evaluation >:(
  try:
    mods_sites_buttons.append([sites[0][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[0][0],sites[0][1],sites[0][2])])
    mods_sites_buttons.append([sites[1][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[1][0],sites[1][1],sites[1][2])])
    mods_sites_buttons.append([sites[2][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[2][0],sites[2][1],sites[2][2])])
    mods_sites_buttons.append([sites[3][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[3][0],sites[3][1],sites[3][2])])
    mods_sites_buttons.append([sites[4][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[4][0],sites[4][1],sites[4][2])])
    mods_sites_buttons.append([sites[5][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[5][0],sites[5][1],sites[5][2])])
    mods_sites_buttons.append([sites[6][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[6][0],sites[6][1],sites[6][2])])
    mods_sites_buttons.append([sites[7][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[7][0],sites[7][1],sites[7][2])])
    mods_sites_buttons.append([sites[8][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[8][0],sites[8][1],sites[8][2])])
    mods_sites_buttons.append([sites[9][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[9][0],sites[9][1],sites[9][2])])
    mods_sites_buttons.append([sites[10][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[10][0],sites[10][1],sites[10][2])])
    mods_sites_buttons.append([sites[11][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[11][0],sites[11][1],sites[11][2])])
    mods_sites_buttons.append([sites[12][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[12][0],sites[12][1],sites[12][2])])
    mods_sites_buttons.append([sites[13][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[13][0],sites[13][1],sites[13][2])])
    mods_sites_buttons.append([sites[14][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[14][0],sites[14][1],sites[14][2])])
    mods_sites_buttons.append([sites[15][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[15][0],sites[15][1],sites[15][2])])
    mods_sites_buttons.append([sites[16][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[16][0],sites[16][1],sites[16][2])])
    mods_sites_buttons.append([sites[17][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[17][0],sites[17][1],sites[17][2])])
    mods_sites_buttons.append([sites[18][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[18][0],sites[18][1],sites[18][2])])
    mods_sites_buttons.append([sites[19][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[19][0],sites[19][1],sites[19][2])])
    mods_sites_buttons.append([sites[20][3], lambda func: \
      set_go_to_atom_chain_residue_atom_name(sites[20][0],sites[20][1],sites[20][2])])
  except IndexError:
    pass
  generic_button_dialog("Possible modified sites:",mods_sites_buttons)

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
