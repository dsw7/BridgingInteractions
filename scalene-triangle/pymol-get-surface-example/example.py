# dsw written example for fetching protein surface

from pymol   import cmd, stored
from inspect import stack
from os      import path
filepath_inspected = stack()[0][1]
root = path.dirname(filepath_inspected)


# input
# -----
pdbfile = '1rcy.pdb'
chain = 'A'
area_SA = 39.0 # I set solvent exposure to a very high 39 A for illustrative purposes


# load file into pymol
# --------------------
file = path.join(root, pdbfile)
cmd.load(file)


# create a pymol temporary object
# -------------------------------
tmpObj = "__tmp"


# isolate protein parts
# ---------------------
# "(%s and polymer) and not resn HOH" # include all chains
cmd.create(tmpObj, "({} and polymer and chain {}) and not resn HOH".format(pdbfile.split('.')[0], chain))


# get the surface using pymol internals
# -------------------------------------
cmd.set("dot_solvent")
cmd.get_area(selection=tmpObj, load_b=1)


# atom exposure threshold in Angstroms squared
# --------------------------------------------
cmd.remove(tmpObj + " and b < " + str(area_SA))


# show the [very] solvent exposed residues
# ----------------------------------------
cmd.show(selection=tmpObj, representation="dots") # show the exposed atoms


# store the exposed data
# ----------------------
ITER_STATE_EXP = "stored.tmp_dict[(chain, resv, name, x, y, z)] = 1"
stored.tmp_dict = {}
cmd.iterate_state(state=-1, selection=tmpObj, expression=ITER_STATE_EXP)
exposed = stored.tmp_dict.keys()
exposed.sort()


# print example data to console
# -----------------------------
for item in exposed:
    print item


# clear buffer and return data
# cmd.delete(tmpObj)   # destroys tmpObj - no idea if this is needed
# cmd.delete('all')  # removes everything in PyMOL viewer



