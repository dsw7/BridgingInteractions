"""
A secondary script to dsw_pymol_bridging.py
Here we need to manually parse bridging_n_n+5000.csv files because
PyMOL's Python 2 installation is not compatible with our Python 3 Pandas
installation.

This script requires an installation of Numpy for Python2.7
    $ python2.7 -m pip install numpy
Which created a Numpy installation in PyMOL's Python2.7 import path which
can be found using the following:
    $ python2.7
    >>> import numpy
    >>> print numpy.__path__
    
I did the same for Pandas:
    $ python2.7 -m pip install pandas
    
    
    
    
Rewrite this script one more time for GitHub SI
===============================================
* The problem with this script is that it's not super readable *

Import PDB file -> PDB FILE OBJECT

PDB FILE OBJECT -> get bridges
PDB FILE OBJECT -> get metals
PDB FILE OBJECT -> use get_solvent_exposed to get surface coordinate

If all three conditions are true then collect stats for the scalene triangle dimension

Dump PDB FILE OBJECT
"""

from pymol import cmd, stored
import os
import numpy as np  # will import version for Python 2 or 3 depending on loc
import pandas as pd
import time
import sys

# path to where individual .csvs should be saved
root_output = r'/Volumes/MSC/Databases/2018 Bridging Databases/Geometric Data'

ITER_STATE_EXP = "stored.tmp_dict[(chain, resv, elem, x, y, z)] = 1"

OUT_COL_LABELS = ['PDB CODE', 
                  'SD_x', 'SD_y', 'SD_z', 
                  'Z_x', 'Z_y', 'Z_z',
                  'S_x', 'S_y', 'S_z', 
                  'MET RES', 'SUR RES', 
                  'ARO A', 'ARO A RES', 'ARO B', 'ARO B RES']

# -------------------------------------------------------

# move to dependency dir?
def error_handler(exception):
    # returns the exception type and line number
    _, _, exc_traceback = sys.exc_info()   
    exc_a = ' {}\n'.format(exception)
    exc_b = ' Is occurring on line {}'.format(exc_traceback.tb_lineno)
    print ' The error:\n' + exc_a + exc_b
    del(exc_traceback)  # delete from Python garbage collector

# move to dependency dir?
def get_solvent_exposed(code, cutoff=2.5):
    # the file, 1rcy.pdb, needs to be within the current working directory
    # prior to calling cmd.load()
    cmd.fetch(code)  # get from ftp client
    time.sleep(0.05)
    
    # search for different file formats in pwd
    for ext in ('.pdb', '.cif'):
        if os.path.exists(os.path.join(os.getcwd(), code + ext)):
            new_path = code + ext
            cmd.load(new_path)  # load into pymol - might not be needed?

    # create a pymol temporary object
    tmpObj = "__tmp"

    # "(%s and polymer) and not resn HOH" # include all chains
    cmd.create(tmpObj, "(%s and polymer and chain A) and not resn HOH" % code)

    # get the surface using pymol internals
    cmd.set("dot_solvent")
    cmd.get_area(selection=tmpObj, load_b=1)

    # atom exposure threshold in Angstroms squared
    cmd.remove(tmpObj + " and b < " + str(cutoff))

    # store the exposed data
    stored.tmp_dict = {}
    # cmd.iterate(tmpObj, "stored.tmp_dict[(chain, resv, elem)] = 1") # maybe here?
    cmd.iterate_state(state=-1, selection=tmpObj, expression=ITER_STATE_EXP)
    exposed = stored.tmp_dict.keys()
    exposed.sort()

    # remove the file when done
    path_to_file = os.path.join(os.getcwd(), new_path)
    if os.path.exists(path_to_file):
        os.remove(path_to_file)

    # clear buffer and return data
    cmd.delete(tmpObj)   # destroys tmpObj - no idea if this is needed
    cmd.delete('all')  # removes everything in PyMOL viewer
    return exposed

# -------------------------------------------------------

root_path = r'/Volumes/MSC/Databases/2018 Bridging Databases'
# root_path = r'G:/Databases/2018 Bridging Databases'

# paths to bridging_n_n+5000.csv datasets
paths = os.listdir(root_path)
paths = [item for item in paths if item[-4:] == '.csv']

# paths to .csvs in bridging_proc
paths_to_bridging = os.path.join(root_path, 'bridging_proc')
paths_to_bridging = os.listdir(paths_to_bridging)

# -------------------------------------------------------

# get all metals
df_METALS = pd.DataFrame()
for p in paths:
    if p[0:2] == '._':  # pass over DS store settings storage files
        pass
    else:
        df_m = pd.read_csv(os.path.join(root_path, p))  
        df_m = df_m[df_m['Z-Metal'].notnull()]
        df_METALS = df_METALS.append(df_m)

df_METALS = df_METALS.reset_index(drop=True)

# get all bridges
df_BRIDGES = pd.DataFrame()
for p in paths_to_bridging:
    if p[0:2] == '._':  # pass over DS store settings storage files
        pass
    else:
        df_b = pd.read_csv(os.path.join(root_path, 'bridging_proc', p))  # angular masked? <- yes
        df_BRIDGES = df_BRIDGES.append(df_b)

# core loop
for idx in range(0, df_METALS.shape[0]):
    # create a temporary df in a fashion similar to Pandas groupby machinery
    df = df_BRIDGES[df_BRIDGES['0'] == df_METALS.iloc[idx, 1]]
    PDB_CODE = df_METALS.iloc[idx, 1][0:4]
    print PDB_CODE
    if df.empty == False:

        """
        This is not the best approach.
        Protein A chains that contain N number of metals will result in
        N number of PDB file downloads to pwd in addition to N number of
        os.remove() calls. In this case I prefer ease of programming and
        reading at the expense of efficiency however. Moreover, there's
        only about 9,000 files to go through.
        """

        try:
            # get entire protein surface
            exposed_data = get_solvent_exposed(PDB_CODE, cutoff=2.50)

            MET_DATA = df.drop_duplicates(('MET RES', 'SD_x', 'SD_y', 'SD_z'))
            MET_DATA = MET_DATA[['MET RES', 'SD_x', 'SD_y', 'SD_z']].values

            for met_data in MET_DATA:  # loop to iterate over multiple SDs
                MET_RES = met_data[0]
                vec_a = met_data[1:4]  # met SD coordinates
                vec_b = df_METALS.iloc[idx, 21:24].values  # metal coordinates

                # export aromatics that are bridging MET_RES
                df_aromatics = df[df['MET RES'] == MET_RES][['ARO', 'ARO RES']].drop_duplicates().values
                ARO_A = df_aromatics[0]
                ARO_B = df_aromatics[1]

                # locate the minimum from protein surface to the MET SD in current iteration
                norms = [(idx_2, np.linalg.norm(coord[3:6] - vec_a)) for idx_2, coord in enumerate(exposed_data)]
                min_norm = min(norms, key = lambda t: t[1])
                vec_c = exposed_data[min_norm[0]][3:6]  # coord for minimum to SD
                SUR_RES = exposed_data[min_norm[0]][1]  # get the residue position for the surface residue

                # export all output data
                output = [PDB_CODE] + list(vec_a) + list(vec_b) + list(vec_c) + \
                         [MET_RES, SUR_RES, ARO_A[0], ARO_A[1], ARO_B[0], ARO_B[1]]                                                            
                df_out = pd.DataFrame(output).transpose()
                df_out.columns = OUT_COL_LABELS

                # a terminator for files with multiple metals
                # prevents files from being overwritten
                metal_termin = str(df_METALS.iloc[idx, 21]).replace('.', '')[0:4]

                # add aromatic labels to filenames
                c1 = str(ARO_A[0]) + str(ARO_A[1])
                c2 = str(ARO_B[0]) + str(ARO_B[1])

                out_filename = 'brdg_geom_{}_{}_{}_{}.csv'.format(PDB_CODE, c1, c2, metal_termin)
                df_out.to_csv(os.path.join(root_output, out_filename))

        except Exception as exception:
            error_handler(exception)
            continue

