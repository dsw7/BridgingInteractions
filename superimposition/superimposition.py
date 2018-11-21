"""
Written by David S. Weber
    
(1) Script downloads a JSON file of all PDB codes from PDB
(2) Script converts JSON to list
(3) Script sets up a SQL database for storing data
(4) Script iterates for entries in list
(5) Script downloads the PDB file corresponding to each iteration
(6) Script applies Met-Aromatic algorithm to each local PDB file
(7) Script finds closely spaced Tyr/Trp residues in each local PDB file
(8) Script compares bridge/chain membership
"""

# TODO: I put a benchmark flag in nn.py -> check this

# ------------------------------------------------------------------------------
# input parameters to be specified by researcher


# search parameters
# --------------
START = 0
END   = 30  # set to -1 to iterate over entire PDB


# PDB parameters
# --------------
CHAIN = 'A'       # input a chain of interest here
CUTOFF = 6.0      # input some cutoff vector v norm
ANGLE = 109.5     # input some cutoff vector a / vector v and vector g / vector v angle
MODEL = 'cp'      # 'cp' or 'rm' -> Cross Product or Rodrigues' Method lone pair interpolation methods
args = (CHAIN, CUTOFF, ANGLE, MODEL)


# ------------------------------------------------------------------------------

from sys                    import path as PATH; PATH.append('libs')
from requests               import get
from json                   import loads
from PDB_filegetter         import PDBFile
from ma_lowlevel            import met_aromatic
from platform               import node
from time                   import sleep
from nn                     import get_nn


# ------------------------------------------------------------------------------

# TODO: prepare SQL database here

# ------------------------------------------------------------------------------

# get current list of all PDB codes from RCSB PDB
# -----------------------------------------------
url = 'https://www.rcsb.org/pdb/json/getCurrent'
current_PDB_files_JSON       = get(url)
current_PDB_files_pyDict     = loads(current_PDB_files_JSON.content)
current_PDB_files_pyList     = current_PDB_files_pyDict.get('idList')
current_PDB_files_entrycount = current_PDB_files_pyDict.get('resultCount')
if END == -1: END = current_PDB_files_entrycount


for iteration, CODE in enumerate(current_PDB_files_pyList[START:END]):
    status = '{} - Analyzed {} of {} PDB entries.'.format(CODE, iteration + 1, current_PDB_files_entrycount)
    print(status)
    
    file_pdb = PDBFile(CODE)
    path_to_file = file_pdb.fetch_from_PDB()
    
    if path_to_file == 'URLError':  # next iteration if file no longer exists
        print('URLError')
        continue
    
    else:
        try:  # catch other exceptions such as "string crashes", missing coordinates, etc., in PDB files
            chains = get_nn(filepath=path_to_file, cutoff=7.4)
            if chains == []:        # next iteration if no chains in file
                file_pdb.clear()
                continue
            
            interactions = met_aromatic(path_to_file, *args)
            if interactions == []:  # next iteration if no interactions in file
                file_pdb.clear()
                continue
        
                    
        except Exception as exception: 
            print(exception)
        file_pdb.clear()
    
