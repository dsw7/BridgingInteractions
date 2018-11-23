"""
Written by David S. Weber
    
(1) Script downloads a JSON file of all PDB codes from PDB
(2) Script converts JSON to list
(3) Script sets up a SQL database for storing data
(4) Script iterates for entries in list
(5) Script downloads the PDB file corresponding to each iteration
(6) Script applies Met-Aromatic algorithm to each local PDB file
(7) Script finds closely spaced Tyr/Trp residues in each local PDB file
(8) Script filters in only 2-bridges from Met-Aromatic output
(9) Chain / bridge data is exported for analysis


Testing:
Previous results.txt:
    
1BS2 : NR : {'TYR362', 'TYR188'} : {'TRP192', 'TRP266'} : 6.1.1.19 : SACCHAROMYCES CEREVISIAE;  
1BS2 : NR : {'TYR362', 'TYR188'} : {'TYR491', 'TYR440'} : 6.1.1.19 : SACCHAROMYCES CEREVISIAE;  
1BS2 : NR : {'TYR362', 'TYR188'} : {'TYR534', 'TYR553', 'TYR134'} : 6.1.1.19 : SACCHAROMYCES CEREVISIAE;  
1BS2 : NR : {'TYR362', 'TYR188'} : {'TRP393', 'TYR359'} : 6.1.1.19 : SACCHAROMYCES CEREVISIAE;  
1BS2 : NR : {'TYR362', 'TYR188'} : {'TYR281', 'TYR288', 'TYR277', 'TYR291'} : 6.1.1.19 : SACCHAROMYCES CEREVISIAE;  
1BS2 : NR : {'TYR362', 'TYR188'} : {'TRP181', 'TYR176'} : 6.1.1.19 : SACCHAROMYCES CEREVISIAE;  
1BS2 : IS : {'TYR362', 'TYR188'} : {'TYR347', 'TYR188'} : 6.1.1.19 : SACCHAROMYCES CEREVISIAE;  
1BS2 : NR : {'TYR362', 'TYR188'} : {'TYR585', 'TYR565'} : 6.1.1.19 : SACCHAROMYCES CEREVISIAE;  

Current results.txt:
    
1BS2 : Chains : TYR288, TYR277, TYR291, TYR281,      y
1BS2 : Chains : TYR553, TYR534, TYR134,              y
1BS2 : Chains : TRP393, TYR359,                      y
1BS2 : Chains : TYR491, TYR440,                      y
1BS2 : Chains : TYR347, TYR188,                      y
1BS2 : Chains : TRP181, TYR176,                      y
1BS2 : Chains : TYR585, TYR565,                      y
1BS2 : Chains : TRP266, TRP192,                      y

1BS2 : Bridges : TYR362, TYR188,     <- match
1BS2 : Bridges : TRP192, PHE252,     <- new
1BS2 : Bridges : TYR369, PHE364,     <- new

"""

# ------------------------------------------------------------------------------
# input parameters to be specified by researcher


# search parameters
# --------------
START = 0
END   = -1  # set to -1 to iterate over entire PDB


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
from operator               import itemgetter
from itertools              import groupby, chain
import networkx as nx


# ------------------------------------------------------------------------------

# create a text file for dumping results
# --------------------------------------
f = open('results_superimposition.txt', 'w')
f.close()


# get current list of all PDB codes from RCSB PDB
# -----------------------------------------------
url = 'https://www.rcsb.org/pdb/json/getCurrent'
current_PDB_files_JSON       = get(url)
current_PDB_files_pyDict     = loads(current_PDB_files_JSON.content)
current_PDB_files_pyList     = current_PDB_files_pyDict.get('idList')
current_PDB_files_entrycount = current_PDB_files_pyDict.get('resultCount')
if END == -1: END = current_PDB_files_entrycount


# iterate over all PDB files
# --------------------------
for iteration, CODE in enumerate(current_PDB_files_pyList[START:END]):
    status = '{} - Analyzed {} of {} PDB entries.'.format(CODE, iteration + 1, current_PDB_files_entrycount)
    print(status)
    
    file_pdb = PDBFile(CODE)
    path_to_file = file_pdb.fetch_from_PDB()
    
    if path_to_file == 'URLError':  # next iteration if file no longer exists
        print('URLError')
        continue
    
    try:  # catch other exceptions such as "string crashes", missing coordinates, etc., in PDB files
        chains = get_nn(filepath=path_to_file, cutoff=7.4)
        if chains == []: file_pdb.clear(); continue        # next iteration if no chains
        
        interactions = met_aromatic(path_to_file, *args)
        if interactions == []: file_pdb.clear(); continue  # next iteration if no met-aromatic interactions
    
    except Exception as exception: 
        print(exception)
        file_pdb.clear()
        continue
        
    # isolate bridging pairs
    pairs = [('{}{}'.format(i[0], i[1]), '{}{}'.format(i[2], i[3])) for i in interactions]
    pairs = list(set(pairs))  # remove identical lines due to multiple vectors v
    G1 = nx.Graph()
    G1.add_edges_from(pairs)
    bridges = []
    for disconnects in list(nx.connected_components(G1)):
        if len(disconnects) == 3: bridges.append(list(disconnects))                
    
    # next iteration if no 2-bridges
    if bridges == []: file_pdb.clear(); continue
    
    # remove MET data
    for nx_sets in bridges:
        for entry in nx_sets.copy():
            if 'MET' in entry: nx_sets.remove(entry)
                    
    # remove inverse bridges
    bridges = list(filter(lambda entry: len(entry) > 1, bridges))
    
    # next iteration if no bridges remain
    if bridges == []: file_pdb.clear(); continue

    # get disconnected components
    G2 = nx.Graph()
    G2.add_edges_from(chains)
    chains = list(nx.connected_components(G2))
    
    # write to file
    f = open('results_superimposition.txt', 'a')
    for c in chains: f.write('{} : Chains : {} \n'.format(CODE, ('{}, ' * len(c)).format(*c)))
    for b in bridges: f.write('{} : Bridges : {} \n'.format(CODE, ('{}, ' * len(b)).format(*b)))
    f.close()
    
    file_pdb.clear()
    
