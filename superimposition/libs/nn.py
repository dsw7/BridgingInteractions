# Written by David Weber
# dsw7@sfu.ca

"""
This script iteratively finds closely spaced Tyr/Trp chains for some
cutoff distance CUTOFF_DIST. Here it is assumed that the PDB file
is already in the script's pwd as a result of some other operation.
"""

from operator     import itemgetter
from itertools    import groupby
from re           import search
from numpy.linalg import norm
from numpy        import array

IDX_ATOM = 0
IDX_CHAIN = 4
IDX_AA = 3
IDX_ATM_LABEL = 2
ATOMS_TYR = r'CG|CZ'  
ATOMS_TRP = r'CD2|CH2'

def groupby_to_midpoint(nested_list):
    # TODO: benchmark this
    # groupby object -> tuple of res/pos, midpoint
    v1 = array(nested_list[0][6:9]).astype(float)
    v2 = array(nested_list[1][6:9]).astype(float)
    return (nested_list[0][3] + nested_list[0][5], 0.5 * (v1 + v2))

def get_nn(filepath, cutoff, chain='A'):
    """
    Function gets a set of Tyr/Trp nearest neighbors from a PDB entry
    Parameters:
        filepath -> path to the PDB file, .ent format
        cutoff   -> the maximum distance deemed neighboring
        chain    -> 'A', 'B', 'C', etc.
    Returns:
        A list of neighboring Tyr/Trp residues
    """
    with open(filepath, 'r') as f:
        data_incoming = f.readlines()
        
    # chunk lines 
    data = [line.split() for line in data_incoming]
    
    # stop at end of first model
    model_first = []
    for line in data:
        if line[IDX_ATOM] != 'ENDMDL':
            model_first.append(line)
        else:
            break
        
    # get only ATOM records
    model_first = [line for line in model_first if line[IDX_ATOM] == 'ATOM']
    
    # get only specific chains
    model_first = [line for line in model_first if line[IDX_CHAIN] == chain]
    
    # strip down to specific residues
    DATA_TYR = [line for line in model_first if line[IDX_AA] == 'TYR']
    DATA_TRP = [line for line in model_first if line[IDX_AA] == 'TRP']
    
    # strip down to specific atoms using regex
    DATA_TYR = [line for line in DATA_TYR if search(ATOMS_TYR, line[IDX_ATM_LABEL]) != None]
    DATA_TRP = [line for line in DATA_TRP if search(ATOMS_TRP, line[IDX_ATM_LABEL]) != None]
    
    # sort along residue position prior to groupby
    TYR_sorted = sorted(DATA_TYR, key=itemgetter(5)) 
    TRP_sorted = sorted(DATA_TRP, key=itemgetter(5))
    
    # groupby
    TYR_grouped = [list(group) for _, group in groupby(TYR_sorted, lambda x: x[5])]
    TRP_grouped = [list(group) for _, group in groupby(TRP_sorted, lambda x: x[5])]
    
    # get midpoints
    mp_TYR = [groupby_to_midpoint(nl) for nl in TYR_grouped]
    mp_TRP = [groupby_to_midpoint(nl) for nl in TRP_grouped]
    
    # find nearest neighbours
    mp_ALL = mp_TYR + mp_TRP
    
    # iteratively get neighbors
    neighbors = []
    for i in mp_ALL:
        for j in mp_ALL:
            euclidean_dist = norm(i[1] - j[1])
            if euclidean_dist == 0:
                continue
            elif euclidean_dist <= cutoff:
                neighbors.append((i[0], j[0]))
            else:
                continue
    
    # sort neighbors
    neighbors = [tuple(sorted(n)) for n in neighbors]
    
    # drop duplicates
    neighbors = list(set(neighbors))
    
    return neighbors
