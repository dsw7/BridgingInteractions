"""
dsw7@sfu.ca
-----------


What this script does:
======================
Step 1. Import PDB file
Step 2. Get bridge
Step 3. Get metal
Step 4. Get surface
Step 5. Check if all conditions are met
Step 6. Export metal / bridge SD / protein surface coordinate data
Step 7. Destroy PDB file


What does the output look like?
===============================
z84tqfr 5agk SD ATOM 2478 SD MET A 614 11.556 -21.347 37.618 1.00 50.77 S 
z84tqfr 5agk MT HETATM 6728 FE HEM A 750 15.724 4.541 24.638 1.00 24.05 FE 
z84tqfr 5agk SF A 614 CE 12.6590003967 -21.5559997559 39.0209999084 
z84tqfr 5agk BR PHE551 TRP553 
t8ee1jp 5agk SD ATOM 1051 SD MET A 438 2.569 7.652 6.208 1.00 39.30 S <--- the methionine SD
t8ee1jp 5agk MT HETATM 6728 FE HEM A 750 15.724 4.541 24.638 1.00 24.05 FE <--- the metal
t8ee1jp 5agk SF A 415 SG 17.4610004425 5.96299982071 25.716999054 <--- the protein surface
t8ee1jp 5agk BR PHE472 PHE462 <--- the bridge members
   ^
   |
Unique key for grouping surface / bridge / metal data (apply groupby operation here)


Installation notes:
===================
This script requires an installation of Numpy for Python2.7 as we are
using the PyMOL command line to execute the script.
    $ python2.7 -m pip install numpy
Which created a Numpy installation in PyMOL's Python2.7 import path.
The import path can be found as follows:
    $ python2.7
    >>> import numpy
    >>> print numpy.__path__

Note that you may need to manually install NetworkX for Python 2.7:
    $ python2.7 -m pip install networkx


How to use:
===========
(1) Start PyMOL
(2) Type the following:
    PyMOL> run /Users/path/to/scalene_main.py

run /Users/davidweber/Desktop/scalene-triangle/scalene_main.py


BENCHMARKING
============

1a65 data {
    OLD
    ---
    Bridges:
    TYR108 : MET9 : TRP65
    PHE245 : MET197 : TYR195
    
    NEW
    ---
    Bridges:
    'MET197', 'PHE245', 'TYR195'
    
    Note that TYR108 : MET9 : TRP65 is not in the new data because
    NetworkX returns {'MET124', 'MET9', 'PHE106', 'PHE97', 'TRP65', 'TYR108'}
    as disconnected unit (this is not a true 2-bridge)
}
  
  
1aet data {
    OLD
    ---
    Bridges:
    PHE158 : MET172 : PHE202
    PHE202 : MET231 : TRP211

    NEW
    ---
    []
    
    Why?
    NetworkX algorithm returns a connected component:
    {'MET163', 'MET172', 'MET231', 'PHE158', 'PHE202', 'TRP211'}
}


8nse data {
    OLD
    ---
    Bridges: 
    PHE233 : MET209 : PHE243
    PHE429 : MET466 : TRP447
    
    NEW
    ---
    ['PHE233', 'MET209', 'PHE243']
    
    Why?
    NetworkX algorithm returns the connected components:
    [{'MET341', 'PHE355', 'PHE475', 'TYR477'},
     {'MET430', 'MET466', 'PHE429', 'TRP447'},
     {'MET385', 'TRP324'},
     {'MET209', 'PHE233', 'PHE243'}] <-- only this one is a true 2-bridge 
}

"""

from inspect                import stack 
from os                     import path
# the second item in stack[0] is the filename with path
# cannot use sys.path.append(os.path.join(os.getcwd(), libs))) here
filepath_inspected = stack()[0][1]
root = path.dirname(filepath_inspected)
path_to_libs = path.join(root, 'libs')
from sys                    import path as PATH; PATH.append(path_to_libs)
from itertools              import chain
from PDB_filegetter         import PDBFile       # need for PDB file I/O
from get_metals             import get_metals    # need for getting metals
from ma_lowlevel            import met_aromatic  # need for getting bridges
from pymol                  import cmd, stored   # pymol API - need for solvent exposed surface coordinates
from time                   import sleep
from get_nearest_coordinate import get_closest   # function for closest surface coordinates
from numpy                  import array         # need to pass ndarray into get_closest()
from string                 import ascii_lowercase
from random                 import choice, randint, shuffle
import networkx             as nx                # need for getting bridges
import os



# some functions needed for this study
# ------------------------------------
def generate_random_key(n = 5):
    """
    Generates a pseudorandom alphanumeric string
    with 2 integers and n number of chars. Need
    for mapping unique keys to each A : Met : B.
    i.e. for:
        yfmiz80: {A1 : MET : B1}
        2rp0axw: {A2 : MET : B2}
    
    Benchmarking:
    n  Number of iterations  Number of unique combinations
    ------------------------------------------------------
    0  100000                100
    1  100000                7800
    2  100000                88650
    3  100000                99727
    4  100000                99990
    5  100000                100000
    6  100000                100000
    6  1000000               999999
    ------------------------------------------------------
    """
    int_A = str(randint(0, 9))
    int_B = str(randint(0, 9))
    chars = [choice(ascii_lowercase) for c in range(0, n)]
    seq = [int_A, int_B] + chars
    shuffle(seq)
    return ''.join(seq)



# create basic text file for loading results
# ------------------------------------------
filepath_results = os.path.join(root, 'results_scalene.txt')
if os.path.exists(filepath_results) == False:
    open(filepath_results, 'w').close()
    


# Met-aromatic parameters
# -----------------------
CHAIN = 'A'       # input a chain of interest here
CUTOFF = 6.0      # input some cutoff vector v norm
ANGLE = 109.5     # input some cutoff vector a / vector v and vector g / vector v angle
MODEL = 'cp'      # 'cp' or 'rm' -> Cross Product or Rodrigues' Method lone pair interpolation methods
args = (CHAIN, CUTOFF, ANGLE, MODEL)



# pymol literals
# --------------
tmpObj = "__tmp"             # create a pymol temporary object
SOLVENT_EXPOSED_CUTOFF = 2.5 # Angstroms squared
ITER_STATE_EXP = "stored.tmp_dict[(chain, resv, name, x, y, z)] = 1"



# get PDB files known to contain both a metal and a bridge
# ========================================================
with open(path.join(root, 'delim.txt')) as f: 
    delim = f.readlines()  # delim.txt from ~/legacy/get_delimtxt.py
delim = [d.strip() for d in delim]  # remove carriage return

"""
We had a stall at idx 368
We had a stall at idx 659
Both these stalls occurred due to network issues.

Had to rerun the program modifying the following:
delim = delim[368:]  # first stall
delim = delim[368 + 659:]  # second stall


I found a total of 11 duplicates in the dataset using the following code:
==================================================
    from itertools   import groupby
    from collections import Counter

    # check for duplicate data
    with open('results_scalene.txt', 'r') as f:
        input_data = f.readlines()

    # dict style nested list pair    
    input_data = [[d[0:7], d[8:-1]] for d in input_data]

    # concatenate unique data into one string
    isolated_data_A = []
    for _, item in groupby(input_data, lambda x: x[0]):
        list_item = list(item)
        a = list_item[0][1]
        b = list_item[1][1]
        c = list_item[2][1]
        d = list_item[3][1]
        isolated_data_A.append(''.join([a, b, c, d]))  # split by $
        
    duplicates = [i for i, count in Counter(isolated_data_A).items() if count > 1]
    # len(duplicates) is 11
==================================================
"""



# ============================================
# iterate over all files in our delimiter list
# ============================================
for idx, CODE in enumerate(delim):  

    # fetch PDB file
    # --------------
    file_pdb = PDBFile(CODE)
    path_to_file = file_pdb.fetch_from_PDB()
    print 'Iteration {}. Analyzing file {}'.format(idx + 1, CODE)
    if path_to_file == 'URLError': 
        # continue to next file if current file was removed from PDB
        print 'URLError'
        continue



    # [1] find all of SC|TI|V|CR|MN|FE|CO|NI|CU|ZN metal coordinates in PDB file
    # ------------------------------------------------------------
    """
    Here we get metal data directly from the PDB.
    Data is of the form:
    ['HETATM', '3762', 'FE1', 'SF4', 'A', '501', '13.056', '-23.824', '51.566', '1.00', '27.63', 'FE']
    ['HETATM', '3763', 'FE2', 'SF4', 'A', '501', '14.589', '-22.474', '49.584', '1.00', '24.66', 'FE']
    ['HETATM', '3764', 'FE3', 'SF4', 'A', '501', '12.059', '-21.191', '50.530', '1.00', '28.44', 'FE']
    ['HETATM', '3765', 'FE4', 'SF4', 'A', '501', '11.954', '-23.869', '48.715', '1.00', '22.39', 'FE']
    """
    metals = get_metals(path_to_file)
    if metals == []: 
        # next iteration if no metal
        file_pdb.clear()
        print 'No metal.'
        continue



    # [2] find all bridge SD coordinates in PDB file
    # ----------------------------------------------
    try:
        interactions, interacting_SD_data = met_aromatic(path_to_file, *args)
    except Exception as exception:
        # some files are buggy (string crashes, etc)
        # continue over buggy files
        file_pdb.clear()
        print exception
        continue

    """
    interactions list from MetAromatic algorithm is of form:
    ['PHE', '462', 'MET', '438', 5.421154443474194, 157.59264543338077, 91.49943982269996]
    ['PHE', '462', 'MET', '438', 5.583273882768066, 145.2269193938586, 103.53009980437312]
    ['PHE', '462', 'MET', '438', 5.148623845261956, 150.13142075786985, 97.12727156477301]
    ['PHE', '462', 'MET', '438', 5.200337513085088, 162.18691268481302, 87.87862427476162]
    ['PHE', '472', 'MET', '438', 4.248106725354249, 134.17538076546128, 24.771342373699675]
    ['PHE', '472', 'MET', '438', 4.0432115020117365, 130.15630271937837, 26.562880376405165]
    ['PHE', '472', 'MET', '438', 5.0805907382901845, 124.66967377031006, 27.81701176453635]
    ['PHE', '472', 'MET', '438', 5.397525567331757, 131.7643830563905, 22.892236120995378]

    We narrow down to:
    ('PHE', '462', 'MET', '438')
    ('PHE', '462', 'MET', '438')
    ('PHE', '462', 'MET', '438')
    ('PHE', '462', 'MET', '438')
    ('PHE', '472', 'MET', '438')
    ('PHE', '472', 'MET', '438')
    ('PHE', '472', 'MET', '438')
    ('PHE', '472', 'MET', '438')

    And remove duplicates:
    ('PHE', '462', 'MET', '438')
    ('PHE', '472', 'MET', '438')
    """
    pairs = [tuple(i[0:4]) for i in interactions]
    pairs = list(set(pairs))

    """
    Next compress data to prepare for nx.Graph()
    ('PHE', '462', 'MET', '438')
    ('PHE', '472', 'MET', '438')
    To:
    ('PHE462', 'MET438')
    ('PHE472', 'MET438')
    """
    pairs = [('{}{}'.format(i[0], i[1]), '{}{}'.format(i[2], i[3])) for i in pairs]

    """
    We can now see that the two example entries above
    share a common node: MET438. We can use NetworkX
    to find 2-bridges of form:
        PHE462 : MET438 : PHE472

    Which is in the form:
        set(['PHE462', 'MET438', 'PHE472'])
    """
    G1 = nx.Graph()
    G1.add_edges_from(pairs)
    connected_components = list(nx.connected_components(G1))  # i.e. A-B-C, D-E-F-G
    bridges = [i for i in connected_components if len(i) == 3]  # set(['MET163', 'PHE141', 'MET160'])
    bridges = [list(i) for i in bridges]  # cast to list as .pop() (used downstream) is polymorphic

    """
    Now I cast sets of form:
        set(['PHE462', 'MET438', 'PHE472'])
    To dictionary entries of form:
        {'MET438': [['PHE462', 'PHE472'], ['ATOM', '...', 'SD', 'MET', 'A', '438', 'x', 'y', 'z', ...]]}
    Note that I have also added raw SD data to these dictionary entries
    """
    dict_bridges = dict()
    for item in bridges:
        idx_brdg = [i for i, s in enumerate(item) if 'MET' in s]  # first get index of MET entry
        key = item.pop(idx_brdg[0])  # remove the MET entry to generate a key

        # here I also append raw SD data alongside the bridging aromatics
        for SD_data in interacting_SD_data:
            if 'MET' + SD_data[5] == key: # compare 'MET360' to '360'
                additional_SD_data = SD_data

        # the value (in key:value) is the remaining list after the MET entry was .pop()'ed
        if 'MET' not in '{}{}'.format(*item):
            # MET163 ['PHE141', 'MET160'] ['ATOM', '1051'... screen out MET:A:MET data
            dict_bridges[key] = [item, additional_SD_data]
        else:
            continue

    # next iteration if no bridges in protein
    if not any(dict_bridges):
        file_pdb.clear()
        print 'No bridge.'
        continue

   

    # [3] find all protein surface coordinates in file based on solvent exposure
    # --------------------------------------------------------------------------
    """
    Here we obtain all surface coordinates using the PyMOL API:
    Data is of the form:
    ('A', 11, 'C', 8.720000267028809, -8.883999824523926, 27.06999969482422)
    ('A', 11, 'CA', 8.746999740600586, -9.67300033569336, 25.760000228881836)
    ('A', 11, 'CG', 9.218000411987305, -9.347000122070312, 23.23699951171875)
    ('A', 11, 'CZ', 7.163000106811523, -8.295999526977539, 20.868000030517578)
    ('A', 11, 'N', 7.614999771118164, -10.586999893188477, 25.663000106811523)
    ('A', 11, 'NH1', 7.546000003814697, -9.354999542236328, 20.166000366210938)
    ('A', 11, 'NH2', 5.9730000495910645, -7.745999813079834, 20.645999908447266)
    ('A', 11, 'O', 7.65500020980835, -8.621999740600586, 27.628999710083008)
    ('A', 12, 'CA', 9.998000144958496, -7.690999984741211, 28.7810001373291)
    (... 
    """

    # rename to .pdb such that I don't have to cmd.fetch(CODE) (reimport from PDB)
    head, tail = os.path.split(path_to_file)
    output_path = os.path.join(head, tail[3:7] + '.pdb')
    os.rename(path_to_file, output_path)  # .ent -> .pdb

    # see https://github.com/dsw7/BridgingInteractions/tree/master/scalene-triangle/pymol-get-surface-example
    cmd.load(output_path)
    cmd.create(tmpObj, "({} and polymer and chain {}) and not resn HOH".format(CODE, CHAIN))
    cmd.set("dot_solvent")
    cmd.get_area(selection=tmpObj, load_b=1)
    cmd.remove(tmpObj + " and b < " + str(SOLVENT_EXPOSED_CUTOFF))
    cmd.show(selection=tmpObj, representation="dots") # show the exposed atoms
    stored.tmp_dict = {}
    cmd.iterate_state(state=-1, selection=tmpObj, expression=ITER_STATE_EXP)
    exposed = stored.tmp_dict.keys()
    exposed.sort()
    cmd.delete(tmpObj)
    cmd.delete('all')
    os.rename(output_path, path_to_file)  # .pdb -> .ent such that file_pdb.clear() can del dir



    # compute closest surface coordinate if protein ends up having both bridge and metal
    # --------------------------------------------------------------------------
    try:
        # must put this in a try / catch block as SOME R3 coords are float dtype but 
        # still string crashed: 2m48 for example, -27.595-106.521-108.313
        for _, value in dict_bridges.iteritems():
            for m in metals:

                # need to convert to float dtype before passing into get_closest()
                array_b_coord = array(value[1][6:9]).astype(float)  # bridge SD
                array_m_coord = array(m[6:9]).astype(float)         # metal
                coord_surf = get_closest(exposed, array_b_coord, array_m_coord)

                # need to convert back to str dtype before writing to text
                coord_surf = [str(e) for e in coord_surf[1]]

                # need to get unique key for each iteration
                key_rowlabel = generate_random_key()

                # outgoing data
                SD_OUT = '{} {} SD {} \n'.format(key_rowlabel, CODE, ' '.join(value[1]))
                MT_OUT = '{} {} MT {} \n'.format(key_rowlabel, CODE, ' '.join(m))
                SF_OUT = '{} {} SF {} \n'.format(key_rowlabel, CODE, ' '.join(coord_surf))
                BRIDGE = '{} {} BR {} \n'.format(key_rowlabel, CODE, ' '.join(value[0]))

                # have to open/write/close per iteration
                # otherwise data is loaded to a buffer and appended to .txt
                # only after the file has been closed
                with open(filepath_results, 'a') as txtfile_results:
                    txtfile_results.write(SD_OUT)
                    txtfile_results.write(MT_OUT)
                    txtfile_results.write(SF_OUT)
                    txtfile_results.write(BRIDGE)

    except Exception as exception:
        file_pdb.clear()
        print exception
        continue



    # end of analysis
    # ---------------
    file_pdb.clear()

print 'Done!'
    

