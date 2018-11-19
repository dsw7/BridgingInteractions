"""
Written by David S. Weber
dsw7@sfu.ca

This script collects bridge anatomy counts from the SQL database. This script
also strips out periodic 2-bridges of form ...A : Met : B : Met : C...

Flow scheme:  
    main.py ---> [SQL DATABASE] ---> counts.py ---> 3x3.png heatmap
"""

import pyodbc
import networkx as nx  
import matplotlib.pyplot as plt

from platform     import node
from itertools    import groupby, chain
from operator     import itemgetter

# ------------------------------------------------------------------------------
# SQL

# connect to SQL database  
# -----------------------     
DRIVER   = 'SQL Server'
SERVER   = '{}\SQLEXPRESS'.format(node())
DATABASE = 'Bridging'
TRUSTED  = 'yes'
str_conn = """
           Driver={};
           Server={};
           Database={};
           Trusted_Connection={}
           """.format(DRIVER, SERVER, DATABASE, TRUSTED)            


# read in all data from SQL database
# ----------------------------------
connection = pyodbc.connect(str_conn)
cursor = connection.cursor()
cursor.execute('SELECT * FROM BridgingResults')
rows = cursor.fetchall()  
connection.close()


# ------------------------------------------------------------------------------
# get bridges of form A : MET : B


# group by PDB code
# -----------------
rows_sorted  = sorted(rows, key=itemgetter(0))
rows_grouped = [(a, list(b)) for a, b in groupby(rows_sorted, lambda x: x[0])]


# apply NetworkX methods for isolating 2-bridges from SQL database
# ----------------------------------------------------------------
bridges = []
for _, GB_OBJ in rows_grouped:
    pairs = [('{}{}'.format(i[1], i[2]), '{}{}'.format(i[3], i[4])) for i in GB_OBJ]
    pairs = list(set(pairs))  # remove identical lines due to multiple vectors v
    G = nx.Graph()
    G.add_edges_from(pairs)
    for disconnects in list(nx.connected_components(G)):
        if len(disconnects) == 3: bridges.append(list(disconnects))


# remove MET data
# ---------------
for nx_sets in bridges:
    for entry in nx_sets.copy():     # work on a copy of a set due to runtime error
        if 'MET' in entry: nx_sets.remove(entry)
 

# remove inverse 2-bridges (Met : Aro : Met) and prep for imshow rendering
# ------------------------------------------------------------------------
bridges = list(filter(lambda entry: len(entry) > 1, bridges))
bridging_pairs = ['{}-{}'.format(i[0][0:3], i[1][0:3]) for i in bridges]


# ------------------------------------------------------------------------------
# heatmap


# diagonals
ENTRY_1 = 'PHE-PHE'
ENTRY_2 = 'TYR-TYR'
ENTRY_3 = 'TRP-TRP'

# lowers
ENTRY_4 = 'PHE-TYR'
ENTRY_5 = 'TYR-PHE'
ENTRY_6 = 'PHE-TRP'
ENTRY_7 = 'TRP-PHE'
ENTRY_8 = 'TYR-TRP'
ENTRY_9 = 'TRP-TYR'

# load diagonals
a0_0 = bridging_pairs.count(ENTRY_1)
a1_1 = bridging_pairs.count(ENTRY_2)
a2_2 = bridging_pairs.count(ENTRY_3)

# load lowers
a0_1 = bridging_pairs.count(ENTRY_4) + bridging_pairs.count(ENTRY_5)
a0_2 = bridging_pairs.count(ENTRY_6) + bridging_pairs.count(ENTRY_7)
a1_2 = bridging_pairs.count(ENTRY_8) + bridging_pairs.count(ENTRY_9)

# format imshow matrix
to_imshow = [[a0_0, a0_1, a0_2],
             [   0, a1_1, a1_2],
             [   0,    0, a2_2]]

# get net number of counted bridges          
total_bridges = sum(list(chain(*to_imshow)))             
     
plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(8, 8))

# add white spaces between pixels
pos = list(range(0, 3))
for item in pos:
    plt.axhline(item - 0.5, lw=5, c='white')
    plt.axvline(item - 0.5, lw=5, c='white')
 
# set tick labels
LIST_LABELS = ('PHE', 'TYR', 'TRP')
plt.tick_params(axis='both', which='both', length=0)
plt.xticks(pos, LIST_LABELS, rotation='vertical', size=28)  
plt.yticks(pos, LIST_LABELS, rotation='horizontal', size=28)
    
# annotate individual pixels
for i in range(len(to_imshow)):
    for j in range(len(to_imshow[i])):
        if (i, j) not in [(1, 0), (2, 0), (2, 1)]: 
            perc = str(round((to_imshow[i][j] / total_bridges), 2))
            plt.text(j, i, str(to_imshow[i][j]) + '\n' + perc, 
                             ha="center", va="center", color='r', size=25)     
    
# get plot
plt.imshow(to_imshow, cmap='binary', aspect='equal', origin='lower')
plt.tight_layout()

# hide the border
plt.gca().set_frame_on(False)

# save the heatmap
bool_save = False
if bool_save: plt.savefig('3x3.png', dpi=1000, bbox_inches='tight')

plt.show()


