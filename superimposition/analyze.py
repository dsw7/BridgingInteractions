""" 
Written by David S. Weber
Get pie chart from results_superimposition.txt 
"""

# -----------------------------------------------------------------------------------
# upstream defs

from itertools   import groupby
from operator    import itemgetter
import matplotlib.pyplot as plt

NR_CONDITION = 'NR : {} : {}'
PM_CONDITION = 'PM : {} : {}'
DS_CONDITION = 'DS : {} : {}'
IS_CONDITION = 'IS : {} : {}'
UK_CONDITION = 'UC'  # unknown condition

def compare_networks(chain_list, bridge_list):
    """
    Compares the relationship between a 
    set of disconnected graph components
    
    list(nx.connected_components(G1))
    
    Parameters
    ----------
    chain_list : a list of closely spaced aromatic chains
    bridge_list : a list A / B bridges where we have A::MET::B
            
    chain_list, bridge_list are of type list(nx.connected_components(G))
    where G is an nx.Graph() object
    Returns
    -------
    A list of any of NR_CONDITION, PM_CONDITION,
                     DS_CONDITION, IS_CONDITION,
                     UK_CONDITION
    """
    # beautiful, syntactic Python logic
    list_cond = []
    for aromatic_chain in chain_list:
        for bridge in bridge_list:    
            if aromatic_chain & bridge == set():                 # NR
                list_cond.append(NR_CONDITION.format(bridge, aromatic_chain))
            elif aromatic_chain - bridge == set():               # 2R
                list_cond.append(PM_CONDITION.format(bridge, aromatic_chain))
            elif aromatic_chain & bridge == bridge:              # DS
                list_cond.append(DS_CONDITION.format(bridge, aromatic_chain))
            elif next(iter(aromatic_chain & bridge)) in bridge:  # IS
                list_cond.append(IS_CONDITION.format(bridge, aromatic_chain))
            else:
                list_cond.append(UK_CONDITION)
    return list_cond

# -----------------------------------------------------------------------------------
# process data and plot pie chart

# filename = r'/Volumes/MSC/DATABASES/2018 Bridging Databases/results_updated_with_AroMetX.txt'
filename = 'results_updated_with_AroMetX.txt'
with open(filename) as f: data = f.readlines()
        
data_sorted   = sorted(data, key=itemgetter(0)) 
data_sorted   = [d.split(':') for d in data]
grouped       = [list(group) for _, group in groupby(data_sorted, lambda x: x[0])]
all_codes     = [keys for keys, _ in groupby(data_sorted, lambda x: x[0])]

number_of_proteins   = 146502             # all proteins in PDB on 23 Nov 2018
analyzed_proteins    = len(all_codes)     # proteins that have both a chain and bridge
discarded_proteins   = number_of_proteins - analyzed_proteins

set_comparisons = []
for PDBFile in grouped:
    chains =  [set(line[2].split(',')[:-1]) for line in PDBFile if 'Chains' in line[1]]
    bridges = [set(line[2].split(',')[:-1]) for line in PDBFile if 'Bridges' in line[1]]
    for comparisons in compare_networks(chains, bridges):
        set_comparisons.append(comparisons)
    
membership = [d[0:2] for d in set_comparisons]
counts_NR = membership.count('NR')
counts_PM = membership.count('PM')
counts_IS = membership.count('IS')
counts_DS = membership.count('DS')
    
sizes1 = [analyzed_proteins, discarded_proteins]
sizes2 = [counts_NR, counts_PM, counts_IS, counts_DS]

def formatter1(pct):
    return '{:0.0f}\n({:0.1f}%)'.format(pct * sum(sizes1) / 100, pct)

dimensions = 1.5
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(dimensions * 5.0, dimensions * 2.5))

labels1 = ('A', 'NA')
labels2 = ('NR', '2R', 'IS', 'DS')

explode1 = (0.0, 0.075)
explode2 = (0.0, 0.075, 0.15, 0.225)

ax1.pie(sizes1, 
        labels=labels1, 
        #autopct='%1.1f%%',
        autopct=formatter1,
        explode=explode1, 
        textprops={'fontsize': 12})
        
ax2.pie(sizes2, 
        labels=labels2,  
        autopct='%1.1f%%', 
        explode=explode2, 
        startangle=-30, 
        pctdistance=0.75, 
        textprops={'fontsize': 12},
        labeldistance=1.15)

ax1.axis('equal')
ax2.axis('equal')

ax1.text(-0.05, 
         -1.25, 
         'PDB Files', 
         ha='center', 
         size=16, 
         va='center', 
         bbox=dict(boxstyle='round', facecolor='white', edgecolor='k'))
         
ax2.text(0.0, 
         1.2, 
         'Bridges',
         ha='center', 
         size=16,
         va='center',  
         bbox=dict(boxstyle='round', facecolor='white', edgecolor='k'))

#plt.savefig('bridge_counts.png', bbox_inches='tight', dpi=1000)
plt.show()


# print some stats to console
print('Total number of proteins analyzed: {}'.format(number_of_proteins))
print('Proteins containing both a chain and a bridge: {}'.format(analyzed_proteins))
print('Number of bridge / chain pairs compared: {}'.format(len(membership)))
print('NR counts: {}'.format(counts_NR))
print('DS counts: {}'.format(counts_DS))
print('IS counts: {}'.format(counts_IS))
print('2R counts: {}'.format(counts_PM))

    
