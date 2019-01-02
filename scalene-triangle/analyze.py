"""
dsw7@sfu.ca
This script processes data in results_scalene.txt into
three histograms corresponding to scalene triangles of
vertices:

PHE : MET(SD) : PHE / METAL / SURFACE
PHE : MET(SD) : ARO / METAL / SURFACE
ARO : MET(SD) : ARO / METAL / SURFACE
"""


from itertools import groupby
from numpy     import array, linalg, cross
import matplotlib.pyplot as plt

size = 7   # plot size in inches
adj = 0.5  # float to pass into subplot adjust
save_plot = False



# get data
# =========================== 
with open('results_scalene.txt', 'r') as f:
    data = f.readlines()
data = [d.split() for d in data]  



# print number of PDB entries studied to console
# =========================== 
codes_PDB = [item[1] for item in data]
codes_PDB = list(set(codes_PDB))
output = 'Number of PDB codes studied: {}'.format(len(codes_PDB))
print(output)
    
    
    
# compute triangle geometries
# ===========================    
"""
Here the data is grouped by the unique alphanumeric key we assigned in scalene_main.py

['pew8zy6', '5xfq', 'SD', 'ATOM', '1351', 'SD', 'MET', 'A', '201', '48.571', '12.902', '-1.699', '1.00', '44.27', 'S']
['pew8zy6', '5xfq', 'MT', 'HETATM', '5546', 'ZN', 'ZN', 'A', '403', '43.771', '14.564', '-7.031', '1.00', '40.94', 'ZN']
['pew8zy6', '5xfq', 'SF', 'A', '195', 'O', '52.827999115', '13.4440002441', '-3.04399991035']
['pew8zy6', '5xfq', 'BR', 'TYR188', 'TRP210']
"""

dict_master = dict()
for key, grouper_obj in groupby(data, lambda x: x[0]):
    # cast iterator to list
    # list_grouper = [list(g) for g in grouper_obj]
    list_grouper = list(grouper_obj)
    
    # isolate the data of interest
    SD = array(list_grouper[0][9:12]).astype(float)
    MT = array(list_grouper[1][9:12]).astype(float)
    SF = array(list_grouper[2][6:9]).astype(float)
    BR = '{}-{}'.format(list_grouper[3][3][0:3], list_grouper[3][4][0:3])
    
    u = MT - SD
    v = SF - SD

    # get triangle faces
    SD_MT = linalg.norm(u) # MT - SD
    SF_SD = linalg.norm(v) # SF - SD
    SF_MT = linalg.norm(SF - MT)
    
    # triangle area
    AREA = 0.5 * linalg.norm(cross(u, v))
    
    # prep an ID key for histogram dict and load dict
    key_for_dict = ':'.join([key, list_grouper[0][1], BR])
    dict_master[key_for_dict] = [SD_MT, SF_MT, SF_SD, AREA]
   
   
   
# separate data
# ===========================    
   
dict_PHE_PHE = dict() 
dict_PHE_ARO = dict() 
dict_ARO_ARO = dict() 

for key, value in dict_master.items():
    AA = key.split(':')[2]
    if AA == 'PHE-PHE':                     # i.e. PHE-PHE bridges
        dict_PHE_PHE[key] = value
    elif 'PHE' in AA and AA != 'PHE-PHE':   # i.e. PHE-TYR, TYR-PHE, etc.
        dict_PHE_ARO[key] = value
    elif 'PHE' not in AA:                   # i.e. TYR-TYR, TYR-TRP, etc.
        dict_ARO_ARO[key] = value
    
len_PHE_PHE = len(dict_PHE_PHE)
len_PHE_ARO = len(dict_PHE_ARO)
len_ARO_ARO = len(dict_ARO_ARO)
len_master  = len(dict_master)

# ensure data is not being lost
if len_PHE_PHE + len_PHE_ARO + len_ARO_ARO != len_master:
    print('Data is being lost.')
    
print('Total number of scalene triangles: {}'.format(len_master))
print('Number of PHE-PHE triangles: {}'.format(len_PHE_PHE))
print('Number of PHE-ARO triangles: {}'.format(len_PHE_ARO))
print('Number of ARO-ARO triangles: {}'.format(len_ARO_ARO))



# plot
# ===========================    

# @DSW: cutoffs and definitions from analyze_Geometric_Data.py (legacy script)
areas_cutoff = 150.0   
vec_a_cutoff = 100.0
vec_b_cutoff = 100.0
vec_c_cutoff = 100.0
bin_choice   = 40

"""
Our definitions (may change order):
VEC_A = surface coordinate - met SD coordinate
VEC_B = metal coordinate - met SD coordinate
VEC_C = surface coordinate - metal coordinate

                    VEC_B  VEC_C  VEC_A  AREAS
dict_master order: [SD_MT, SF_MT, SF_SD, AREA]
"""



# PHE / PHE case
# ***************************

# @DSW: "connect" new data with legacy code
VEC_A = []
VEC_B = []
VEC_C = []
AREAS = []

for _, entry in dict_PHE_PHE.items():
    VEC_A.append(entry[2])  # see above chart for this order
    VEC_B.append(entry[0])
    VEC_C.append(entry[1])
    AREAS.append(entry[3])

# now eliminate data above a certain cutoff for visualization
VEC_A = [i for i in VEC_A if i <= vec_a_cutoff]  # surface / SD
VEC_B = [i for i in VEC_B if i <= vec_b_cutoff]  # metal / SD
VEC_C = [i for i in VEC_C if i <= vec_c_cutoff]  # surface / metal
AREAS = [i for i in AREAS if i <= areas_cutoff]
    
f, axarr = plt.subplots(2, 2, figsize=(size, size))
f.suptitle('PHE : MET : PHE', size=15)
f.subplots_adjust(wspace=adj, hspace=adj)
axarr[0, 0].hist(AREAS, color='b', edgecolor='k', bins=bin_choice)
axarr[0, 0].set_title('AREAS', size=16)
axarr[0, 0].set_xlabel(r'Area / $\AA^2$', size=14)
axarr[0, 1].hist(VEC_A, color='r', edgecolor='k', bins=bin_choice)
axarr[0, 1].set_title('SF / SD', size=16)
axarr[0, 1].set_xlabel(r'Length / $\AA$', size=14)
axarr[1, 0].hist(VEC_B, color='r', edgecolor='k', bins=bin_choice)
axarr[1, 0].set_title('MT / SD', size=16)
axarr[1, 0].set_xlabel(r'Length / $\AA$', size=14)
axarr[1, 1].hist(VEC_C, color='r', edgecolor='k', bins=bin_choice)
axarr[1, 1].set_title('SF / MT', size=16)
axarr[1, 1].set_xlabel(r'Length / $\AA$', size=14)

for a, b in zip([0, 0, 1, 1], [0, 1, 0, 1]):
    axarr[a, b].spines['right'].set_visible(False)
    axarr[a, b].spines['top'].set_visible(False)
    axarr[a, b].set_ylabel('Frequency', size=14)
    axarr[a, b].set_ylim(0, 3500)
    
if save_plot:
    plt.savefig('all_phe.png', dpi=800, bbox_inches='tight')
    
    
    
    
# PHE / ARO case
# ***************************  
    
# @DSW: "connect" new data with legacy code
VEC_A = []
VEC_B = []
VEC_C = []
AREAS = []

for _, entry in dict_PHE_ARO.items():
    VEC_A.append(entry[2])  # see above chart for this order
    VEC_B.append(entry[0])
    VEC_C.append(entry[1])
    AREAS.append(entry[3])

# now eliminate data above a certain cutoff for visualization
VEC_A = [i for i in VEC_A if i <= vec_a_cutoff]
VEC_B = [i for i in VEC_B if i <= vec_b_cutoff]
VEC_C = [i for i in VEC_C if i <= vec_c_cutoff]
AREAS = [i for i in AREAS if i <= areas_cutoff]
    
f, axarr = plt.subplots(2, 2, figsize=(size, size))
f.suptitle('PHE : MET : ARO', size=15)
f.subplots_adjust(wspace=adj, hspace=adj)
axarr[0, 0].hist(AREAS, color='b', edgecolor='k', bins=bin_choice)
axarr[0, 0].set_title('AREAS', size=16)
axarr[0, 0].set_xlabel(r'Area / $\AA^2$', size=14)
axarr[0, 1].hist(VEC_A, color='r', edgecolor='k', bins=bin_choice)
axarr[0, 1].set_title('SF / SD', size=16)
axarr[0, 1].set_xlabel(r'Length / $\AA$', size=14)
axarr[1, 0].hist(VEC_B, color='r', edgecolor='k', bins=bin_choice)
axarr[1, 0].set_title('MT / SD', size=16)
axarr[1, 0].set_xlabel(r'Length / $\AA$', size=14)
axarr[1, 1].hist(VEC_C, color='r', edgecolor='k', bins=bin_choice)
axarr[1, 1].set_title('SF / MT', size=16)
axarr[1, 1].set_xlabel(r'Length / $\AA$', size=14)

for a, b in zip([0, 0, 1, 1], [0, 1, 0, 1]):
    axarr[a, b].spines['right'].set_visible(False)
    axarr[a, b].spines['top'].set_visible(False)
    axarr[a, b].set_ylabel('Frequency', size=14)
    axarr[a, b].set_ylim(0, 4500)
    
if save_plot:
    plt.savefig('phe_aro.png', dpi=800, bbox_inches='tight')    
    
    
    
    
# ARO / ARO case
# ***************************  
    
# @DSW: "connect" new data with legacy code
VEC_A = []
VEC_B = []
VEC_C = []
AREAS = []

for _, entry in dict_ARO_ARO.items():
    VEC_A.append(entry[2])  # see above chart for this order
    VEC_B.append(entry[0])
    VEC_C.append(entry[1])
    AREAS.append(entry[3])

# now eliminate data above a certain cutoff for visualization
VEC_A = [i for i in VEC_A if i <= vec_a_cutoff]
VEC_B = [i for i in VEC_B if i <= vec_b_cutoff]
VEC_C = [i for i in VEC_C if i <= vec_c_cutoff]
AREAS = [i for i in AREAS if i <= areas_cutoff]
   

f, axarr = plt.subplots(2, 2, figsize=(size, size))
f.suptitle('ARO : MET : ARO', size=15)
f.subplots_adjust(wspace=adj, hspace=adj)
axarr[0, 0].hist(AREAS, color='b', edgecolor='k', bins=bin_choice)
axarr[0, 0].set_title('AREAS', size=16)
axarr[0, 0].set_xlabel(r'Area / $\AA^2$', size=14)
axarr[0, 1].hist(VEC_A, color='r', edgecolor='k', bins=bin_choice)
axarr[0, 1].set_title('SF / SD', size=16)
axarr[0, 1].set_xlabel(r'Length / $\AA$', size=14)
axarr[1, 0].hist(VEC_B, color='r', edgecolor='k', bins=bin_choice)
axarr[1, 0].set_title('MT / SD', size=16)
axarr[1, 0].set_xlabel(r'Length / $\AA$', size=14)
axarr[1, 1].hist(VEC_C, color='r', edgecolor='k', bins=bin_choice)
axarr[1, 1].set_title('SF / MT', size=16)
axarr[1, 1].set_xlabel(r'Length / $\AA$', size=14)

for a, b in zip([0, 0, 1, 1], [0, 1, 0, 1]):
    axarr[a, b].spines['right'].set_visible(False)
    axarr[a, b].spines['top'].set_visible(False)
    axarr[a, b].set_ylabel('Frequency', size=14)
    axarr[a, b].set_ylim(0, 2000)

if save_plot:
    plt.savefig('aro_aro.png', dpi=800, bbox_inches='tight')    


    
    

    
   
    

