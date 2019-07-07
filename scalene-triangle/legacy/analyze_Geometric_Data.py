import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

# ----------------------------------------------------
# all functions and other definitions

areas_cutoff = 150.0   
vec_a_cutoff = 10.0
vec_b_cutoff = 100.0
vec_c_cutoff = 100.0

bin_choice = 40
saveBoolean = False

def area(p1, p2, p3):
    # a function for calculating area between three coordinates
    cross = np.cross(p2 - p1, p3 - p1)
    return 0.5 * np.linalg.norm(cross)

# ----------------------------------------------------

# get data
root = '/Volumes/MSC/DATABASES/2018 Bridging Databases/GeometricData_compressed.csv'
df = pd.read_csv(root).drop('Unnamed: 0', axis=1)

# start off with a count of how many files were analyzed:
pre_count = df.shape[0]
ddr_count = df[['PDB CODE']].drop_duplicates().shape[0]

print('Size of input:', pre_count)
print('Number of PDB files:', ddr_count)

# next we split data into two streams: TYR/TRP bridges and non-TYR/TRP bridges
df_PHE = df[(df['ARO A'] == 'PHE') | (df['ARO B'] == 'PHE')]
df_RAD = df[(df['ARO A'] != 'PHE') & (df['ARO B'] != 'PHE')]

if df_PHE.shape[0] + df_RAD.shape[0] != pre_count:
    sys.exit('Some sort of counting error occurring!')
    
# ----------------------------------------------------
# PHE containing bridge stream

# strip down to parts
df_SDS_PHE = df_PHE[['SD_x', 'SD_y', 'SD_z']]  # all SD coordinates
df_MTL_PHE = df_PHE[['Z_x', 'Z_y', 'Z_z']]  # all metal coordinates
df_SUR_PHE = df_PHE[['S_x', 'S_y', 'S_z']]  # all surface coordinates

array_SDS_PHE = df_SDS_PHE.values
array_MTL_PHE = df_MTL_PHE.values
array_SUR_PHE = df_SUR_PHE.values

"""
Our definitions (may change order):

VEC_A = surface coordinate - met SD coordinate
VEC_B = metal coordinate - met SD coordinate
VEC_C = surface coordinate - metal coordinate
"""

# get distances between various points
VEC_A = []
VEC_B = []
VEC_C = []

for index in range(0, df_PHE.shape[0]):
    VEC_A.append(np.linalg.norm(array_SUR_PHE[index] - array_SDS_PHE[index]))
    VEC_B.append(np.linalg.norm(array_MTL_PHE[index] - array_SDS_PHE[index]))
    VEC_C.append(np.linalg.norm(array_SUR_PHE[index] - array_MTL_PHE[index]))
 
# get area of scalene triangles  
AREAS = []
for index in range(0, df_PHE.shape[0]):
    ret_val = area(array_SDS_PHE[index], array_MTL_PHE[index], array_SUR_PHE[index])
    AREAS.append(ret_val)
    
# cut out data above a certain threshold as it results in poor axes dimensions
VEC_A = [i for i in VEC_A if i <= vec_a_cutoff]
VEC_B = [i for i in VEC_B if i <= vec_b_cutoff]
VEC_C = [i for i in VEC_C if i <= vec_c_cutoff]
AREAS = [i for i in AREAS if i <= areas_cutoff]
    
f, axarr = plt.subplots(2, 2, figsize=(10, 10))
axarr[0, 0].hist(AREAS, color='b', edgecolor='k', bins=bin_choice)
axarr[0, 0].set_title('AREAS', size=16)
axarr[0, 0].set_xlabel(r'Area / $\AA^2$', size=14)
axarr[0, 1].hist(VEC_A, color='r', edgecolor='k', bins=bin_choice)
axarr[0, 1].set_title('S / SD', size=16)
axarr[0, 1].set_xlabel(r'Length / $\AA$', size=14)
axarr[1, 0].hist(VEC_B, color='r', edgecolor='k', bins=bin_choice)
axarr[1, 0].set_title('M / SD', size=16)
axarr[1, 0].set_xlabel(r'Length / $\AA$', size=14)
axarr[1, 1].hist(VEC_C, color='r', edgecolor='k', bins=bin_choice)
axarr[1, 1].set_title('S / M', size=16)
axarr[1, 1].set_xlabel(r'Length / $\AA$', size=14)

for a, b in zip([0, 0, 1, 1], [0, 1, 0, 1]):
    axarr[a, b].spines['right'].set_visible(False)
    axarr[a, b].spines['top'].set_visible(False)
    axarr[a, b].set_ylabel('Frequency', size=14)
    axarr[a, b].set_ylim(0, 5500)
    
f.subplots_adjust(hspace=0.3, wspace=0.2)  # adjust subplot dists from each other
if saveBoolean:
    plt.savefig('triangle_geometry_PHE.png', dpi=1000, bbox_inches='tight')
plt.show()
    


# ----------------------------------------------------
# TYR / TRP containing bridge stream

# strip down to parts
df_SDS_RAD = df_RAD[['SD_x', 'SD_y', 'SD_z']]  # all SD coordinates
df_MTL_RAD = df_RAD[['Z_x', 'Z_y', 'Z_z']]  # all metal coordinates
df_SUR_RAD = df_RAD[['S_x', 'S_y', 'S_z']]  # all surface coordinates

array_SDS_RAD = df_SDS_RAD.values
array_MTL_RAD = df_MTL_RAD.values
array_SUR_RAD = df_SUR_RAD.values

"""
Our definitions (may change order):

VEC_A = surface coordinate - met SD coordinate
VEC_B = metal coordinate - met SD coordinate
VEC_C = surface coordinate - metal coordinate
"""

# get distances between various points
VEC_A = []
VEC_B = []
VEC_C = []

for index in range(0, df_RAD.shape[0]):
    VEC_A.append(np.linalg.norm(array_SUR_RAD[index] - array_SDS_RAD[index]))
    VEC_B.append(np.linalg.norm(array_MTL_RAD[index] - array_SDS_RAD[index]))
    VEC_C.append(np.linalg.norm(array_SUR_RAD[index] - array_MTL_RAD[index]))

# get area of scalene triangles  
AREAS = []
for index in range(0, df_RAD.shape[0]):
    ret_val = area(array_SDS_RAD[index], array_MTL_RAD[index], array_SUR_RAD[index])
    AREAS.append(ret_val)
    
# cut out data above a certain threshold as it results in poor axes dimensions
VEC_A = [i for i in VEC_A if i <= vec_a_cutoff]
VEC_B = [i for i in VEC_B if i <= vec_b_cutoff]
VEC_C = [i for i in VEC_C if i <= vec_c_cutoff]
AREAS = [i for i in AREAS if i <= areas_cutoff]
    
f, axarr = plt.subplots(2, 2, figsize=(10, 10))
axarr[0, 0].hist(AREAS, color='b', edgecolor='k', bins=bin_choice)
axarr[0, 0].set_title('AREAS', size=16)
axarr[0, 0].set_xlabel(r'Area / $\AA^2$', size=14)
axarr[0, 1].hist(VEC_A, color='r', edgecolor='k', bins=bin_choice)
axarr[0, 1].set_title('S / SD', size=16)
axarr[0, 1].set_xlabel(r'Length / $\AA$', size=14)
axarr[1, 0].hist(VEC_B, color='r', edgecolor='k', bins=bin_choice)
axarr[1, 0].set_title('M / SD', size=16)
axarr[1, 0].set_xlabel(r'Length / $\AA$', size=14)
axarr[1, 1].hist(VEC_C, color='r', edgecolor='k', bins=bin_choice)
axarr[1, 1].set_title('S / M', size=16)
axarr[1, 1].set_xlabel(r'Length / $\AA$', size=14)

for a, b in zip([0, 0, 1, 1], [0, 1, 0, 1]):
    axarr[a, b].spines['right'].set_visible(False)
    axarr[a, b].spines['top'].set_visible(False)
    axarr[a, b].set_ylabel('Frequency', size=14)
    axarr[a, b].set_ylim(0, 1750)

f.subplots_adjust(hspace=0.3, wspace=0.2)  # adjust subplot dists from each other
if saveBoolean:
    plt.savefig('triangle_geometry_RAD.png', dpi=1000, bbox_inches='tight')
plt.show()
    





