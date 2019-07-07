"""
dsw7@sfu.ca

Testing output of Winter 2018 scalene triangle study against that of
Summer 2018 study.
"""

from itertools import groupby
from numpy import array, cross, linalg

test_code = '3ii2'

# summer 2018 data
# ===================================================

path_summer_data = '/Volumes/MSC/DATABASES/2018 Bridging Databases/GeometricData_compressed.csv'

with open(path_summer_data) as f:
    data = f.readlines()
    
data = [d.split(',') for d in data]
data = [d for d in data if d[6] == test_code]

"""
,           0
ARO A,      1
ARO A RES,  2
ARO B,      3
ARO B RES,  4
MET RES,    5
PDB CODE,   6
SD_x,       7
SD_y,       8
SD_z,       9
SUR RES,    10
S_x,        11
S_y,        12
S_z,        13
Z_x,        14
Z_y,        15
Z_z         16
"""

print('Summer 2018 data:')
print('-' * 50)
for i in data:
    print()
    brg = 'Bridge: {}{}-{}{}'.format(*i[1:5])
    met  = 'MET: {}'.format(i[5])
    sur  = 'Surface: {}'.format(i[10])
    
    SD = array(i[7:10]).astype(float)
    SF = array(i[11:14]).astype(float)
    MT = array(i[14:17]).astype(float)
    
    SD_SF = 'SD-SF Distance: {}'.format(linalg.norm(SD - SF))
    SD_MT = 'SD-MT Distance: {}'.format(linalg.norm(SD - MT))
    MT_SF = 'MT-SF Distance: {}'.format(linalg.norm(SF - MT))
    
    AREA = 0.5 * linalg.norm(cross(MT - SF, SD - MT))
    area = 'Area: {}'.format(AREA)
    
    print(brg)
    print(met)
    print(sur)
    print(area)
    print(SD_SF)
    print(SD_MT)
    print(MT_SF)
    
    print('*' * 10)
    
# winter 2018 data
# ===================================================

print()
print('Winter 2018 data:')
print('-' * 50)
path_winter_data = '/Users/davidweber/Desktop/scalene-triangle/results_scalene.txt'

with open(path_winter_data) as f:
    data2 = f.readlines()
    
data2 = [d.split() for d in data2]
data2 = [d for d in data2 if d[1] == test_code]

for key, value in groupby(data2, lambda k: k[0]):
    print()
    value = list(value)
    
    brg = 'Bridge: {}-{}'.format(*value[3][3:5])
    met  = 'Methionine: {}'.format(value[0][8])
    sur = 'Surface: {} \ {}'.format(*value[2][4:6])
    
    SD = array(value[0][9:12]).astype(float)
    SF = array(value[2][6:9]).astype(float)
    MT = array(value[1][9:12]).astype(float)
    
    AREA = 0.5 * linalg.norm(cross(MT - SF, SD - MT))
    area = 'Area: {}'.format(AREA)
    
    print(brg)
    print(met)    
    print(sur)
    print(area)
    
    SD_SF = 'SD-SF Distance: {}'.format(linalg.norm(SD - SF))
    SD_MT = 'SD-MT Distance: {}'.format(linalg.norm(SD - MT))
    MT_SF = 'MT-SF Distance: {}'.format(linalg.norm(SF - MT))
    
    print(SD_SF)
    print(SD_MT)
    print(MT_SF)
    
    print('*' * 10)
    
"""
Files tested:

1lko - identical
4dy1 - identical except A:Met:B:Met... data is being removed
3nvl - 

// not in new data - right off the bat I can tell why: not a 2-bridge
Bridge: PHE531-PHE537
MET: 535
Surface: 14.0
Area: 45.62220235877811
**********
Bridge: PHE531-PHE537
MET: 535
Surface: 14.0
Area: 55.35523238052146
**********
Bridge: PHE531-TRP364
MET: 415
Surface: 417.0
Area: 29.124735233719786
**********
Bridge: PHE531-TRP364
MET: 415
Surface: 417.0
Area: 24.658555922878154
**********

{
Bridge: PHE292-TYR310
MET: 288
Surface: 288.0
Area: 21.113220846469698
**********
Bridge: TYR310-PHE292
Methionine: 288
Surface: 288 \ CG
Area: 21.113220846433798
**********
}

{
Bridge: PHE47-PHE503
MET: 43
Surface: 43.0
Area: 3.4164676622511783e-06
**********
Bridge: PHE47-PHE503
Methionine: 43
Surface: 43 \ SD
Area: 3.4163753484221425e-06
**********
}

{
Bridge: PHE292-TYR310
MET: 288
Surface: 288.0
Area: 16.967451859130414
**********
Bridge: TYR310-PHE292
Methionine: 288
Surface: 288 \ CG
Area: 16.967451859097856
**********
}

{
Bridge: PHE47-PHE503
MET: 43
Surface: 43.0
Area: 1.6388880330057552e-06
**********
Bridge: PHE47-PHE503
Methionine: 43
Surface: 43 \ SD
Area: 1.6386662255191633e-06
**********
}

This looks all good.
Next confirm that the right data is being piped into histograms. Confirmed. Good.


It does indeed appear that the closest coordinates to either MT or SD are being chosen
test_code = '2qrb' <- appears correct coordinate being selected as closest to MT
test_code = '4kki' <- appears correct
test_code = '5ef7' <- appears correct
test_code = '3ii2' <- appears correct
"""



    
