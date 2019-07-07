# ignore this file
# for dsw use only

import os
import pandas as pd

root_path = r'/Volumes/MSC/Databases/2018 Bridging Databases'

# paths to bridging_n_n+5000.csv datasets - these have the metals
paths = os.listdir(root_path)
paths = [item for item in paths if item[0:8] == 'bridging' and os.path.splitext(item)[1] == '.csv']

# paths to .csvs in bridging_proc
paths_to_bridging = os.path.join(root_path, 'bridging_proc')
paths_to_bridging = os.listdir(paths_to_bridging)

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
        
pdb_codes = []        
for idx in range(0, df_METALS.shape[0]):
    # create a temporary df in a fashion similar to Pandas groupby machinery
    df = df_BRIDGES[df_BRIDGES['0'] == df_METALS.iloc[idx, 1]]
    PDB_CODE = df_METALS.iloc[idx, 1][0:4]
    if df.empty == False: 
        pdb_codes.append(PDB_CODE + '\n')
    
pdb_codes = list(set(pdb_codes))  # length should be 8,743

with open('delim.txt', 'w') as f: f.writelines(pdb_codes)