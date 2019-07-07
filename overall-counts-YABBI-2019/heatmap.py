"""
Written by David S. Weber
dsw7@sfu.ca

Script creates heatmaps for all 2-bridges in MongoDB collections "non_redundant_no_ang_limit" and
"non_redundant_1095_ang_limit"
"""

from numpy import arange
from pandas import DataFrame
from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams['figure.dpi'] = 200

NUM_COLS = 7
NUM_ROWS = 8


# ------------------------------------------------------------------------------------------------------
# starting data from analyze.py

angular_cutoff = True

# bridges from ma.non_redundant_no_ang_limit
data_no_angular_cutoff = {
    '0': {'PHE-PHE': 1952, 'TYR-TYR': 753, 'TRP-TRP': 172, 'TYR-TRP': 774, 'PHE-TYR': 2373, 'PHE-TRP': 1078},
    '1': {'PHE-PHE': 423, 'TYR-TYR': 170, 'TRP-TRP': 56, 'TYR-TRP': 175, 'PHE-TYR': 480, 'PHE-TRP': 267},
    '2': {'PHE-PHE': 560, 'TYR-TYR': 258, 'TRP-TRP': 66, 'TYR-TRP': 177, 'PHE-TYR': 693, 'PHE-TRP': 281},
    '3': {'PHE-PHE': 645, 'TYR-TYR': 332, 'TRP-TRP': 109, 'TYR-TRP': 368, 'PHE-TYR': 833, 'PHE-TRP': 430},
    '4': {'PHE-PHE': 148, 'TYR-TYR': 102, 'TRP-TRP': 8, 'TYR-TRP': 60, 'PHE-TYR': 196, 'PHE-TRP': 78},
    '5': {'PHE-PHE': 149, 'TYR-TYR': 42, 'TRP-TRP': 6, 'TYR-TRP': 26, 'PHE-TYR': 114, 'PHE-TRP': 45},
    '6': {'PHE-PHE': 103, 'TYR-TYR': 50, 'TRP-TRP': 12, 'TYR-TRP': 31, 'PHE-TYR': 157, 'PHE-TRP': 56}
}

# bridges from ma.non_redundant_1095_ang_limit
data_angular_cutoff = {
    '0': {'PHE-PHE': 1938, 'TYR-TYR': 740, 'TRP-TRP': 170, 'TYR-TRP': 753, 'PHE-TYR': 2335, 'PHE-TRP': 1059},
    '1': {'PHE-PHE': 425, 'TYR-TYR': 160, 'TRP-TRP': 54, 'TYR-TRP': 156, 'PHE-TYR': 470, 'PHE-TRP': 261},
    '2': {'PHE-PHE': 568, 'TYR-TYR': 258, 'TRP-TRP': 61, 'TYR-TRP': 175, 'PHE-TYR': 680, 'PHE-TRP': 275},
    '3': {'PHE-PHE': 653, 'TYR-TYR': 323, 'TRP-TRP': 110, 'TYR-TRP': 357, 'PHE-TYR': 840, 'PHE-TRP': 426},
    '4': {'PHE-PHE': 146, 'TYR-TYR': 101, 'TRP-TRP': 8, 'TYR-TRP': 58, 'PHE-TYR': 201, 'PHE-TRP': 78},
    '5': {'PHE-PHE': 151, 'TYR-TYR': 40, 'TRP-TRP': 7, 'TYR-TRP': 25, 'PHE-TYR': 111, 'PHE-TRP': 44},
    '6': {'PHE-PHE': 106, 'TYR-TYR': 48, 'TRP-TRP': 12, 'TYR-TRP': 32, 'PHE-TYR': 161, 'PHE-TRP': 56}
}

if angular_cutoff:
    data = data_angular_cutoff
else:
    data = data_no_angular_cutoff


# ------------------------------------------------------------------------------------------------------
# preprocess data into dataframe

df = DataFrame(data).transpose()
df['Sums'] = df.sum(axis=1)
colsum = DataFrame(df.sum(axis=0)).transpose()
df = df.append(colsum)
ECs = ['N/A'] + ['EC' + str(i) for i in range(1, 7)] + ['Sums']
df.index = ECs


# ------------------------------------------------------------------------------------------------------
# plot

fig, ax = plt.subplots(figsize=(4, 4))
ax.imshow(df, cmap='Pastel1')
ax.set_xticks(arange(len(df.columns)))
ax.set_yticks(arange(len(df.index)))
ax.set_xticklabels(df.columns)
ax.set_yticklabels(df.index)
ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)  # set bridge types to top
plt.setp(ax.get_xticklabels(), rotation=45, ha="left", va="center", rotation_mode="anchor")  # rotate bridge types
# plt.gca().set_frame_on(False)  # remove border

# add grid lines
for line in range(0, NUM_COLS): plt.axvline(line - 0.5, lw=3, c='white')
for line in range(0, NUM_ROWS): plt.axhline(line - 0.5, lw=3, c='white')
plt.axvline(NUM_COLS - 1.5, lw=1, c='k')
plt.axhline(NUM_ROWS - 1.5, lw=1, c='k')

# annotate pixels
for i in range(NUM_COLS):
    for j in range(NUM_ROWS):
        ax.text(i, j, df.iloc[j, i], ha="center", va="center", color="k", size=8)

plt.tight_layout()
plt.show()


