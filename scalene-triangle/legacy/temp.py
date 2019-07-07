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

    
    
"""
    
# remove duplicates from data
isolated_data_B = list(set(isolated_data_A))

# split data back out into individual lines using $ char
isolated_data_B = [d.split('$') for d in isolated_data_B]

print(len(isolated_data_A))
print(len(isolated_data_B))

# split operation loads output into lists - flatten
isolated_data_B = list(chain(*isolated_data_B))

# add carriage return to each entry
isolated_data_B = [d + '\n' for d in isolated_data_B]

# send data back to new txt file
with open('results_scalene_cleaned.txt', 'w') as f:
    f.writelines(isolated_data_B)

"""