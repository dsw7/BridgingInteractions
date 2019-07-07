"""
dsw7@sfu.ca
This namespace houses a function for getting metals from PDB files.
"""

from re import search
pat = 'SC|TI|V|CR|MN|FE|CO|NI|CU|ZN'
# pat = 'MN|FE|CO|NI|CU' # isolate to metals that can actually undergo redox
                         # jeff didn't want this


def get_metals(path_to_pdb_file, chain='A'):
    """ Function extracts metal coordinates from PDB file. """

    # open the PDB file which should be in pwd
    with open(path_to_pdb_file) as f: data = f.readlines()

    # split lines by whitespace
    whitesplit = [lines.split() for lines in data]

    # get only first model in multimodel entries
    firstmodel = []
    for line in whitesplit:
        firstmodel.append(line)
        if firstmodel[0] == 'ENDMDL': break

    # get only HETATM records
    hetatm_records = [line for line in firstmodel if line[0] == 'HETATM']

    # get only a single chain
    firstchain = [line for line in hetatm_records if line[4] == chain]

    # get only metal containing lines
    return [line for line in firstchain if search(pat, line[2])]


    