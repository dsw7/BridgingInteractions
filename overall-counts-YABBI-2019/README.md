## Met-aromatic low level routine + bridging interactions
---  
A low level object oriented package for running the Met-aromatic algorithm and collecting statistics on bridging interactions. See DSW thesis for a theoretical description. Met-aromatic data is first collected using ```runner.py``` which loads data into a MongoDB database. The PDB codes of interest were low redundancy (high resolution) structures stored in the file ```low_redundancy_delimiter_list.txt```. The script ```analyze.py``` runs a set of queries on the MongoDB database and returns meaningful data which can be redirected for storage. The script ```heatmap.py``` creates a heatmap representation of bridging interactions broken down by EC classifiers.

### Contents
```
./runner.py                              -- See above
./analyze.py                             -- See above
./heatmap.py                             -- See above
./low_redundancy_delimiter_list.txt      -- Delimiter list of PDB codes used in bridging interaction study
./utils/filegetter.py                    -- Fetches PDB files over ftp
./utils/ma.py                            -- Contains Met-aromatic class
./utils/utils.py                         -- Contains Met-aromatic helper functions
./utils/apply_angular_limit_to_no_ang.py -- Contains a method of applying angular limit to an existing MongoDB collection
./tests/utils_init/                      -- Contains some of the first ever Met-aromatic implementations
./tests/randomized_pdb_codes.csv         -- A .csv containing random PDB test codes
./tests/test.py                          -- Unit tests executed here
./figures/no_angular_cutoff.png          -- Figure obtained from heatmap.png - no angular cutoff applied to starting data
./figures/1095_angular_cutoff.png        -- Figure obtained from heatmap.png - 109.5 degree cutoff applied to starting data
```

### Usage: runner.py
---
To run the program:
```
$ python runner.py <args>
```
The program requires, at bare minimum, either a valid PDB code or a path to a text file containing a list of PDB codes to analyze. Here the arbitrary PDB code 1rcy is analyzed:
```
$ python runner.py --code 1rcy
```
A batch job can be performed as follows:
```
$ python runner.py --batch /path/to/low_redundancy_delimiter_list.txt
```
The PDB codes in the batch job text file should be separated by newline characters. *NOTE:* Both code and batch parameters cannot be passed simultaneously. Next come the Met-aromatic algorithm constraints:
```
$ python runner.py --code 1rcy --cutoff 4.9 --angle 90.0 --model cp
```
Here the cutoff has been set to 4.9 Angstroms (the max norm of vector *v*) and the maximum angle of either Met-theta or Met-phi cannot exceed 90.0 degrees. The model used to interpolate lone pair positions is cp or Cross Product. These parameters do not have to be passed. Default values are used if these values are not specified. Defaults can be obtained by reading:
```
$ python runner.py --help
```
Console output is normally suppressed. Suppression can be lifted by passing the verbose parameter:
```
$ python runner.py --code 1rcy --verbose
```
There are several options available for working with output data. Data can either be exported to a .csv file or a MongoDB database. Export to a MongoDB database is recommended. Data cannot be exported to both a .csv file and a MongoDB database simultaneously. To save to a .csv file:
```
$ python runner.py --code 1rcy --export-csv /path/to/output.csv
```
Or a MongoDB database:
```
$ python runner.py --code 1rcy --export-mongo
```
MongoDB export can be modified as follows:
```
$ python runner.py --code 1rcy --export-mongo --mongoport 27017 --mongohost localhost --database my_database --collection my_collection
```
Default MongoDB parameters are passed if no export parameters are specified. No data is saved if no export parameter is passed. As always, defaults can be obtained using:
```
$ python runner.py --help
```

### Usage: analyze.py
---
To run the script:
```
$ python analyze.py
```
Script interpretation will be terminated if data was not previously collected using ```runner.py```. Otherwise, data will be printed to the console. Output can be redirected to a file for storage:
```
$ python analyze.py > /path/to/results/results.txt
```

### Usage: heatmap.py
---
There's really little to this script. I literally copy pasted data from ```analyze.py``` into this script. To run:
```
$ python heatmap.py
```
Which yields (no angular cutoff applied):
<p align="center">
  <img width="420" height="400" src="https://github.com/dsw7/BridgingInteractions/tree/master/overall-counts-YABBI-2019/figures/no_angular_cutoff.png">
</p>
And (109.5 degree angular cutoff applied):
<p align="center">
  <img width="420" height="400" src="https://github.com/dsw7/BridgingInteractions/tree/master/overall-counts-YABBI-2019/figures/1095_angular_cutoff.png">
</p>

