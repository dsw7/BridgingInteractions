## Script for collecting bridging statistics for all proteins in PDB  
Written by David S. Weber

---  
## Instructions  

### Step 1
Ensure both Python3 and SQL* are installed.  

### Step 2  
Edit the main.py script. Assign a range of PDB files to iterate over:

    # find START and END variables and assign the following:
    START = 0
    END = -1
    # to iterate over entire PDB.
    # Another example: START = 1000 and END = 2000 will iterate over 1000 files.

### Step 3  
Set Met-Aromatic algorithm parameters. The default settings are the recommended settings.

    CHAIN = 'A'       <- i.e. PDB file A delimited chains
    CUTOFF = 6.0      <- vector v magnitude of 6.0 Angstroms
    ANGLE = 109.5     <- Met-theta or Met-phi angles <= 109.5 degrees
    MODEL = 'cp'      <- Cross product or Rodrigues' method models. Contact dsw7 for more details

### Step 4
Execute the script. Mining the entire PDB should take about a day (dependent on network bit rate).

### Step 5
Execute counts.py to analyze data in the SQL* database once the mining job is complete.

### Step 6
The counts.py assumes that the NetworkX library (https://networkx.github.io/) is available. 
NetworkX can be be installed with pip.  

A user fluent in the Python programming language can modify the script to dump output
locally should the research group be uninterested in working with SQL databases.*

---  
## Flow scheme
                                                --------------
    ~\overall-counts\libs ==>  ( main.py ) --> | SQL database | --> ( counts.py ) --> HEATMAP
                                                --------------                                           
---  
## Known issues

<p align="justify">
Iterating over pseudorandomly chosen PDB files has revealed that some files throw
exceptions.   
</p>

Examples:

    1A7S -> Missing methionine SD data (?) will throw an exception in ma_lowlevel.py   
    1EK1 -> "String crash" - Occasional exception that results from being unable to tease apart 
             coordinates of form '-44.860-108.842' into two separate floats.

<p align="justify">
The script may stall if the network connection is unstable. The script may be
restarted at the stall position by updating the START variable with the integer
representing the location of the stall.
</p>

Example:  

    # input
    START = 0  
    END = -1  

    1ABC - Analyzed 516 of 146502 PDB entries.
    1BBC - Analyzed 517 of 146502 PDB entries.
    1CBC - Analyzed 518 of 146502 PDB entries.
    ... stall ...

Kill the prompt and reset START = 518. Re-execute the script.  

    # input
    START = 518  
    END = -1  
