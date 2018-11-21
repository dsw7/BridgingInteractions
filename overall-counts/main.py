"""
Written by David S. Weber

The main .py script for the bridging interaction study.
Basic overview:
    
(1) Script downloads a JSON file of all PDB codes from PDB
(2) Script converts JSON to list
(3) Script sets up a SQL database for storing data
(4) Script iterates for entries in list
(5) Script downloads the PDB file corresponding to each iteration
(6) Script applies Met-Aromatic algorithm to each local PDB file
(7) Script sends Met-Aromatic output to the SQL database
"""


# ------------------------------------------------------------------------------
# input parameters to be specified by researcher

# search parameters
# --------------
START = 0
END   = 10  # set to -1 to iterate over entire PDB


# PDB parameters
# --------------
CHAIN = 'A'       # input a chain of interest here
CUTOFF = 6.0      # input some cutoff vector v norm
ANGLE = 109.5     # input some cutoff vector a / vector v and vector g / vector v angle
MODEL = 'cp'      # 'cp' or 'rm' -> Cross Product or Rodrigues' Method lone pair interpolation methods


# ------------------------------------------------------------------------------

from sys                    import path as PATH; PATH.append('libs')
from requests               import get
from json                   import loads
from PDB_filegetter         import PDBFile
from ma_lowlevel            import met_aromatic
from platform               import node
from time                   import sleep


# ------------------------------------------------------------------------------
# prepare sql database

import pyodbc

DRIVER   = 'SQL Server'
SERVER   = '{}\SQLEXPRESS'.format(node())
DATABASE = 'master'
TRUSTED  = 'yes'

str_conn = """
           Driver={};
           Server={};
           Database={};
           Trusted_Connection={}
           """.format(DRIVER, SERVER, DATABASE, TRUSTED)

connection = pyodbc.connect(str_conn, autocommit=True)
cursor = connection.cursor()

try:  # create new database if database doesn't exist
    sqlcommand = " CREATE DATABASE Bridging "
    cursor.execute(sqlcommand)
    sqlcommand = """
    CREATE TABLE Bridging.dbo.BridgingResults (
        PDBCODE varchar(4),
        ARO varchar(3),
        AROResPos int,
        MET varchar(3),
        METResPos int,
        Norm decimal(18, 3),
        MetTheta decimal(18, 3),
        MetPhi decimal(18, 3)
    ); """ 
    cursor.execute(sqlcommand)
except Exception as exception:
    status = """ Database requirement satisfied. Here is the error code: {} """.format(exception)
    print(status)

# close existing connection and reconnect to the new database
# fix for SSMS error 42S22
connection.close()
sleep(0.05)
DATABASE = 'Bridging'
str_conn = """
           Driver={};
           Server={};
           Database={};
           Trusted_Connection={}
           """.format(DRIVER, SERVER, DATABASE, TRUSTED)

connection = pyodbc.connect(str_conn)
cursor = connection.cursor()


# ------------------------------------------------------------------------------
# core loop


# get current list of all PDB codes from RCSB PDB
# -----------------------------------------------
url = 'https://www.rcsb.org/pdb/json/getCurrent'
current_PDB_files_JSON       = get(url)
current_PDB_files_pyDict     = loads(current_PDB_files_JSON.content)
current_PDB_files_pyList     = current_PDB_files_pyDict.get('idList')
current_PDB_files_entrycount = current_PDB_files_pyDict.get('resultCount')
if END == -1: END = current_PDB_files_entrycount


# apply met aromatic to each entry in JSON
# send results to sql database
# ----------------------------
for iteration, CODE in enumerate(current_PDB_files_pyList[START:END]):
    status = '{} - Analyzed {} of {} PDB entries.'.format(CODE, iteration + 1, current_PDB_files_entrycount)
    print(status)
    
    file_pdb = PDBFile(CODE)
    path_to_file = file_pdb.fetch_from_PDB()
    
    if path_to_file == 'URLError':
        print('URLError')
        continue
    else:
        try:  # catch other exceptions such as "string crashes", missing coordinates, etc., in PDB files
            data_retrieved = met_aromatic(CHAIN=CHAIN, 
                                          CUTOFF=CUTOFF, 
                                          ANGLE=ANGLE, 
                                          MODEL=MODEL, 
                                          filepath=path_to_file)
        except Exception as exception:
            print(exception)
        file_pdb.clear()
    
    if len(data_retrieved) == 0:  # next iteration if no interactions
        continue
    else:                         # write to sql
        for line in data_retrieved:
            vals = (("'{}', " * 7) + "'{}'").format(CODE, *line)
            outgoing = """ INSERT INTO BridgingResults (
                           PDBCODE, ARO, AROResPos, MET, METResPos, Norm, MetTheta, MetPhi) 
                           VALUES ({}) 
                       """.format(vals)
            cursor.execute(outgoing)
            connection.commit()
            
connection.close()
