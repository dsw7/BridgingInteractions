"""
dsw7@sfu.ca
The "main" script that accepts PDB code input from user.
"""

import os
from sys import path; path.append("utils")
from ma import MetAromatic
from pprint import pprint
from argparse import ArgumentParser, RawTextHelpFormatter
from pymongo import MongoClient
from hashlib import md5
from time import time

COLUMNS = ["ARO", "ARO POS", "MET", "MET POS", "NORM", "MET-THETA", "MET-PHI"]
DEFAULT_PORT = 27017
DEFAULT_HOST = "localhost"
DB = "ma"
COL = "ma"

msg_code = 'Process a pdb code. \nUsage: $ python runner.py --code <1abc>'
msg_batch = 'Process a batch of pdb codes. \nUsage: $ python runner.py --batch /path/to/foo.txt'
msg_cutoff = 'Set a Euclidean cutoff. \nDefault = 6.0 Angstroms. \nUsage: $ python runner.py --cutoff <float>'
msg_angle = 'Set Met-theta/Met-phi angle. \nDefault = 109.5 degrees. \nUsage: $ python runner.py --angle <float>'
msg_model = 'Set a lone pair interpolation model. \nDefault = cp. \nUsage: $ python runner.py --model <cp|rm>'
msg_verbosity = 'Set output verbosity. \nUsage: $ python runner.py --verbose'
msg_export_csv = 'Export results to csv. \nUsage: $ python runner.py --export-csv /path/to/bar.txt'
msg_export_mongo = 'Export results to MongoDB. \nUsage: $ python runner.py --export-mongo'
msg_port = 'Set a MongoDB port. \nDefault = 27017. \nUsage: $ python runner.py --mongoport <port>'
msg_host = 'Set a MongoDB host. \nDefault = localhost. \nUsage: $ python runner.py --mongohost <host>'
msg_db = 'Choose a MongoDB export database name. \nDefault = ma. \nUsage: $ python runner.py --database <name>'
msg_col = 'Choose a MongoDB export collection name. \nDefault = ma. \nUsage: $ python runner.py --collection <name>'

parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('--code', help=msg_code, default='0', type=str)
parser.add_argument('--batch', help=msg_batch, default='0', type=str)
parser.add_argument('--cutoff', help=msg_cutoff, default=6.0, type=float)
parser.add_argument('--angle', help=msg_angle, default=109.5, type=float)
parser.add_argument('--model', help=msg_model, default='cp')
parser.add_argument('--verbose', help=msg_verbosity, action='store_true')
parser.add_argument('--export-csv', help=msg_export_csv, default='False', dest='export_csv')
parser.add_argument('--export-mongo', help=msg_export_mongo, action='store_true', dest='export_mongo')
parser.add_argument('--mongoport', help=msg_port, default=DEFAULT_PORT, type=int)
parser.add_argument('--mongohost', help=msg_host, default=DEFAULT_HOST, type=str)
parser.add_argument('--database', help=msg_db, default=DB, type=str)
parser.add_argument('--collection', help=msg_col, default=COL, type=str)

code = parser.parse_args().code
path = parser.parse_args().batch
cutoff = parser.parse_args().cutoff
angle = parser.parse_args().angle
model = parser.parse_args().model
verbose = parser.parse_args().verbose
export_csv = parser.parse_args().export_csv
export_mongo = parser.parse_args().export_mongo
mongoport = parser.parse_args().mongoport
mongohost = parser.parse_args().mongohost
database = parser.parse_args().database
collection = parser.parse_args().collection


def verify_user_input():
    if (code == '0') and (path == '0'):
        exit('No PDB code or path to .txt given.')
    elif (len(code) != 4) and (code != '0'):
        exit('Invalid pdb code: {}'.format(code))
    elif (code != '0') and (path != '0'):
        exit('Cannot choose between .txt file and pdb code.')
    elif model not in ('cp', 'rm'):
        exit("Invalid model. Valid models are: cp (Cross Product) or rm (Rodrigues' method).")
    elif (angle < 0.0) or (angle > 360.00):
        exit('Angle must be between 0 and 360 degrees.')
    elif cutoff < 0:
        exit('Cutoff must be greater than or equal to 0.0 Angstroms.')
    elif export_mongo and (export_csv != 'False'):
        exit('Cannot export to both MongoDB and a .csv document simultaneously.')
    else:
        pass


def print_args():
    print("Analyzing: {}".format(code))
    print("Cutoff: {}".format(cutoff))
    print("Angle: {}".format(angle))
    print("Model: {}".format(model))
    print("Mongo Port: {}".format(mongoport))
    print("Mongo Host: {}".format(mongohost))
    print("Database Name: {}".format(database))
    print("Collection Name: {}".format(collection))
    print("Export to MongoDB: {}".format(export_mongo))
    print("Export to csv: {}\n".format(export_csv))


def read_pdb_code_txt_file(filepath):
    with open(filepath, 'r') as f:
        codes = f.read().splitlines()

    if not codes:
        exit('Empty .txt file.')
    else:
        return codes


def run_met_aromatic(pdbcode):
    try:
        ma = MetAromatic(pdbcode, cutoff=cutoff, angle=angle, model=model)
        return ma.met_aromatic(), ma.get_ec_classifier()
    except Exception as exception:
        print('An exception has occurred:')
        print(exception)


def mapper(result, pdbcode):
    # a function for adapting Met-aromatic results to MongoDB
    outgoing, ec = [], result[1]
    for item in result[0]:
        doc = {
            "code": pdbcode,
            "aro": item[0],
            "arores": item[1],
            "met": item[3],
            "norm": item[4],
            "met-theta": item[5],
            "met-phi": item[6],
            "ec": ec
        }

        # overwrite MongoDB _id with custom _id to prevent writing duplicate data into database
        _id = ''.join([str(i) for i in doc.values()])
        doc['_id'] = md5(_id.encode()).hexdigest()
        outgoing.append(doc)
    return outgoing


def print_progress(pdbcode, current_count, overall_count):
    percent = round(current_count * 100 / overall_count, 2)
    msg = '{}. Iteration {} out of {}. {} % complete.'.format(pdbcode, current_count, overall_count, percent)
    print(msg)


if __name__ == '__main__':
    verify_user_input()
    print_args()

    if export_mongo:
        client = MongoClient(mongohost, mongoport)
        db = client[database]
        col = db[collection]

    if (code != '0') and (path == '0'):  # user inputs a valid pdb code but no path to batch file
        results = run_met_aromatic(code)

        if results is None:
            print('NoneType object was returned from MetAromatic algorithm.')
        elif not results[0]:
            print('No interactions.')
        else:
            if verbose:
                pprint(mapper(results, code))

            if export_mongo:
                col.insert_many(mapper(results, code))
            # TODO: else export to csv... might remove this -> .csvs are really not a good way to work with data

    elif (code == '0') and (path != '0'):  # user inputs no pdb code but valid path to batch file
        if not os.path.exists(path):
            exit('Path to file does not exist.')
        else:
            pdb_codes = read_pdb_code_txt_file(path)
            overall = len(pdb_codes)
            start = time()

        for u, code in enumerate(pdb_codes, 1):
            print_progress(code, u, overall)
            results = run_met_aromatic(code)

            if results is None:
                print('NoneType object was returned from MetAromatic algorithm.')
            elif not results[0]:
                continue
            else:
                if verbose:
                    pprint(mapper(results, code))

                if export_mongo:
                    col.insert_many(mapper(results, code))
                # TODO: else export to csv...

        print('\n' + '-' * 50)
        print('Total processing time: {} s'.format(time() - start))
