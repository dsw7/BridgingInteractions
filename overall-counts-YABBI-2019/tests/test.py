"""
dsw7@sfu.ca
Unit tests for new refactored program. I am comparing against the program I wrote
in grad school.

Run with command:
    $ python -m pytest -v -s test.py

NOTES:
    [1] Dropbox syncing will lag system and leave files in pwd. Disable if possible
    or move the package out of Dropbox.
    [2] String crashes leave files in pwd (the standard does this - not the new code).

"""

import pytest
from sys import path; path.append(r"../utils")
from ma import MetAromatic; path.append(r"utils_init")
from ma_lowlevel import met_aromatic
from PDB_filegetter import PDBFile
from time import sleep
from pandas import DataFrame, read_csv, testing

COLUMNS = ["ARO", "ARO POS", "MET", "MET POS", "NORM", "MET-THETA", "MET-PHI"]
CHAIN = "A"
CUTOFF = 6.0
ANGLE = 109.5
MODEL = "cp"
START = 3000
END = 4000
ERRORS_TO_IGNORE = (ValueError, PermissionError)

df = read_csv("randomized_pdb_codes.csv")
df = df.iloc[START:END]
df.columns = ["CODE"]


def check_results_equality(result_a, result_b, sort_by="NORM", col=COLUMNS):
    if (result_a == []) or (result_b == []):
        assert result_a == result_b, "One result is non-empty!"
        print("No aromatic interactions.")
    else:
        df_result_a = DataFrame(result_a)
        df_result_a.columns = col
        df_result_a = df_result_a.sort_values(by=[sort_by]).reset_index(drop=True)

        df_result_b = DataFrame(result_b)
        df_result_b.columns = col
        df_result_b = df_result_b.sort_values(by=[sort_by]).reset_index(drop=True)

        testing.assert_frame_equal(df_result_a, df_result_b)


@pytest.mark.parametrize("code", df.CODE.tolist())
def test_random_pdb_codes(code):
    try:
        # reference program
        file_pdb = PDBFile(code)
        path_to_file = file_pdb.fetch_from_PDB()
        standard = met_aromatic(CHAIN=CHAIN, CUTOFF=CUTOFF, ANGLE=ANGLE, MODEL=MODEL, filepath=path_to_file)
        file_pdb.clear()

        # refactored program
        sleep(0.50)  # sleep to deal with dropbox lag if enabled
        test = MetAromatic(code, chain=CHAIN, cutoff=CUTOFF, angle=ANGLE, model=MODEL).met_aromatic()

        check_results_equality(standard, test)
    except ERRORS_TO_IGNORE as error:
        print(error)
        print("Try temporarily disabling Dropbox in the event of PermissionError.")
        pass

