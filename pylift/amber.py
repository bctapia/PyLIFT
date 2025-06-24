"""
pylift.amber module

License: The MIT License (MIT)

Copyright (c) 2024 Brandon C. Tapia
"""

import subprocess
import os
from typing import Optional
from pylift import utilities

ANTECHAMBER_EXEC = os.environ.get("ANTECHAMBER_EXEC")


def test_antechamber_exec():
    """
    Tests if Antechamber executable can be found and run from the environment variable.

    Returns:
        bool for success/fail
    """

    if not ANTECHAMBER_EXEC:
        print(
            """
====================ANTECHAMBER TEST RESULT====================
ANTECHAMBER_EXEC environment variable not set. 
Please add the location of Antechamber to your environment variables by running the following in your terminal:
    echo "export ANTECHAMBER_EXEC=PATH/TO/YOUR/ANTECHAMBER_EXEC" >> ~/.bashrc
    source ~/.bashrc
If you don't know where Antechamebr is, try running:
    whereis antechamber
====================ANTECHAMBER TEST RESULT====================
        """
        )
        return False

    if not os.path.isfile(ANTECHAMBER_EXEC):
        print(
            """
====================ANTECHAMBER TEST RESULT====================
ANTECHAMBER_EXEC environment variariable does not point to a file. 
Please add the location of Antechamber to your environment variables by running the following in your terminal:
    echo "export VMD_EXEC='PATH/TO/YOUR/ANTECHAMBER_EXEC'" >> ~/.bashrc
    source ~/.bashrc
If you don't know where Antechamebr is, try running:
    whereis antechamber
====================ANTECHAMBER TEST RESULT====================
        """
        )
        return False

    try:
        result = subprocess.run(
            [ANTECHAMBER_EXEC, "-L"], capture_output=True, text=True, check=True
        )
        if result.returncode == 0:
            print("====================ANTECHAMBER TEST RESULT====================")
            print(f"Antechamber successfully found at {ANTECHAMBER_EXEC}")
            print("====================ANTECHAMBER TEST RESULT====================")
            return True
        else:
            print("====================ANTECHAMBER TEST RESULT====================")
            print(
                f"Antechamber exited with error code: {result.returncode}\n{result.stderr}"
            )
            print("====================ANTECHAMBER TEST RESULT====================")
            return False
    except subprocess.CalledProcessError as e:
        print(f"Antechamber executable failed with error: {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False


def antechamber(
    mol2_in: str,
    mol2_out: str,
    forcefield: str,
    charge_method: Optional[str] = None,
    missing_search: Optional[str] = "parmchk2",
    verbose: Optional[bool] = True,
):
    """
    Calls Antechamber from AmberTools to apply a forcefield and find missing parameters

    Arguments:
        mol2_in (str): mol2 input file
        mol2_out (str): mol2 output file with forcefield applied
        forcefield: which forcefield to use
        charge_method (str): charge method to use, default is None
            Recommended to use QM calculations rather than rely on Antechamber for charges
        missing_search (str): program used to find missing parameters, default is 'parmchk2'
            Old versions of AmberTools might require 'parmchk' instead of 'parmchk2'
        verbose (bool): prints stderr and stdout from Antechamber, default is verbose
    """
    command = f"{ANTECHAMBER_EXEC} -i {mol2_in} -fi mol2 -o {mol2_out} -fo mol2 -at {forcefield}"

    if charge_method:
        command += f" -c {charge_method}"

    utilities.cleanup_pylift(temp=False, verbose=False)

    result = subprocess.run(command.split(), capture_output=True, text=True, check=True)

    if verbose:
        if result.stdout:
            print("Antechamber stdout:")
            print(result.stdout)
        if result.stderr:
            print("Antechamber stderr:")
            print(result.stderr)

    if missing_search:
        if os.path.exists("missing_ff_params.frcmod"):
            os.remove("missing_ff_params.frcmod")

        command = f"{missing_search} -i {mol2_out} -f mol2 -o missing_ff_params.frcmod -s {forcefield}"
        result = subprocess.run(
            command.split(), capture_output=True, text=True, check=True
        )

        if verbose:
            if result.stdout:
                print(f"{missing_search} stdout:")
                print(result.stdout)
            if result.stderr:
                print(f"{missing_search} stderr:")
                print(result.stderr)

    print("[antechamber] completed")
