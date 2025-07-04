"""
pylift.utilities module

License: The MIT License (MIT)

Copyright (c) 2024 Brandon C. Tapia
"""

import fnmatch
import os
import json
# import importlib
import importlib.resources as pkg_resources

# optional_modules = ["rdkit", "openbabel", "pymol"]
imported_modules = {}

# for module in optional_modules:
#    try:
#        imported_modules[module] = importlib.import_module(module)
#    except ImportError:
#        imported_modules[module] = None


def cleanup_pylift(
    user_files: list[str] = None,
    temp: bool = True,
    antechamber: bool = True,
    verbose: bool = True,
) -> dict:
    """
    pylift.utilities.cleanup_pylift

    Cleans up temporary files generated during program execution.

    Arguments:
        user_files (list[str]): List of additional file paths to be removed. Default is None.
        temp (bool) : Removes all files with tmp in their name
        antechamber (bool) : Removes all extraneous files generated with Antechamber

    Returns:
        dict: information regarding the success/failure of deletion
    """
    default_files = []

    file_list = default_files.copy()
    if user_files:
        file_list.extend(user_files)

    for root, _, files in os.walk("."):
        if temp:
            for filename in fnmatch.filter(files, "*tmp*"):
                file_list.append(os.path.join(root, filename))
        if antechamber:
            for filename in fnmatch.filter(files, "ANTECHAMBER*"):
                file_list.append(os.path.join(root, filename))
            for filename in ["ATOMTYPE.INF", "sqm.in", "sqm.out", "sqm.pdb"]:
                file_list.append(os.path.join(root, filename))

    result = {}
    for file_path in file_list:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
                result[file_path] = "Removed"
            else:
                result[file_path] = "File does not exist"
        except Exception as e:
            result[file_path] = f"Failed to remove: {e}"

    if verbose:
        print("[cleanup_pylift] cleanup complete")

    return result


def read_json(in_json):
    """
    pylift.utilities.read_json

    Reads a json file from pylift/ff_data directory
    """
    with pkg_resources.path("pylift.ff_data", in_json) as file_path:
        with open(file_path, encoding="utf-8") as file:
            dict_file = json.load(file)

            print(f"[read_json] read {in_json}")
            return dict_file


def write_json(dict_loc: dict, out_json: str) -> None:
    """ 
    *
    """

    with open(dict_loc, "w", encoding="utf-8") as file:
        json.dump(out_json, file, indent=4)
