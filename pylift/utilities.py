'''
******************************************************************************
pylift.utilities module
*******************************************************************************

*******************************************************************************
License
*******************************************************************************
The MIT License (MIT)

Copyright (c) 2024 Brandon C. Tapia

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

import fnmatch
import os
import json

def cleanup_pylift(user_files: list[str] = None,
                       temp: bool = True,
                       antechamber: bool = True) -> dict:
    """
    PyLIFT.utilities.cleanup_pylift

    Cleans up temporary files generated during program execution
    
    Parameters:
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

    for root, _, files in os.walk('.'):
        if temp:
            for filename in fnmatch.filter(files, '*tmp*'):
                file_list.append(os.path.join(root, filename))
        if antechamber:
            for filename in fnmatch.filter(files, 'ANTECHAMBER*'):
                file_list.append(os.path.join(root, filename))
            file_list.append(os.path.join(root, 'ATOMTYPE.INF'))

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
    print('[md_builder] cleanup complete')

    return result

def read_json(in_json: str) -> dict:
    '''
    '''
    with open(in_json, 'r', encoding='utf-8') as file:
        dict_file = json.load(file)

    return dict_file

def write_json(dict_loc: dict,
               out_json: str) -> None:
    '''
    '''

    with open(dict_loc, 'w', encoding='utf-8') as file:
        json.dump(out_json, file, indent=4)
