'''
******************************************************************************
pylift.vmd module
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

import subprocess
import os
from typing import Optional

VMD_EXEC = os.environ.get('VMD_EXEC')

def test_vmd_exec():
    '''
    Tests if VMD executable can be found and run from the environment variable.
    '''
    if not VMD_EXEC:
        print('''
=========================VMD TEST RESULT=========================
VMD_EXEC environment variable not set. 
Please add the location of VMD to your environment variables by running the following in your terminal:
    echo "export VMD_EXEC='PATH/TO/YOUR/VMD_EXEC'" >> ~/.bashrc
    source ~/.bashrc
If you don't know where VMD is, try running:
    whereis vmd
=========================VMD TEST RESULT=========================
        ''')
        return

    try:
        _, stderr, returncode = run_vmd_commands(commands = None, verbose = False)

        if returncode == 0:
            print('=========================VMD TEST RESULT=========================')
            print(f'VMD successfully found at {VMD_EXEC}')
            print('=========================VMD TEST RESULT=========================')
        else:
            print('=========================VMD TEST RESULT=========================')
            print(f"VMD found but exited with an error code: {returncode}\n{stderr}")
            print('=========================VMD TEST RESULT=========================')

    except FileNotFoundError:
        print(f'''
=========================VMD TEST RESULT=========================
FileNotFoundError: VMD executable not found at {VMD_EXEC}.
Please ensure the path is correct and the program exists.
You can add the location of VMD to your environment variables by running the following in your terminal:
    echo "export VMD_EXEC='PATH/TO/YOUR/VMD_EXEC'" >> ~/.bashrc
    source ~/.bashrc
If you don't know where VMD is, try running:
    whereis vmd
=========================VMD TEST RESULT=========================
            ''')
    except Exception as e:
        print(f'An error occurred while trying to execute VMD: {str(e)}')

def run_vmd_commands(commands: str,
                     verbose: Optional[bool] = True) -> tuple[str, str, int]:
    '''
    pylift.vmd.run_vmd_commands
    '''
    try:
        process = subprocess.Popen([VMD_EXEC, '-dispdev', 'text'],
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   text=True)

        stdout, stderr = process.communicate(commands)

        if verbose:
            print("Output:")
            print(stdout)
            print("Errors:")
            print(stderr)

        if process.returncode != 0:
            print(f"VMD exited with an error code: {process.returncode}")

    except FileNotFoundError:
        if verbose:
            print('-----------------------------------------------------------------------')
        print("VMD is not installed or not found in the environment")
        if verbose:
            print('Run test_vmd_exec() from pysimm.apps.vmd for help')
            print('------------------------------------------------------------------------')
    except TypeError:
        if verbose:
            print('------------------------------------------------------------------------')
        print("The location of the VMD program has not been added to the env. variable VMD_EXEC")
        if verbose:
            print('Run test_vmd_exec() from pysimm.apps.vmd for help')
            print('------------------------------------------------------------------------')

    return stdout, stderr, process.returncode

def topo_write(molecule_in: str,
                    lammps_out: str,
                    bonds: Optional[bool] = True,
                    angles: Optional[bool] = True,
                    dihedrals: Optional[bool] = True,
                    impropers: Optional[bool] = True,
                    verbose: Optional[bool] = True) -> None:
    '''
    pysimm.vmd.writelammpsdata
    '''
    cmd = f"""
    mol new {molecule_in}
    package require topotools
    """
    if bonds:
        cmd += "topo retypebonds\n"
    if angles:
        cmd += "topo guessangles\n"
    if dihedrals:
        cmd += "topo guessdihedrals \n "
    if impropers:
        cmd += "topo guessimpropers\n"
    cmd += f"""
    topo writelammpsdata {lammps_out}
    exit
    """
    run_vmd_commands(cmd, verbose = verbose)
    if verbose:
        print(f"""
-----------------------------------------------------------
Used TopoTools to write a skeleton lammps data file...
Input: {molecule_in}
Output: {lammps_out}
The user is cautioned to inspect the file
Ensure proper information supplied in {lammps_out}
-----------------------------------------------------------
            """)
