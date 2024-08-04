'''
pylift.reader module

License: The MIT License (MIT)

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
import os
import fnmatch
import re
import importlib
import json
import copy
from typing import Optional

def read_mol2(input_file: str,
              out_json: Optional[str] = None) -> dict:
    '''
    PyLIFT.reader.read_mol2

    Reads a MOL2 file and extracts the molecule, atom, bond, and substructure sections
    into separate dictionaries, which are then combined into one dictionary.
    '''
    with open(input_file, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    molecule_section = False
    atom_section = False
    bond_section = False
    substructure_section = False
    
    remaining_atoms = False
    remaining_end = False

    atom_dict = {}
    bond_dict = {}
    auxinfo_dict = {}

    molecule_iter = 0

    res = None
    atom_num = None
    bond_num = None
    ident_1 = None
    ident_2 = None
    mol_type = None
    charge_method = None
    substructure = None

    for i, line in enumerate(lines):
        if not line:
            continue

        line = line.strip()
        columns = line.split()

        if line.startswith('@<TRIPOS>MOLECULE'):
            molecule_section = True
            atom_section = False
            bond_section = False
            substructure_section = False
            molecule_init = line

        elif line.startswith('@<TRIPOS>ATOM'):
            molecule_section = False
            atom_section = True
            bond_section = False
            substructure_section = False
            atom_init = line

        elif line.startswith('@<TRIPOS>BOND'):
            molecule_section = False
            atom_section = False
            bond_section = True
            substructure_section = False
            bond_init = line

        elif line.startswith('@<TRIPOS>SUBSTRUCTURE'):
            molecule_section = False
            atom_section = False
            bond_section = False
            substructure_section = True
            substructure_init = line

        elif molecule_section:
            if molecule_iter == 0:
                res = line
            if len(columns) > 0:
                if columns[0].isdigit() and columns[1].isdigit():
                    atom_num = int(columns[0])
                    bond_num = int(columns[1])
                    ident_1 = int(columns[3])
                    ident_2 = int(columns[4])
            if molecule_iter == 2:
                mol_type = line
            if molecule_iter == 3:
                charge_method = line
            molecule_iter += 1

        elif atom_section:
            atom_dict[columns[0]] = {
                'atom_name': columns[1],
                'x': float(columns[2]),
                'y': float(columns[3]),
                'z': float(columns[4]),
                'atom_type': columns[5],
                'molecule_num': int(columns[6]),
                'res': columns[7],
                'charge': float(columns[8]),
            }

        elif bond_section:
            bond_dict[columns[0]] = {
                'atom_1': int(columns[1]),
                'atom_2': int(columns[2]),
                'bond_type': columns[3]
            }

        elif substructure_section:
            substructure = line

    auxinfo_dict['header_info'] = {
        'molecule_init': molecule_init,
        'atom_init': atom_init,
        'bond_init': bond_init,
        'substructure_init': substructure_init
    }

    auxinfo_dict['molecule_info'] = {
        'res': res,
        'atom_num': atom_num,
        'bond_num': bond_num,
        'ident_1': ident_1,
        'ident_2': ident_2,
        'mol_type': mol_type,
        'charge_method': charge_method
    }

    auxinfo_dict['substructure_info'] = {
        'substructure': substructure
    }

    mol2_dict = {
        'auxinfo_dict': auxinfo_dict,
        'atom_dict': atom_dict,
        'bond_dict': bond_dict
    }

    overall_charge = sum(float(atom_info['charge']) for atom_info in atom_dict.values())
    mol2_dict['auxinfo_dict']['molecule_info'].update({'total_charge': overall_charge})

    if out_json:
        with open(out_json, 'w', encoding='utf-8') as file:
            json.dump(mol2_dict, file, indent=4)
        print(f'printed .json file ({out_json}) containing information from {input_file}')

    print(f'[read_mol2] read {input_file}')
    return mol2_dict

def read_topo(pseudoatoms: str,
            input_file: str = 'topo.tmp.lmps',
            linker_identifier: Optional[str] = 'L',
            out_json: Optional[str] = None) -> dict:
    '''
    '''
    if pseudoatoms is None:
        print('Variable psuedoatoms must be set as True or False')
        return

    with open(input_file, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    pair_section = False
    bond_section = False
    angle_section = False
    dihedral_section = False
    improper_section = False
    mass_section = False
    atom_section = False
    
    remaining = False
    atom_remaining = False

    header_section = True

    mass_dict = {}
    pair_dict = {}
    bond_dict = {}
    angle_dict = {}
    atom_dict = {}
    dihedral_dict = {}
    improper_dict = {}

    header = {}
    footer = {}

    for i, line in enumerate(lines):

        line = line.strip()
        pattern = r'[ \-]+'
        columns = re.split(pattern, line)

        if not line:
            continue

        if columns[0] == 'Masses':
            pair_section = False
            bond_section = False
            angle_section = False
            dihedral_section = False
            improper_section = False
            header_section = False
            mass_section = True
            atom_section = False
        
        elif columns[0] == 'Atoms':
            pair_section = False
            bond_section = False
            angle_section = False
            dihedral_section = False
            improper_section = False
            header_section = False
            mass_section = False
            atom_section = True
        
        elif columns[0] == 'Bonds':
            pair_section = False
            bond_section = False
            angle_section = False
            dihedral_section = False
            improper_section = False
            header_section = False
            mass_section = False
            atom_section = False
        
        elif mass_section:
            mass_dict[columns[0]] = {
                    'atom': columns[3],
                    'mass': columns[1]
                }
        
        elif atom_section:
            atom_dict[columns[0]] = {'molecule': int(columns[1]),
                                     'atom_type': int(columns[2]),
                                     'charge': float(columns[3]),
                                     'x': float(columns[4]),
                                     'y': float(columns[5]),
                                     'z': float(columns[6]),
                                     'comment': columns[8]
            }

        elif len(columns) > 1:
            if columns[1] == 'Pair':
                pair_section = True
                bond_section = False
                angle_section = False
                dihedral_section = False
                improper_section = False

            elif columns[1] == 'Bond':
                pair_section = False
                bond_section = True
                angle_section = False
                dihedral_section = False
                improper_section = False

            elif columns[1] == 'Angle':
                pair_section = False
                bond_section = False
                angle_section = True
                dihedral_section = False
                improper_section = False

            elif columns[1] == 'Dihedral':
                pair_section = False
                bond_section = False
                angle_section = False
                dihedral_section = True
                improper_section = False

            elif columns[1] == 'Improper':
                pair_section = False
                bond_section = False
                angle_section = False
                dihedral_section = False
                improper_section = True

            elif pair_section:
                pair_dict[columns[1]] = {
                    'atom': columns[2]
                }

            elif bond_section:
                bond_dict[columns[1]] = {
                    'atom_1': columns[2],
                    'atom_2': columns[3]
                }

            elif angle_section:
                angle_dict[columns[1]] = {
                    'atom_1': columns[2],
                    'atom_2': columns[3],
                    'atom_3': columns[4]
                }

            elif dihedral_section:
                dihedral_dict[columns[1]] = {
                    'atom_1': columns[2],
                    'atom_2': columns[3],
                    'atom_3': columns[4],
                    'atom_4': columns[5]
                }

            elif improper_section:
                improper_dict[columns[1]] = {
                    'atom_1': columns[2],
                    'atom_2': columns[3],
                    'atom_3': columns[4],
                    'atom_4': columns[5]
                }

        if atom_dict and any('x' in atom for atom in atom_dict.values()):
            max_x = max(atom['x'] for atom in atom_dict.values())
            min_x = min(atom['x'] for atom in atom_dict.values())
            max_y = max(atom['x'] for atom in atom_dict.values())
            min_y = min(atom['x'] for atom in atom_dict.values())
            max_z = max(atom['x'] for atom in atom_dict.values())
            min_z = min(atom['x'] for atom in atom_dict.values())
            print(f"The maximum value of x is: {max_x}")
        else:
            print("atom_dict is empty or contains no 'x' values.")
        #print(atom_dict)
        #max_x = max(atom['x'] for atom in atom_dict.values())
        #min_x = min(atom['x'] for atom in atom_dict.values())
        #max_y = max(atom['y'] for atom in atom_dict.values())
        #min_y = min(atom['y'] for atom in atom_dict.values())
        #max_z = max(atom['z'] for atom in atom_dict.values())
        #min_z = min(atom['z'] for atom in atom_dict.values())


        if not any([pair_section, bond_section, angle_section, dihedral_section, improper_section, mass_section, atom_section]):
            line_strip = line.strip()
            line_split = line_strip.split()
            if header_section:
                try:
                    float(line_split[0])
                    if float(line_split[0]) > 0:
                        header[i] = {'info': line}
                    else:
                        header[i] = {line_split[2]: line_split[0],
                                 line_split[3]: line_split[1]}
                except ValueError:
                    header[i] = {'info': line}    
            else:
                footer[i] = {'info': line}

    lammps_dict = {
        'mass_dict': mass_dict,
        'pair_dict': pair_dict,
        'bond_dict': bond_dict,
        'angle_dict': angle_dict,
        'dihedral_dict': dihedral_dict,
        'improper_dict': improper_dict,
        'atom_dict': atom_dict
    }

    if pseudoatoms:
        lammps_dict_ff_form = {}
        for key, value in lammps_dict.items():
            if isinstance(value, dict):
                lammps_dict_ff_form[key] = {}
                for sub_key, sub_value in value.items():
                    if isinstance(sub_value, dict):
                        lammps_dict_ff_form[key][sub_key] = {}
                        for item_key, item_value in sub_value.items():
                            lammps_dict_ff_form[key][sub_key][item_key] = \
                            item_value[:-1] if item_value else item_value
                    else:
                        lammps_dict_ff_form[key][sub_key] = \
                        sub_value[:-1] if sub_value else sub_value
            else:
                lammps_dict_ff_form[key] = value

    else:
        lammps_dict_ff_form = copy.deepcopy(lammps_dict)

    lammps_dict_full = {'lammps_dict': lammps_dict,
                        'lammps_dict_ff_form': lammps_dict_ff_form,
                        'header': header,
                        'footer': footer}

    if linker_identifier is not None:
        for key, value in lammps_dict_full['lammps_dict_ff_form'].items():
            if isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    if isinstance(value, dict):
                        for item_key, item_value in sub_value.items():
                            if str(item_value).startswith(str(linker_identifier)):
                                lammps_dict_full['lammps_dict_ff_form'][key][sub_key][item_key] = \
                                item_value[1:]

    if out_json:
        with open(out_json, 'w', encoding='utf-8') as file:
            json.dump(lammps_dict_full, file, indent=4)
        print(f'printed .json file ({out_json}) containing information from {input_file}')

    print(f'[read_topo] read {input_file}')

    return lammps_dict_full

def read_gaff2(gaff2_source: str,
               verbose: Optional[bool] = True,
               out_json: Optional[str] = None) -> dict:
    '''
    Parse a GAFF2 file to extract atom, bond, angle, and other data sections.
        Format of returned gaff2_dict:
    gaff2_dict
        |--Information
                |--Header \n
                 --Forcefield \n
                 --Version \n
                 --Date \n
                 --Pair_Type \n
                 --Bond_Type \n
                 --Angle_Type \n
                 --Dihedral_Type \n
                 --Improper_Type \n
        |--Atoms
            |--line index
                |--Atom \n
                 --Mass \n
                 --Unknown Param (not used) \n
                 --Information \n
        |--Pairs
            |--line index
                |--Pair \n
                 --r_Pair \n
                 --Atom \n
                 --sigma \n
                 --eps \n
        |--Bonds
            |--line index
                |--Bond \n
                 --r_Bond \n
                 --Atom_1 \n
                 --Atom_2 \n
                 --K \n
                 --r \n
                 --Information \n
        |--Angles
            |--line index
                |--Angle \n
                 -- r_Angle \n
                 --Atom_1 \n
                 --Atom_2 \n
                 --Atom_3 \n
                 --Data_1 \n
                 --Data_2 \n
                 --Info \n
        |--Dihedrals
            |--line index
                |--Atom_1 \n
                 --Atom_2 \n
                 --Atom_3 \n
                 --Atom_4 \n
                 --Data_1 \n
                 --Data_2 \n
                 --Info \n
        |--Impropers
            |--line index
                |--Atom_1 \n
                 --Atom_2 \n
                 --Atom_3 \n
                 --Atom_4 \n
                 --Data_1 \n
                 --Data_2 \n
                 --Info \n
    '''
    gaff2 = {}
    Information = []
    Pairs = {}
    Atoms = {}
    Bonds = {}
    Angles = {}
    Dihedrals = {}
    Impropers = {}

    skipped = False
    problem = False
    please_check = False

    section_counter = 1

    # regex identifiers for these sections because of the wildcard X containing a space before it
    bond_pattern = re.compile(r'([a-zA-Z0-9+]+)\s*[-]\s*([a-zA-Z0-9+]+)\s*(.*)')
    angle_pattern = re.compile(r'([a-zA-Z0-9+]+)\s*[- ]\s*([a-zA-Z0-9+]+)\s*[- ]\s*([a-zA-Z0-9+]+)\s*(.*)')
    dihedral_pattern = re.compile(r'([a-zA-Z0-9+]+)\s*[- ]\s*([a-zA-Z0-9+]+)\s*[- ]\s*([a-zA-Z0-9+]+)\s*[- ]\s*([a-zA-Z0-9+]+)\s*(.*)')
    improper_pattern = re.compile(r'([a-zA-Z0-9+]+)\s*[- ]\s*([a-zA-Z0-9+]+)\s*[- ]\s*([a-zA-Z0-9+]+)\s*[- ]\s*([a-zA-Z0-9+]+)\s*(.*)')

    with open(gaff2_source, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        line = line.strip()
        columns = line.split()

        # iterating through each section recognizing that each section has a line break between them
        # if we find that these breaks change between versions, we can switch to a keyword search
        if len(columns) == 0:
            section_counter += 1
            continue

        if i == 0:  # header info
            Information = {'Header': line,
                           'Forcefield': 'gaff2',
                           'Version': columns[columns.index('(Version')+1],
                           'Date': ' '.join(map(str, columns[columns.index('(Version')+2:])),
                           'Pair_Type': 'lj/charmmfsw/coul/long', # this is meant to replace lj/charmm/coul/long which could also be used if desired
                           'Bond_Type': 'Harmonic',
                           'Angle_Type': 'Harmonic',
                           'Dihedral_Type': 'Fourier',
                           'Improper_Type': 'cvff'}

        if section_counter == 1 and i != 0 :  # atom info
            Atoms[int(i+1)] = {'Atom': columns[0],
                             'Mass': float(columns[1]),
                             'Unknown Param': float(columns[2]), # unsure what this parameter is -> have not found a use for it
                             'Information': ' '.join(map(str, columns[3:]))}

        elif section_counter == 2:  # bond info
            match = bond_pattern.match(line)
            if match:
                bond_data = match.group(3).split()
                Atom_1 = match.group(1)
                Atom_2 = match.group(2)
                Bonds[int(i+1)] = {'Bond': Atom_1+'-'+Atom_2,
                              'r_Bond': Atom_2+'-'+Atom_1,
                                'Atom_1': Atom_1,
                                 'Atom_2': Atom_2,
                                 'K': float(bond_data[0]),
                                 'r': float(bond_data[1]),
                                 'Information': ' '.join(map(str, bond_data[2:]))}
            else:
                if problem is True:
                    print(problem)
                    print(f'Matching problem with {line}')
                    problem = True

        elif section_counter == 3:  # angle info
            match = angle_pattern.match(line)
            if match:
                angle_data = match.group(4).split()
                Atom_1 = match.group(1)
                Atom_2 = match.group(2)
                Atom_3 = match.group(3)
                Angles[int(i+1)] = {'Angle': Atom_1+'-'+Atom_2+'-'+Atom_3,
                               'r_Angle': Atom_3+'-'+Atom_2+'-'+Atom_1,
                                'Atom_1': Atom_1,
                                 'Atom_2': Atom_2,
                                 'Atom_3': Atom_3,
                                 'K': float(angle_data[0]),
                                 'r': float(angle_data[1]),
                                 'Information': ' '.join(map(str, angle_data[2:]))}
            else:
                print(f'Matching problem with {line}')

        elif section_counter == 4:  # dihedral info
            match = dihedral_pattern.match(line)
            if match:
                dihedral_data = match.group(5).split() 
                # The method we have employed to convert Amber style dihedrals into LAMMPS is through the Fourier style
                # m (the number of summations is  always 1)
                # the first two parameters in  gaff2.dat correspond to K=param_1/param_2
                K = float(dihedral_data[1])/float(dihedral_data[0])
                Atom_1 = match.group(1)
                Atom_2 = match.group(2)
                Atom_3 = match.group(3)
                Atom_4 = match.group(4)
                Dihedrals[int(i+1)] = {'Dihedral': Atom_1+'-'+Atom_2+'-'+Atom_3+'-'+Atom_4,
                                  'r_Dihedral': Atom_4+'-'+Atom_3+'-'+Atom_2+'-'+Atom_1,
                                'Atom_1': Atom_1,
                                 'Atom_2': Atom_2,
                                 'Atom_3': Atom_3,
                                 'Atom_4': Atom_4,
                                 'm': 1,
                                 'K': K,
                                 'n': int(float(dihedral_data[3])),
                                 'd': float(dihedral_data[2]),
                                 'Information': ' '.join(map(str, dihedral_data[4:]))}
            else:
                print(f'Matching problem with {line}')

        elif section_counter == 5:  # improper info

            match = improper_pattern.match(line)
            if match:
                improper_data = match.group(5).split()
                # The method we have employed to convert Amber style impropers into LAMMPS is through the cvff style
                # Amber provides the d-angle as either 0 or 180 which we must convert into cvff_style via 
                # where d = cos(d-angle) (e.g., -1 for 180 degrees and 1 for 0 degrees)
                if float(improper_data[1]) == float(180.0):
                    d = int(-1)
                elif float(improper_data[1]) == float(0):
                    d = int(1)
                else:
                    print('unknown d-value, writing as None')
                    d = None

                Atom_1 = match.group(1)
                Atom_2 = match.group(2)
                Atom_3 = match.group(3)
                Atom_4 = match.group(4)

                Impropers[int(i+1)] = {'Improper': Atom_1+'-'+Atom_2+'-'+Atom_3+'-'+Atom_4,
                                  'r_Improper': Atom_4+'-'+Atom_3+'-'+Atom_2+'-'+Atom_1,
                                 'Atom_1':Atom_1,
                                 'Atom_2': Atom_2,
                                 'Atom_3': Atom_3,
                                 'Atom_4': Atom_4,
                                 'K': float(improper_data[0]),
                                 'n': int(float(improper_data[2])),
                                 'd': d,
                                 'Information': ' '.join(map(str, dihedral_data[3:]))}
            else:
                print(f'Matching problem with {line}')
             
        elif section_counter == 8: # pair info
        
            if len(columns) > 2:
                # To convert Amber style pairs into LAMMPS, the sigma definition must be adjusted
                # Amber defined Lennard-Jones potentials as U(r) = eps*((sigma/r)^12-2*(sigma/r)^6)
                # while LAMMPS defines them as 4*eps*((sigma/r)^12-(sigma/r)^6
                # Amber also provides the radius instead of the diameter for sigma
                # To address these two issues we must multiple the Amber value by 2*2**(-1/6)
                # Credit for this descovery goes to Andrew Jewett of MolTemplate
                # https://github.com/jewettaij/moltemplate/blob/master/moltemplate/amber2lt.py
                Atom_1= columns[0]
                Atom_2 = columns[0]
                sigma = float(columns[1])*2*2**(-1/6)
                Pairs[int(i+1)] = {'Pair': Atom_1+'-'+Atom_2,
                                   'r_Pair': Atom_2+'-'+Atom_1,
                                   'Atom_1': Atom_1,
                                   'Atom_2': Atom_2,
                                   'eps': float(columns[2]),
                                   'sigma': sigma}
                
    gaff2 = {'Information': Information,
             'Atoms': Atoms,
             'Pairs': Pairs,
             'Bonds': Bonds,
             'Angles': Angles,
             'Dihedrals': Dihedrals,
             'Impropers': Impropers}

    if verbose:
        for index, dicts in enumerate(gaff2.keys()):
            keys = list(gaff2[dicts].keys())

            if index != 0:
                for j in range(1, len(keys)):
                    expected_key = int(keys[j-1]) + 1
                    while expected_key < int(keys[j]):
                        expected_key += 1
                        print(f'Missing {int(keys[j])-1} from {gaff2_source} in {dicts} Section')
                        please_check = True
        if please_check:
            print('----------------------------------------------------------')
            print(f'Please check {gaff2_source} file for errors in these line')
            print('----------------------------------------------------------')

        else:
            print(f'No errors found while reading {gaff2_source}!')
            for i, dicts in enumerate(gaff2.keys()):
                if i != 0:
                    keys = list(gaff2[dicts].keys())
                    print(f'Section {dicts} in lines [{keys[0]}-{keys[-1]}]')

        if out_json:
            with open(out_json, 'w', encoding='utf-8') as file:
                json.dump(gaff2, file, indent=4)
            print(f'printed .json file ({out_json}) containing information from {gaff2_source}')

    return gaff2

def read_frcmod(frcmod_file: str,
                out_json: Optional[str] = None) -> dict:
    '''
    Format of returned frcmod_dict:
    frcmod_dict
        |--mass_dict
        |--bond_dict
        |--angle_dict
        |--dihedral_dict
                |--line index
                    |--atom_1 
                    --atom_2
                    --atom_3
                    --atom_4
                    --m (Fourier info 1)
                    --K (Fourier info 2)
                    --n (Fourier info 3)
                    --d (Fourier info 4)
                    --r_atom_1 (replacement atom 1)
                    --r_atom_2 (replacement atom 2)
                    --r_atom_3 (replacement atom 3)
                    --r_atom_4 (replacement atom 4)
                    --penalty_score
        |--improper_dict
                |--line index
                    |--atom_1 
                    --atom_2
                    --atom_3
                    --atom_4
                    --K (Fourier info 1)
                    --d (Fourier info 2)
                    --n (Fourier info 3)
                    --r_atom_1 (replacement atom 1 if provided)
                    --r_atom_2 (replacement atom 2 if provided)
                    --r_atom_3 (replacement atom 3 if provided)
                    --r_atom_4 (replacement atom 4 if provided)
                    --penalty_score (if provided)
                    --Information (if provided)
        |--nonbon_dict
    '''
    with open(frcmod_file, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    mass_section = False
    bond_section = False
    angle_section = False
    dihedral_section = False
    improper_section = False
    nonbon_section = False

    mass_check = False
    bond_check = False
    angle_check = False
    nonbon_check = False

    mass_dict = {}
    bond_dict = {}
    angle_dict = {}
    dihedral_dict = {}
    improper_dict = {}
    nonbon_dict = {}
    frcmod_dict = {}

    for i, line in enumerate(lines):

        if not line:
            continue

        line = line.strip()
        columns = line.split()

        dihedral_pattern = re.compile(r'([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(same)\s+(as)\s+([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+)\s*,\s*(penalty)\s+(score=)\s*([\d.]+)')

        improper_pattern = re.compile(
        r'([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*'
        r'(.*?)(?:\s*([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+)\s*[-\s]\s*([a-zA-Z0-9+]+),\s*penalty score=\s*([\d.]+))?$'
    )


        if line.startswith('MASS'):
            mass_section = True
            bond_section = False
            angle_section = False
            dihedral_section = False
            improper_section = False
            nonbon_section = False

        elif line.startswith('BOND'):
            mass_section = False
            bond_section = True
            angle_section = False
            dihedral_section = False
            improper_section = False
            nonbon_section = False

        elif line.startswith('ANGLE'):
            mass_section = False
            bond_section = False
            angle_section = True
            dihedral_section = False
            improper_section = False
            nonbon_section = False

        elif line.startswith('DIHE'):
            mass_section = False
            bond_section = False
            angle_section = False
            dihedral_section = True
            improper_section = False
            nonbon_section = False

        elif line.startswith('IMPROPER'):
            mass_section = False
            bond_section = False
            angle_section = False
            dihedral_section = False
            improper_section = True
            nonbon_section = False

        elif line.startswith('NONBON'):
            mass_section = False
            bond_section = False
            angle_section = False
            dihedral_section = False
            improper_section = True
            nonbon_section = False

        elif mass_section:
            if not line:
                continue
            mass_dict = {'message': 'NOT YET IMPLEMENTED'}
            if mass_check is False:
                print('mass calculations not yet supported')
                print('please send frcmod file to bctapia@mit.edu to update this!')
                mass_check = True

        elif bond_section:
            if not line:
                continue
            bond_dict = {'message': 'NOT YET IMPLEMENTED'}
            if bond_check is False:
                print('bond calculations not yet supported')
                print('please send frcmod file to bctapia@mit.edu to update this!')
                bond_check = True

        elif angle_section:
            if not line:
                continue
            angle_dict = {'message': 'NOT YET IMPLEMENTED'}
            if angle_check is False:
                print('angle calculations not yet supported')
                print('please send frcmod file to bctapia@mit.edu to update this!')
                angle_check = True

        elif dihedral_section:
            if not line:
                continue
            match = dihedral_pattern.match(line)
            if match:
                K = float(match.group(6))/float(match.group(5))
                Atom_1 = match.group(1)
                Atom_2 = match.group(2)
                Atom_3 = match.group(3)
                Atom_4 = match.group(4)
                rep_Atom_1 = match.group(11)
                rep_Atom_2 = match.group(12)
                rep_Atom_3 = match.group(13)
                rep_Atom_4 = match.group(14)
                dihedral_dict[i+1] = {'Dihedral': Atom_1+'-'+Atom_2+'-'+Atom_3+'-'+Atom_4,
                                  'r_Dihedral': Atom_4+'-'+Atom_3+'-'+Atom_2+'-'+Atom_1,
                                    'atom_1': Atom_1,
                                    'atom_2': Atom_2,
                                    'atom_3': Atom_3,
                                    'atom_4': Atom_4,
                                    'm': 1,
                                    'K': K,
                                    'n': float(match.group(8)),
                                    'd': float(match.group(7)),
                                    'rep_Dihedral': rep_Atom_1+'-'+rep_Atom_2+'-'+rep_Atom_3+'-'+rep_Atom_4,
                                    'rep_r_Dihedral': rep_Atom_4+'-'+rep_Atom_3+'-'+rep_Atom_2+'-'+rep_Atom_1,
                                    'rep_atom_1': rep_Atom_1,
                                    'rep_atom_2': rep_Atom_2,
                                    'rep_atom_3': rep_Atom_3,
                                    'rep_atom_4': rep_Atom_4,
                                    'penalty_score': float(match.group(match.lastindex))}
            else:
                print(f'Matching problem with {line}')


        elif improper_section:
            if not line:
                continue
            match = improper_pattern.match(line)
            if match:

                if float(match.group(6)) == float(180.0):
                    d = int(-1)
                elif float(match.group(6)) == float(0):
                    d = int(1)
                else: 
                    print('unknown d-value, writing as None')
                    d = None
                Atom_1 = match.group(1)
                Atom_2 = match.group(2)
                Atom_3 = match.group(3)
                Atom_4 = match.group(4)

                improper_dict[i+1] = {
                    'Improper': Atom_1+'-'+Atom_2+'-'+Atom_3+'-'+Atom_4,
                    'r_Improper': Atom_4+'-'+Atom_3+'-'+Atom_2+'-'+Atom_1,
                    'atom_1': Atom_1,
                    'atom_2': Atom_2,
                    'atom_3':Atom_3,
                    'atom_4': Atom_4,
                    'K': float(match.group(5)),
                    'd': d,
                    'n': int(float(match.group(7))),
                    'rep_Improper': None,
                    'rep_r_Improper': None,
                    'rep_atom_1': None,
                    'rep_atom_2': None,
                    'rep_atom_3': None,
                    'rep_atom_4': None,
                    'penalty_score': None,
                    'Information': match.group(8)
                    }
                default_check = ['Using the default value', None]

                if match.group(8) not in default_check:
                    text_part = match.group(8).strip().rstrip(')')
                    text_parts = re.split(r'[-,\s]', text_part)
                    text_parts = [part for part in text_parts if part]  # Remove empty strings
                    first = True
                    for index, string in enumerate(text_parts):
                        if len(string) < 3 and first:
                            first_index = index
                            first = False
                    rep_Atom_1 = text_parts[first_index]
                    rep_Atom_2 = text_parts[first_index+1]
                    rep_Atom_3 = text_parts[first_index+2]
                    rep_Atom_4 = text_parts[first_index+3]
                    improper_dict[i+1].update({
                        'rep_Improper': rep_Atom_1+'-'+rep_Atom_2+'-'+rep_Atom_3+'-'+rep_Atom_4,
                        'rep_r_Improper': rep_Atom_4+'-'+rep_Atom_3+'-'+rep_Atom_2+'-'+rep_Atom_1,
                        'rep_atom_1': rep_Atom_1,
                        'rep_atom_2': rep_Atom_2,
                        'rep_atom_3': rep_Atom_3,
                        'rep_atom_4': rep_Atom_4,
                        'penalty_score': float(text_parts[-1]),
                        'Information': ' '.join(text_parts[0:(first_index-1)])
                    })

            else:
                print(f'Matching problem with dihedral line: {line}')

        elif nonbon_section:
            if not line:
                continue
            nonbon_dict = {'NOT YET IMPLEMENTED'}
            if nonbon_check is False:
                print('nonbond calculations not yet supported')
                print('please send frcmod file to bctapia@mit.edu to update this!')
                nonbon_check = True

    frcmod_dict = {'mass_dict': mass_dict,
                   'bond_dict': bond_dict,
                   'angle_dict': angle_dict,
                   'dihedral_dict': dihedral_dict,
                   'improper_dict': improper_dict,
                   'nonbon_dict': nonbon_dict}

    if out_json:
        with open(out_json, 'w', encoding='utf-8') as file:
            json.dump(frcmod_dict, file, indent=4)
        print(f'printed .json file ({out_json}) containing information from {frcmod_file}')

    print('[read_frcmod] completed')
    return frcmod_dict







