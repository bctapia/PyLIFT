'''
******************************************************************************
pylift.writer module
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
import os
import fnmatch
import re
import importlib
import json
import copy
from typing import Optional

def write_mol2(mol2_dict: dict, output_file: dict) -> None:
    '''
    '''
    atom_dict = mol2_dict.get('atom_dict')
    bond_dict = mol2_dict.get('bond_dict')
    auxinfo_dict = mol2_dict.get('auxinfo_dict')

    if atom_dict is None or bond_dict is None or auxinfo_dict is None:
        if atom_dict is None:
            print('Missing atom information')
        if bond_dict is None:
            print('Missing bond information')
        if auxinfo_dict is None:
            print('Missing header or footer information')
        print('write_mol2_md_builder is specifically for use \
              with pysimm.apps.md_builder.read_mol2_md_builder')
        print('for the system mol2 read and write please use \
              readmol2 and writemol2 from the pysimm.system module')
        print('exiting...')
        return None

    with open(output_file, 'w', encoding='utf=8') as file:
        molecule_init = auxinfo_dict['header_info'].get('molecule_init', '')
        res = auxinfo_dict['molecule_info'].get('res', '')
        atom_num = auxinfo_dict['molecule_info'].get('atom_num', '0')
        bond_num = auxinfo_dict['molecule_info'].get('bond_num', '0')
        ident_1 = auxinfo_dict['molecule_info'].get('ident_1', '')
        ident_2 = auxinfo_dict['molecule_info'].get('ident_2', '')
        mol_type = auxinfo_dict['molecule_info'].get('mol_type', '')
        charge_method = auxinfo_dict['molecule_info'].get('charge_method', '')

        file.write(
            f"{molecule_init} \n"
            f"{res} \n"
            f"{int(atom_num)} {int(bond_num):4} {int(ident_1):4} {int(ident_2):4}\n"
            f"{mol_type} \n"
            f"{charge_method} \n\n"
            f"{auxinfo_dict['header_info']['atom_init']}\n")

        for key, value in atom_dict.items():
            file.write(f'{key} {value["atom_name"]} \
                        {value["x"]:2} {value["y"]:2} \
                        {value["z"]:2} {value["atom_type"]:2} \
                        {value["molecule_num"]:2} {value["res"]:2} \
                        {value["charge"]:2}\n')

        file.write(f"{auxinfo_dict['header_info']['bond_init']}\n")

        for key, value in bond_dict.items():
            file.write(f'{key} {value["atom_1"]} {value["atom_2"]:2} {value["bond_type"]:2}\n')

        file.write(f"{auxinfo_dict['header_info']['substructure_init']}\n")
        file.write(f"{auxinfo_dict['substructure_info']['substructure']}")

        return None

def write_lammps(lammps_dict, lammps_out, comment_style=None):
    '''
    '''
    lammps_dict_user_form = lammps_dict.get('lammps_dict')
    lammps_dict_ff_form = lammps_dict.get('lammps_dict_ff_form')
    header_dict = lammps_dict.get('header')
    footer_dict = lammps_dict.get('footer')


    pair_dict = lammps_dict_ff_form.get('pair_dict')
    bond_dict = lammps_dict_ff_form.get('bond_dict')
    angle_dict = lammps_dict_ff_form.get('angle_dict')
    dihedral_dict = lammps_dict_ff_form.get('dihedral_dict')
    improper_dict = lammps_dict_ff_form.get('improper_dict')

    pair_dict_user = lammps_dict_user_form.get('pair_dict')
    bond_dict_user = lammps_dict_user_form.get('bond_dict')
    angle_dict_user = lammps_dict_user_form.get('angle_dict')
    dihedral_dict_user = lammps_dict_user_form.get('dihedral_dict')
    improper_dict_user = lammps_dict_user_form.get('improper_dict')

    with open(lammps_out, 'w', encoding='utf-8') as file:

        for key, value in header_dict.items():
            file.write(f"{value}\n")

        for key, value in pair_dict.items():
            file.write(f"{key}")
            if comment_style is not None:
                a = pair_dict_user[key].get('atom')
                file.write(f" # {a}\n")
            else:
                file.write('\n')

        file.write('\n Bond Coeffs')

        for key, value in bond_dict.items():
            K = bond_dict[key].get('K')
            r = bond_dict[key].get('r')
            file.write(f"{key} {K} {r}")
            if comment_style is not None:
                a_1 = bond_dict_user[key].get('atom_1')
                a_2 = bond_dict_user[key].get('atom_2')
                file.write(f" # {a_1}{comment_style}{a_2}\n")
            else:
                file.write('\n')

        file.write('\n Angle Coeffs')

        for key, value in angle_dict.items():
            K = angle_dict[key].get('K')
            r = angle_dict[key].get('r')
            file.write(f"{key} {K} {r}")
            if comment_style is not None:
                a_1 = angle_dict_user[key].get('atom_1')
                a_3 = angle_dict_user[key].get('atom_2')
                a_3 = angle_dict_user[key].get('atom_3')
                a_4 = angle_dict_user[key].get('atom_4')
                file.write(f" # {a_1}{comment_style}{a_2}{comment_style}{a_3}\n")
            else:
                file.write('\n')

        file.write('\n Dihedral Coeffs')

        for key, value in dihedral_dict.items():
            m = dihedral_dict[key].get('m')
            K = dihedral_dict[key].get('K')
            n = dihedral_dict[key].get('n')
            d = dihedral_dict[key].get('d')
            file.write(f"{key} {m} {K} {n} {d}")
            if comment_style is not None:
                a_1 = dihedral_dict_user[key].get('atom_1')
                a_2 = dihedral_dict_user[key].get('atom_2')
                a_3 = dihedral_dict_user[key].get('atom_3')
                a_4 = dihedral_dict_user[key].get('atom_4')
                file.write(f" # {a_1}{comment_style}{a_2}{comment_style}{a_3}{comment_style}{a_4}\n")
            else:
                file.write('\n')

        file.write('\n Improper Coeffs')

        for key, value in improper_dict.items():
            K = dihedral_dict[key].get('K')
            d = dihedral_dict[key].get('d')
            n = dihedral_dict[key].get('n')
            file.write(f"{key} {K} {d} {n}")
            if comment_style is not None:
                a_1 = improper_dict_user[key].get('atom_1')
                a_2 = improper_dict_user[key].get('atom_2')
                a_3 = improper_dict_user[key].get('atom_3')
                a_4 = improper_dict_user[key].get('atom_4')
                file.write(f" # {a_1}{comment_style}{a_2}{comment_style}{a_3}{comment_style}{a_4}\n")
            else:
                file.write('\n')

        for key, value in footer_dict.items():
            file.write(f"{value}\n")
