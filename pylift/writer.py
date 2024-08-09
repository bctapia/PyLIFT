'''
pylift.writer module

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
            file.write(f'{key} {value["atom_name"]} {value["x"]:2} {value["y"]:2} {value["z"]:2} {value["atom_type"]:2} {value["molecule_num"]:2} {value["res"]:2} {value["charge"]:2}\n')

        file.write(f"{auxinfo_dict['header_info']['bond_init']}\n")

        for key, value in bond_dict.items():
            file.write(f'{key} {value["atom_1"]} {value["atom_2"]:2} {value["bond_type"]:2}\n')

        file.write(f"{auxinfo_dict['header_info']['substructure_init']}\n")
        file.write(f"{auxinfo_dict['substructure_info']['substructure']}")

        print(f"[write_mol2] wrote {output_file}")

        return None

def write_lammps(lammps_dict, lammps_out, comment_style=None):
    '''
    '''
    lammps_dict_user_form = lammps_dict.get('lammps_dict')
    lammps_dict_ff_form = lammps_dict.get('lammps_dict_ff_form')
    header_dict = lammps_dict.get('header')
    footer_dict = lammps_dict.get('footer')

    mass_dict = lammps_dict_ff_form.get('mass_dict')
    pair_dict = lammps_dict_ff_form.get('pair_dict')
    bond_dict = lammps_dict_ff_form.get('bond_dict')
    angle_dict = lammps_dict_ff_form.get('angle_dict')
    dihedral_dict = lammps_dict_ff_form.get('dihedral_dict')
    improper_dict = lammps_dict_ff_form.get('improper_dict')
    atom_dict = lammps_dict_ff_form.get('atom_dict')

    mass_dict_user = lammps_dict_user_form.get('mass_dict')
    pair_dict_user = lammps_dict_user_form.get('pair_dict')
    bond_dict_user = lammps_dict_user_form.get('bond_dict')
    angle_dict_user = lammps_dict_user_form.get('angle_dict')
    dihedral_dict_user = lammps_dict_user_form.get('dihedral_dict')
    improper_dict_user = lammps_dict_user_form.get('improper_dict')

    with open(lammps_out, 'w', encoding='utf-8') as file:

        file.write(f"{header_dict.get('info')}\n")
        file.write(f"{header_dict.get('num_atoms')} atoms\n")
        file.write(f"{header_dict.get('num_bonds')} bonds\n")
        file.write(f"{header_dict.get('num_angles')} angles\n")
        file.write(f"{header_dict.get('num_dihedrals')} dihedrals\n")
        file.write(f"{header_dict.get('num_impropers')} impropers\n")
        file.write(f"{header_dict.get('num_type_atom')} atom types\n")
        file.write(f"{header_dict.get('num_type_bond')} bond types\n")
        file.write(f"{header_dict.get('num_type_angle')} angle types\n")
        file.write(f"{header_dict.get('num_type_dihedral')} dihderal types\n")
        file.write(f"{header_dict.get('num_type_improper')} improper types\n")
        file.write(f"{header_dict.get('xlo')} {header_dict.get('xhi')} xlo xhi\n")
        file.write(f"{header_dict.get('ylo')} {header_dict.get('yhi')} ylo yhi\n")
        file.write(f"{header_dict.get('zlo')} {header_dict.get('zhi')} zlo zhi\n")
        
        file.write("\nMasses\n")
        for key, value in mass_dict.items():
            file.write(f"{key} {value['mass']}")
            if comment_style is not None:
                a = mass_dict_user[key].get('atom')
                file.write(f" # {a}\n")
            else:
                file.write('\n')        

        file.write('\nPair Coeffs\n')
        for key, value in pair_dict.items():
            file.write(f"{key} {value['eps']} {value['sigma']}")
            if comment_style is not None:
                a = pair_dict_user[key].get('atom')
                file.write(f" # {a}\n")
            else:
                file.write('\n')

        file.write('\nBond Coeffs\n')

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

        file.write('\nAngle Coeffs\n')

        for key, value in angle_dict.items():
            K = angle_dict[key].get('K')
            theta = angle_dict[key].get('theta')
            file.write(f"{key} {K} {theta}")
            if comment_style is not None:
                a_1 = angle_dict_user[key].get('atom_1')
                a_2 = angle_dict_user[key].get('atom_2')
                a_3 = angle_dict_user[key].get('atom_3')
                file.write(f" # {a_1}{comment_style}{a_2}{comment_style}{a_3}\n")
            else:
                file.write('\n')

        file.write('\nDihedral Coeffs\n')

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

        file.write('\nImproper Coeffs\n')

        for key, value in improper_dict.items():
            K = improper_dict[key].get('K')
            d = improper_dict[key].get('d')
            n = improper_dict[key].get('n')
            file.write(f"{key} {K} {d} {n}")
            if comment_style is not None:
                a_1 = improper_dict_user[key].get('atom_1')
                a_2 = improper_dict_user[key].get('atom_2')
                a_3 = improper_dict_user[key].get('atom_3')
                a_4 = improper_dict_user[key].get('atom_4')
                file.write(f" # {a_1}{comment_style}{a_2}{comment_style}{a_3}{comment_style}{a_4}\n")
            else:
                file.write('\n')

        file.write('\nAtoms # full\n')

        for key, value in atom_dict.items():
            molecule = atom_dict[key].get('molecule')
            atom_type = atom_dict[key].get('atom_type')
            charge = atom_dict[key].get('charge')
            x = atom_dict[key].get('x')
            y = atom_dict[key].get('y')
            z = atom_dict[key].get('z')
            comment = atom_dict[key].get('comment')
            file.write(f"{key} {molecule} {atom_type} {charge} {x} {y} {z} # {comment}\n")

        for key, value in footer_dict.items():
            if value['info'].split()[0].isdigit() is False:
                file.write('\n')
            file.write(f"{value['info']}\n")
        
    print(f"[write_lammps] wrote {lammps_out}")
