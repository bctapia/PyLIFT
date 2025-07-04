"""
pylift.writer

License: The MIT License (MIT)

Copyright 2025 Brandon C. Tapia
"""


def write_mol2(mol2_dict: dict, mol2_out: dict) -> None:
    """
    pylift.writer.write_mol2

    Writes out a MOL2 file from PyLIFT-readable dictionaries

    Arguments:
        mol2_dict (dict): MOL2 information stored from pylift.reader.read_mol2()
        mol2_out (str): name of output mol2 file

    Returns:
        None
    """
    atom_dict = mol2_dict.get("atom_dict")
    bond_dict = mol2_dict.get("bond_dict")
    auxinfo_dict = mol2_dict.get("auxinfo_dict")

    if atom_dict is None or bond_dict is None or auxinfo_dict is None:
        if atom_dict is None:
            print("Missing atom information")
        if bond_dict is None:
            print("Missing bond information")
        if auxinfo_dict is None:
            print("Missing header or footer information")
        print(
            "write_mol2_md_builder is specifically for use \
              with pysimm.apps.md_builder.read_mol2_md_builder"
        )
        print(
            "for the system mol2 read and write please use \
              readmol2 and writemol2 from the pysimm.system module"
        )
        print("exiting...")
        return None

    with open(mol2_out, "w", encoding="utf=8") as file:
        molecule_init = auxinfo_dict["header_info"].get("molecule_init", "")
        res = auxinfo_dict["molecule_info"].get("res", "")
        atom_num = auxinfo_dict["molecule_info"].get("atom_num", "0")
        bond_num = auxinfo_dict["molecule_info"].get("bond_num", "0")
        ident_1 = auxinfo_dict["molecule_info"].get("ident_1", "")
        ident_2 = auxinfo_dict["molecule_info"].get("ident_2", "")
        mol_type = auxinfo_dict["molecule_info"].get("mol_type", "")
        charge_method = auxinfo_dict["molecule_info"].get("charge_method", "")

        file.write(
            f"{molecule_init} \n"
            f"{res} \n"
            f"{int(atom_num)} {int(bond_num):4} {int(ident_1):4} {int(ident_2):4}\n"
            f"{mol_type} \n"
            f"{charge_method} \n\n"
            f"{auxinfo_dict['header_info']['atom_init']}\n"
        )

        for key, value in atom_dict.items():
            file.write(
                f'{key} {value["atom_name"]} {value["x"]:2} {value["y"]:2} {value["z"]:2} {value["atom_type"]:2} {value["molecule_num"]:2} {value["res"]:2} {value["charge"]:2}\n'
            )

        file.write(f"{auxinfo_dict['header_info']['bond_init']}\n")

        for key, value in bond_dict.items():
            file.write(
                f'{key} {value["atom_1"]} {value["atom_2"]:2} {value["bond_type"]:2}\n'
            )

        file.write(f"{auxinfo_dict['header_info']['substructure_init']}\n")
        file.write(f"{auxinfo_dict['substructure_info']['substructure']}")

        print(f"[write_mol2] wrote {mol2_out}")

        return None


def write_lammps(lammps_dict, lammps_out, comment_style=None):
    """
    pylift.writer.write_lammps

    Writes out a LAMMPS file from PyLIFT-readable dictionaries

    Arguments:
        lammps_dict (dict): LAMMPS information stored from pylift.reader.read_lammps()
        lammps_out (str): name of output LAMMPS file
        comment_style (str): if/how parameter types are appended as comments after parameters
            (e.g., comment_style = ',' may append # c3,c3 to a c3-c3 bonding parameter)
            comment_style = ',' is especially useful if using pysimm

    Returns:
        None
    """
    lammps_dict_user_form = lammps_dict.get("lammps_dict")
    lammps_dict_ff_form = lammps_dict.get("lammps_dict_ff_form")
    header_dict = lammps_dict.get("header")
    footer_dict = lammps_dict.get("footer")

    mass_dict = lammps_dict_ff_form.get("mass_dict")
    pair_dict = lammps_dict_ff_form.get("pair_dict")
    bond_dict = lammps_dict_ff_form.get("bond_dict")
    angle_dict = lammps_dict_ff_form.get("angle_dict")
    dihedral_dict = lammps_dict_ff_form.get("dihedral_dict")
    improper_dict = lammps_dict_ff_form.get("improper_dict")
    atom_dict = lammps_dict_ff_form.get("atom_dict")

    mass_dict_user = lammps_dict_user_form.get("mass_dict")
    pair_dict_user = lammps_dict_user_form.get("pair_dict")
    bond_dict_user = lammps_dict_user_form.get("bond_dict")
    angle_dict_user = lammps_dict_user_form.get("angle_dict")
    dihedral_dict_user = lammps_dict_user_form.get("dihedral_dict")
    improper_dict_user = lammps_dict_user_form.get("improper_dict")

    with open(lammps_out, "w", encoding="utf-8") as file:

        file.write(f"{header_dict.get('info')}\n")
        file.write(f"{header_dict.get('num_atoms')} atoms\n")
        file.write(f"{header_dict.get('num_bonds')} bonds\n")
        file.write(f"{header_dict.get('num_angles')} angles\n")
        file.write(f"{header_dict.get('num_dihedrals')} dihedrals\n")
        file.write(f"{header_dict.get('num_impropers')} impropers\n")
        file.write(f"{header_dict.get('num_type_atom')} atom types\n")
        file.write(f"{header_dict.get('num_type_bond')} bond types\n")
        file.write(f"{header_dict.get('num_type_angle')} angle types\n")
        file.write(f"{header_dict.get('num_type_dihedral')} dihedral types\n")
        file.write(f"{header_dict.get('num_type_improper')} improper types\n")
        file.write(f"{header_dict.get('xlo')} {header_dict.get('xhi')} xlo xhi\n")
        file.write(f"{header_dict.get('ylo')} {header_dict.get('yhi')} ylo yhi\n")
        file.write(f"{header_dict.get('zlo')} {header_dict.get('zhi')} zlo zhi\n")

        if len(mass_dict) != 0:
            file.write("\nMasses\n\n")

        for key, value in mass_dict.items():
            file.write(f"{key} {value['mass']}")
            if comment_style is not None:
                a = mass_dict_user[key].get("atom")
                file.write(f" # {a}\n")
            else:
                file.write("\n")

        if len(pair_dict) != 0:
            file.write("\nPair Coeffs\n\n")

        for key, value in pair_dict.items():
            file.write(f"{key} {value['eps']} {value['sigma']}")
            if comment_style is not None:
                a = pair_dict_user[key].get("atom")
                file.write(f" # {a}\n")
            else:
                file.write("\n")

        if len(bond_dict) != 0:
            file.write("\nBond Coeffs\n\n")

        for key, value in bond_dict.items():
            K = bond_dict[key].get("K")
            r = bond_dict[key].get("r")
            file.write(f"{key} {K} {r}")
            if comment_style is not None:
                a_1 = bond_dict_user[key].get("atom_1")
                a_2 = bond_dict_user[key].get("atom_2")
                file.write(f" # {a_1}{comment_style}{a_2}\n")
            else:
                file.write("\n")

        if len(angle_dict) != 0:
            file.write("\nAngle Coeffs\n\n")

        for key, value in angle_dict.items():
            K = angle_dict[key].get("K")
            theta = angle_dict[key].get("theta")
            file.write(f"{key} {K} {theta}")
            if comment_style is not None:
                a_1 = angle_dict_user[key].get("atom_1")
                a_2 = angle_dict_user[key].get("atom_2")
                a_3 = angle_dict_user[key].get("atom_3")
                file.write(f" # {a_1}{comment_style}{a_2}{comment_style}{a_3}\n")
            else:
                file.write("\n")

        if len(dihedral_dict) != 0:
            file.write("\nDihedral Coeffs\n\n")

        for key, value in dihedral_dict.items():
            m = dihedral_dict[key].get("m")
            K = dihedral_dict[key].get("K")
            n = dihedral_dict[key].get("n")
            d = dihedral_dict[key].get("d")
            file.write(f"{key} {m} {K} {n} {d}")
            if comment_style is not None:
                a_1 = dihedral_dict_user[key].get("atom_1")
                a_2 = dihedral_dict_user[key].get("atom_2")
                a_3 = dihedral_dict_user[key].get("atom_3")
                a_4 = dihedral_dict_user[key].get("atom_4")
                file.write(
                    f" # {a_1}{comment_style}{a_2}{comment_style}{a_3}{comment_style}{a_4}\n"
                )
            else:
                file.write("\n")

        if len(improper_dict) != 0:
            file.write("\nImproper Coeffs\n\n")

        for key, value in improper_dict.items():
            K = improper_dict[key].get("K")
            d = improper_dict[key].get("d")
            n = improper_dict[key].get("n")
            file.write(f"{key} {K} {d} {n}")
            if comment_style is not None:
                a_1 = improper_dict_user[key].get("atom_1")
                a_2 = improper_dict_user[key].get("atom_2")
                a_3 = improper_dict_user[key].get("atom_3")
                a_4 = improper_dict_user[key].get("atom_4")
                file.write(
                    f" # {a_1}{comment_style}{a_2}{comment_style}{a_3}{comment_style}{a_4}\n"
                )
            else:
                file.write("\n")

        file.write("\nAtoms # full\n\n")

        for key, value in atom_dict.items():
            molecule = atom_dict[key].get("molecule")
            atom_type = atom_dict[key].get("atom_type")
            charge = atom_dict[key].get("charge")
            x = atom_dict[key].get("x")
            y = atom_dict[key].get("y")
            z = atom_dict[key].get("z")
            comment = atom_dict[key].get("comment")
            file.write(
                f"{key} {molecule} {atom_type} {charge} {x} {y} {z} # {comment}\n"
            )

        for key, value in footer_dict.items():
            if value["info"].split()[0].isdigit() is False:
                file.write("\n")

            file.write(f"{value['info']}\n")

            if value["info"].split()[0] in [
                "Bonds",
                "Angles",
                "Dihedrals",
                "Impropers",
            ]:
                file.write("\n")

    print(f"[write_lammps] wrote {lammps_out}")
