'''
pylift.builder

License: The MIT License (MIT)

Copyright 2025 Brandon C. Tapia
'''
from typing import Optional

def convert_to_pseudo(mol2_dict: dict,
                      h_identifiers: Optional[list[str]] = ['h', 'H']) -> dict:
    '''
    pylift.builder.convert_to_pseudo

    Converts all-atom molecules type names to be able to identify united-atom molecules. 
    Atom types are updated by appending the number of hydrogens onto the end 
    (e.g., a c3 gaff2 type is retyped as c33 if it is connected to three hydrogens)
    
    Arguments:
        mol2_dict (dict) : contains molecule information generated with PyLIFT.reader.read_mol2
        h_identifiers (list[str]) : contains all identifiers hydrogen atoms may be referenced as

    Returns: 
        mol2_dict similar dictionary
    '''
    atoms = []
    bonds = []
    heavy_atom_indices = set()
    hydrogen_indices = set()
    atom_to_bonds = {}
    atom_types = {}

    atom_dict = mol2_dict.get('atom_dict')
    bond_dict = mol2_dict.get('bond_dict')

    for key, value in atom_dict.items():
        atom_id = int(key)
        atom_type = value['atom_type']
        atom_to_bonds[atom_id] = []
        atom_types[atom_id] = atom_type
        atoms.append((atom_id, atom_type))
        if not any(value['atom_type'].startswith(identifier) for identifier in h_identifiers):
            heavy_atom_indices.add(int(key))
        else:
            hydrogen_indices.add(int(key))

    for key, value in bond_dict.items():
        atom1 = int(value['atom_1'])
        atom2 = int(value['atom_2'])
        if atom1 not in atom_to_bonds:
            atom_to_bonds[atom1] = []
        if atom2 not in atom_to_bonds:
            atom_to_bonds[atom2] = []
        atom_to_bonds[atom1].append(atom2)
        atom_to_bonds[atom2].append(atom1)
        bonds.append((atom1, atom2))

    hydrogen_counts = {atom: 0 for atom in heavy_atom_indices}

    for heavy_atom in heavy_atom_indices:
        if heavy_atom in atom_to_bonds:
            for bonded_atom in atom_to_bonds[heavy_atom]:
                if bonded_atom in hydrogen_indices:
                    hydrogen_counts[heavy_atom] += 1

    for key, value in atom_dict.items():
        if int(key) in hydrogen_counts:
            h_count = hydrogen_counts.get(int(key))
            value['atom_type'] = value['atom_type'] +str(h_count)


    mol2_dict['atom_dict'] = atom_dict


    return mol2_dict

def types_to_names(mol2_dict: dict) -> dict:
    '''
    pylift.builder.types_to_names

    Mol2 files have unique columns for atom types and atom names. 
    This function copies the forcefield descriptive atom types to the atom names
    as TopoTools uses atom names to generate LAMMPS skeleton files. 

    Arguments:
        mol2_dict (dict) : contains molecule information generated with PyLIFT.reader.read_mol2
    
    Returns: 
        mol2_dict similar dictionary
    '''

    atom_dict = mol2_dict.get('atom_dict')

    if atom_dict is None:
        print('[types_2_names] No atom dictionary available to parse')
        return None

    for _, value in atom_dict.items():
        if 'atom_type' in value:
            value['atom_name'] = value['atom_type']

    print("[types_to_names] completed")

    return mol2_dict

def remove_h_old(mol2_dict: dict,
             specific_atoms: Optional[list[int]] = None,
             num_delete: Optional[list[int]] = None,
             charge_distribution: Optional[str] = 'heavy',
             h_identifiers: list[str] = ['h', 'H']):
    '''
    pylift.builder.remove_h

    Removes hydrogen atoms. If specific_atoms is not specified, all hydrogens are removed.
    If specific_atoms is specified:
        num_delete[i] amount of hydrogens removed from specific_atoms[i]
        len(num_delete) therefore must equal len(specific_atoms)
    
    Arguments:
        mol2_dict (dict) : contains molecule information generated with PyLIFT.reader.read_mol2
        specific_atoms (list[int]) : if specified, heavy atoms to remove hydrogens from. 
            If None, all hydrogens removed
        charge_distribution (str) : method to distribute the charge of removed hydrogens 
            Options: 'heavy', 'uniform', None  
            'heavy': Charge from removed hydrogens is consolidated into the corresponding heavy atom
                recommended if specific_atoms = None or if adjust_polymer_charges will be used
            'uniform': Charge from removed hydrogens is distributed across all remaining atoms
                recommended if specific_atoms != None and adjust_polymer_charges will not be used
             None: charge is not adjusted (will lead to a molecule with a different overall charge)
                not recommended
        h_identifiers (list[str]) : contains all identifiers hydrogen atom types may start with
        num_delete (list[int]) : if specific_atoms is specified, num_delete specifies 
            the number of hydrogens to remove from each specified atom
        
    Returns: 
        mol2_dict similar dictionary
    '''

    atom_dict = mol2_dict.get('atom_dict')
    bond_dict = mol2_dict.get('bond_dict')

    filtered_atom_dict = {}
    filtered_bond_dict = {}
    hydrogen_indices = set()
    total_hydrogen_charge = 0  # Track total charge of removed hydrogens

    original_overall_charge = mol2_dict['auxinfo_dict']['molecule_info'].get('total_charge')

    if original_overall_charge is None:
        original_overall_charge = sum(float(atom_info['charge']) for atom_info in mol2_dict['atom_dict'].values())

    if charge_distribution not in [None, 'heavy', 'uniform']:
        print('charge_distribution must be None, "heavy", or "uniform"')
        print('exiting...')
        return
    if specific_atoms and not num_delete:
        print('by specifying specific_atoms, you must also specify how many hydrogens to remove via num_delete')
        print('exiting...')
        return
    elif num_delete and not specific_atoms:
        print('by specifying num_delete, you must also specify what heavy atoms to remove them from via specific_atoms')
        print('exiting...')
        return
    elif specific_atoms and num_delete:
        if len(specific_atoms) != len(num_delete):
            print('len(specific_atoms) must equal len(num_delete)')
            print('exiting...')
            return

        atom_hydrogens = {atom: [] for atom in specific_atoms}
        hydrogen_to_heavy = {}

        for bond_id, bond_info in bond_dict.items():
            atom_1 = int(bond_info['atom_1'])
            atom_2 = int(bond_info['atom_2'])

            if atom_1 in specific_atoms and any(atom_dict.get(str(atom_2), {}).get('atom_type', '').startswith(h) for h in h_identifiers):
                atom_hydrogens[atom_1].append(atom_2)
                hydrogen_to_heavy[atom_2] = atom_1
            if atom_2 in specific_atoms and any(atom_dict.get(str(atom_1), {}).get('atom_type', '').startswith(h) for h in h_identifiers):
                atom_hydrogens[atom_2].append(atom_1)
                hydrogen_to_heavy[atom_1] = atom_2

        # Check for excessive hydrogen deletion
        for heavy_atom, hydrogens in atom_hydrogens.items():
            if len(hydrogens) < num_delete[specific_atoms.index(heavy_atom)]:
                print(f'Warning: Attempting to delete more hydrogens ({num_delete[specific_atoms.index(heavy_atom)]}) than available ({len(hydrogens)}) for atom {heavy_atom}')
                num_delete[specific_atoms.index(heavy_atom)] = len(hydrogens)

        hydrogens_to_remove = []
        for heavy_atom, hydrogens in atom_hydrogens.items():
            hydrogens_to_remove.extend(hydrogens[:num_delete[specific_atoms.index(heavy_atom)]])

        if hydrogens_to_remove:
            if charge_distribution == 'heavy':
                # Transfer hydrogen charges to heavy atoms
                for hydrogen in hydrogens_to_remove:
                    heavy = hydrogen_to_heavy[hydrogen]
                    hydrogen_charge = float(atom_dict[str(hydrogen)]['charge'])
                    atom_dict[str(heavy)]['charge'] = str(float(atom_dict[str(heavy)]['charge']) + hydrogen_charge)
                    total_hydrogen_charge += hydrogen_charge

            for key, value in atom_dict.items():
                if not any(value['atom_type'].startswith(identifier) for identifier in h_identifiers) or int(key) not in hydrogens_to_remove:
                    filtered_atom_dict[key] = value
                else:
                    hydrogen_indices.add(int(key))
                    if charge_distribution == 'uniform':
                        total_hydrogen_charge += float(value['charge'])

            for key, value in bond_dict.items():
                if value['atom_1'] not in hydrogen_indices and value['atom_2'] not in hydrogen_indices:
                    filtered_bond_dict[key] = value

            old_to_new_index = {}
            new_index = 1
            for old_index in sorted(filtered_atom_dict.keys(), key=int):
                old_to_new_index[int(old_index)] = new_index
                new_index += 1

            renumbered_atom_dict = {}
            for old_index, atom_info in filtered_atom_dict.items():
                new_index = old_to_new_index[int(old_index)]
                renumbered_atom_dict[str(new_index)] = atom_info

            renumbered_bond_dict = {}
            for bond_id, bond_info in filtered_bond_dict.items():
                atom_1 = int(bond_info['atom_1'])
                atom_2 = int(bond_info['atom_2'])
                new_atom_1 = old_to_new_index[atom_1]
                new_atom_2 = old_to_new_index[atom_2]
                bond_type = bond_info.get('bond_type', '')
                renumbered_bond_dict[bond_id] = {'atom_1': new_atom_1, 'atom_2': new_atom_2, 'bond_type': bond_type}

            mol2_dict['atom_dict'] = renumbered_atom_dict
            mol2_dict['bond_dict'] = renumbered_bond_dict
            mol2_dict['auxinfo_dict']['molecule_info']['atom_num'] = len(renumbered_atom_dict)
            mol2_dict['auxinfo_dict']['molecule_info']['bond_num'] = len(renumbered_bond_dict)

            if charge_distribution == 'uniform' and len(renumbered_atom_dict) > 0:
                uniform_charge_increment = total_hydrogen_charge / len(renumbered_atom_dict)
                for atom in renumbered_atom_dict.values():
                    atom['charge'] = str(float(atom['charge']) + uniform_charge_increment)
            
    else:
        for key, value in atom_dict.items():
            if not any(value['atom_type'].startswith(identifier) for identifier in h_identifiers):
                filtered_atom_dict[key] = value
            else:
                hydrogen_indices.add(int(key))
                if charge_distribution == 'uniform':
                    total_hydrogen_charge += float(value['charge'])

        for key, value in bond_dict.items():
            if value['atom_1'] not in hydrogen_indices and value['atom_2'] not in hydrogen_indices:
                filtered_bond_dict[key] = value

        old_to_new_index = {}
        new_index = 1
        for old_index in sorted(filtered_atom_dict.keys(), key=int):
            old_to_new_index[int(old_index)] = new_index
            new_index += 1

        renumbered_atom_dict = {}
        for old_index, atom_info in filtered_atom_dict.items():
            new_index = old_to_new_index[int(old_index)]
            renumbered_atom_dict[str(new_index)] = atom_info

        renumbered_bond_dict = {}
        for bond_id, bond_info in filtered_bond_dict.items():
            atom_1 = int(bond_info['atom_1'])
            atom_2 = int(bond_info['atom_2'])
            new_atom_1 = old_to_new_index[atom_1]
            new_atom_2 = old_to_new_index[atom_2]
            bond_type = bond_info.get('bond_type', '')
            renumbered_bond_dict[bond_id] = {'atom_1': new_atom_1, 'atom_2': new_atom_2, 'bond_type': bond_type}

        mol2_dict['atom_dict'] = renumbered_atom_dict
        mol2_dict['bond_dict'] = renumbered_bond_dict
        mol2_dict['auxinfo_dict']['molecule_info']['atom_num'] = len(renumbered_atom_dict)
        mol2_dict['auxinfo_dict']['molecule_info']['bond_num'] = len(renumbered_bond_dict)

        if charge_distribution == 'uniform' and len(renumbered_atom_dict) > 0:
            uniform_charge_increment = total_hydrogen_charge / len(renumbered_atom_dict)
            for atom in renumbered_atom_dict.values():
                atom['charge'] = str(float(atom['charge']) + uniform_charge_increment)
    
    # Calculate and update overall charge
    overall_charge = sum(float(atom_info['charge']) for atom_info in mol2_dict['atom_dict'].values())
    mol2_dict['auxinfo_dict']['molecule_info'].update({'total_charge': overall_charge})

    print(original_overall_charge)
    print(overall_charge)
    if original_overall_charge - overall_charge != 0:
        charge_percent_diff = abs(abs(original_overall_charge - overall_charge)/((original_overall_charge + overall_charge)/2)*100)
        print(f'[remove_h] Overall molecule charge differs by {charge_percent_diff:.2e}% from original charge')
        
        if charge_percent_diff < 1E-7:
            print('[remove_h] Rounding error due to floating point arithmetic; Small changes can safely be disregarded')
        
    print("[remove_h] completed")
           
    return mol2_dict

def remove_h(mol2_dict: dict,
             specific_atoms: Optional[list[int]] = None,
             num_delete: Optional[list[int]] = None,
             charge_distribution: Optional[str] = 'heavy',
             h_identifiers: list[str] = ['h', 'H']):
    '''
    pylift.builder.remove_h

    Removes hydrogen atoms. If specific_atoms is not specified, all hydrogens are removed.
    If specific_atoms is specified:
        num_delete[i] amount of hydrogens removed from specific_atoms[i]
        len(num_delete) therefore must equal len(specific_atoms)
    
    Arguments:
        mol2_dict (dict) : contains molecule information generated with PyLIFT.reader.read_mol2
        specific_atoms (list[int]) : if specified, heavy atoms to remove hydrogens from. 
            If None, all hydrogens removed
        charge_distribution (str) : method to distribute the charge of removed hydrogens 
            Options: 'heavy', 'uniform', None  
            'heavy': Charge from removed hydrogens is consolidated into the corresponding heavy atom
                recommended if specific_atoms = None or if adjust_polymer_charges will be used
            'uniform': Charge from removed hydrogens is distributed across all remaining atoms
                recommended if specific_atoms != None and adjust_polymer_charges will not be used
             None: charge is not adjusted (will lead to a molecule with a different overall charge)
                not recommended
        h_identifiers (list[str]) : contains all identifiers hydrogen atom types may start with
        num_delete (list[int]) : if specific_atoms is specified, num_delete specifies 
            the number of hydrogens to remove from each specified atom
        
    Returns: 
        mol2_dict similar dictionary
    '''

    atom_dict = mol2_dict.get('atom_dict')
    bond_dict = mol2_dict.get('bond_dict')

    filtered_atom_dict = {}
    filtered_bond_dict = {}
    hydrogen_indices = set()
    total_hydrogen_charge = 0  # Track total charge of removed hydrogens

    original_overall_charge = mol2_dict['auxinfo_dict']['molecule_info'].get('total_charge')

    if original_overall_charge is None:
        original_overall_charge = sum(float(atom_info['charge']) for atom_info in mol2_dict['atom_dict'].values())

    if charge_distribution not in [None, 'heavy', 'uniform']:
        print('charge_distribution must be None, "heavy", or "uniform"')
        print('exiting...')
        return
    if specific_atoms and not num_delete:
        print('by specifying specific_atoms, you must also specify how many hydrogens to remove via num_delete')
        print('exiting...')
        return
    elif num_delete and not specific_atoms:
        print('by specifying num_delete, you must also specify what heavy atoms to remove them from via specific_atoms')
        print('exiting...')
        return
    elif specific_atoms and num_delete:
        if len(specific_atoms) != len(num_delete):
            print('len(specific_atoms) must equal len(num_delete)')
            print('exiting...')
            return

        atom_hydrogens = {atom: [] for atom in specific_atoms}
        hydrogen_to_heavy = {}

        for bond_id, bond_info in bond_dict.items():
            atom_1 = int(bond_info['atom_1'])
            atom_2 = int(bond_info['atom_2'])

            if atom_1 in specific_atoms and any(atom_dict.get(str(atom_2), {}).get('atom_type', '').startswith(h) for h in h_identifiers):
                atom_hydrogens[atom_1].append(atom_2)
                hydrogen_to_heavy[atom_2] = atom_1
            if atom_2 in specific_atoms and any(atom_dict.get(str(atom_1), {}).get('atom_type', '').startswith(h) for h in h_identifiers):
                atom_hydrogens[atom_2].append(atom_1)
                hydrogen_to_heavy[atom_1] = atom_2

        # Check for excessive hydrogen deletion
        for heavy_atom, hydrogens in atom_hydrogens.items():
            if len(hydrogens) < num_delete[specific_atoms.index(heavy_atom)]:
                print(f'Warning: Attempting to delete more hydrogens ({num_delete[specific_atoms.index(heavy_atom)]}) than available ({len(hydrogens)}) for atom {heavy_atom}')
                num_delete[specific_atoms.index(heavy_atom)] = len(hydrogens)

        hydrogens_to_remove = []
        for heavy_atom, hydrogens in atom_hydrogens.items():
            hydrogens_to_remove.extend(hydrogens[:num_delete[specific_atoms.index(heavy_atom)]])

        if hydrogens_to_remove:
            if charge_distribution == 'heavy':
                # Transfer hydrogen charges to heavy atoms
                for hydrogen in hydrogens_to_remove:
                    heavy = hydrogen_to_heavy[hydrogen]
                    hydrogen_charge = float(atom_dict[str(hydrogen)]['charge'])
                    atom_dict[str(heavy)]['charge'] = str(float(atom_dict[str(heavy)]['charge']) + hydrogen_charge)
                    total_hydrogen_charge += hydrogen_charge

            for key, value in atom_dict.items():
                if not any(value['atom_type'].startswith(identifier) for identifier in h_identifiers) or int(key) not in hydrogens_to_remove:
                    filtered_atom_dict[key] = value
                else:
                    hydrogen_indices.add(int(key))
                    if charge_distribution == 'uniform':
                        total_hydrogen_charge += float(value['charge'])

            for key, value in bond_dict.items():
                if value['atom_1'] not in hydrogen_indices and value['atom_2'] not in hydrogen_indices:
                    filtered_bond_dict[key] = value

            old_to_new_index = {}
            new_index = 1
            for old_index in sorted(filtered_atom_dict.keys(), key=int):
                old_to_new_index[int(old_index)] = new_index
                new_index += 1

            renumbered_atom_dict = {}
            for old_index, atom_info in filtered_atom_dict.items():
                new_index = old_to_new_index[int(old_index)]
                renumbered_atom_dict[str(new_index)] = atom_info

            renumbered_bond_dict = {}
            for bond_id, bond_info in filtered_bond_dict.items():
                atom_1 = int(bond_info['atom_1'])
                atom_2 = int(bond_info['atom_2'])
                new_atom_1 = old_to_new_index[atom_1]
                new_atom_2 = old_to_new_index[atom_2]
                bond_type = bond_info.get('bond_type', '')
                renumbered_bond_dict[bond_id] = {'atom_1': new_atom_1, 'atom_2': new_atom_2, 'bond_type': bond_type}

            mol2_dict['atom_dict'] = renumbered_atom_dict
            mol2_dict['bond_dict'] = renumbered_bond_dict
            mol2_dict['auxinfo_dict']['molecule_info']['atom_num'] = len(renumbered_atom_dict)
            mol2_dict['auxinfo_dict']['molecule_info']['bond_num'] = len(renumbered_bond_dict)

            if charge_distribution == 'uniform' and len(renumbered_atom_dict) > 0:
                uniform_charge_increment = total_hydrogen_charge / len(renumbered_atom_dict)
                for atom in renumbered_atom_dict.values():
                    atom['charge'] = str(float(atom['charge']) + uniform_charge_increment)
            
    else:
        # Handle the case when specific_atoms is None
        atom_hydrogens = {}
        hydrogen_to_heavy = {}

        for bond_id, bond_info in bond_dict.items():
            atom_1 = int(bond_info['atom_1'])
            atom_2 = int(bond_info['atom_2'])

            if any(atom_dict.get(str(atom_2), {}).get('atom_type', '').startswith(h) for h in h_identifiers):
                atom_hydrogens[atom_2] = atom_1
                hydrogen_to_heavy[atom_2] = atom_1
            if any(atom_dict.get(str(atom_1), {}).get('atom_type', '').startswith(h) for h in h_identifiers):
                atom_hydrogens[atom_1] = atom_2
                hydrogen_to_heavy[atom_1] = atom_2

        hydrogens_to_remove = list(atom_hydrogens.keys())

        if hydrogens_to_remove:
            if charge_distribution == 'heavy':
                # Transfer hydrogen charges to heavy atoms
                for hydrogen in hydrogens_to_remove:
                    heavy = hydrogen_to_heavy[hydrogen]
                    hydrogen_charge = float(atom_dict[str(hydrogen)]['charge'])
                    atom_dict[str(heavy)]['charge'] = str(float(atom_dict[str(heavy)]['charge']) + hydrogen_charge)
                    total_hydrogen_charge += hydrogen_charge

            for key, value in atom_dict.items():
                if not any(value['atom_type'].startswith(identifier) for identifier in h_identifiers) or int(key) not in hydrogens_to_remove:
                    filtered_atom_dict[key] = value
                else:
                    hydrogen_indices.add(int(key))
                    if charge_distribution == 'uniform':
                        total_hydrogen_charge += float(value['charge'])

            for key, value in bond_dict.items():
                if value['atom_1'] not in hydrogen_indices and value['atom_2'] not in hydrogen_indices:
                    filtered_bond_dict[key] = value

            old_to_new_index = {}
            new_index = 1
            for old_index in sorted(filtered_atom_dict.keys(), key=int):
                old_to_new_index[int(old_index)] = new_index
                new_index += 1

            renumbered_atom_dict = {}
            for old_index, atom_info in filtered_atom_dict.items():
                new_index = old_to_new_index[int(old_index)]
                renumbered_atom_dict[str(new_index)] = atom_info

            renumbered_bond_dict = {}
            for bond_id, bond_info in filtered_bond_dict.items():
                atom_1 = int(bond_info['atom_1'])
                atom_2 = int(bond_info['atom_2'])
                new_atom_1 = old_to_new_index[atom_1]
                new_atom_2 = old_to_new_index[atom_2]
                bond_type = bond_info.get('bond_type', '')
                renumbered_bond_dict[bond_id] = {'atom_1': new_atom_1, 'atom_2': new_atom_2, 'bond_type': bond_type}

            mol2_dict['atom_dict'] = renumbered_atom_dict
            mol2_dict['bond_dict'] = renumbered_bond_dict
            mol2_dict['auxinfo_dict']['molecule_info']['atom_num'] = len(renumbered_atom_dict)
            mol2_dict['auxinfo_dict']['molecule_info']['bond_num'] = len(renumbered_bond_dict)

            if charge_distribution == 'uniform' and len(renumbered_atom_dict) > 0:
                uniform_charge_increment = total_hydrogen_charge / len(renumbered_atom_dict)
                for atom in renumbered_atom_dict.values():
                    atom['charge'] = str(float(atom['charge']) + uniform_charge_increment)
    
    # Calculate and update overall charge
    overall_charge = sum(float(atom_info['charge']) for atom_info in mol2_dict['atom_dict'].values())
    mol2_dict['auxinfo_dict']['molecule_info'].update({'total_charge': overall_charge})

    if original_overall_charge - overall_charge != 0:
        charge_percent_diff = abs(abs(original_overall_charge - overall_charge)/((original_overall_charge + overall_charge)/2)*100)
        print(f'[remove_h] Overall molecule charge differs by {charge_percent_diff:.2e}% from original charge')
        
        if charge_percent_diff < 1E-7:
            print('[remove_h] Rounding error due to floating point arithmetic; Small changes can safely be disregarded')
        
    print("[remove_h] completed")
           
    return mol2_dict

def assign_linkers(mol2_dict: dict,
                   linker_atoms: list[int],
                   linker_identifier: str = 'L'):

    '''
    pylift.builder.assign_linkers

    Adds an identifier to specific atoms to prepare them for linking in Polymatic
        
    Arguments:
        mol2_dict (dict) : contains molecule information generated with PyLIFT.reader.read_mol2
        linker_atoms (list[int]) : id of atoms to assign identifier to
        linker_identifier (str) : what the identifier is

    Returns: 
        mol2_dict similar dictionary
    '''

    for linker_atom in linker_atoms:
        try:
            atom_name = mol2_dict['atom_dict'][str(linker_atom)].get('atom_name')
            atom_type = mol2_dict['atom_dict'][str(linker_atom)].get('atom_type')
            mol2_dict['atom_dict'][str(linker_atom)]['atom_name'] = \
            str(linker_identifier) + atom_name
            mol2_dict['atom_dict'][str(linker_atom)]['atom_type'] = \
            str(linker_identifier) + atom_type
        except KeyError:
            print('[assign_linkers] Error in assign_linkers:')
            print(f'Linker atom {linker_atom} not found in atom dictionary')
            return None
        
    print("[assign_linkers] completed")

    return mol2_dict

# *need to deal with the issue of adjusting charges if using an amber chargestyle like abcg2
def adjust_charges(monomer_dict: dict,
                          xmer_dict: dict,
                          monomer_atoms: list[int],
                          xmer_atoms: list[int],
                          forced_charge: Optional[float] = float(0),
                          verbose: Optional[bool] = True) -> dict:
    '''
    pylift.builder.adjust_charges

    Adjust charge on atoms involved/near to linkers by using an xmer (e.g., dimer) for some charges 
        (e.g., adjust the charges on the heavy atoms involved in polymerization 
        as well as the hydrogens attached)
    Any charge discrpency is uniformly distributed across all atoms 
        such that overall_charge = forced_charge

    Arguments:
        monomer_dict (dict): monomer read in through reader.read_mol2
        xmer_dict (dict): xmer read in through reader.read_mol2
        monomer_atoms (list[int]): atoms where charge will be adjusted in monomer_dict
        dimer_atoms (list[int]): atoms where charge will be read from in xmer_dict
        forced_charge (float): overall monomer charge after adjusting polymer charges.
            Charge is uniformly distributed across all atoms including monomer_atoms
            If forced_charge=None, forced_charge is the overall_charge existing prior to adjustments

    Returns:
        monomer_dict similar dictionary
    '''
   
    monomer_overall_charge = monomer_dict['auxinfo_dict']['molecule_info'].get('total_charge')
    overall_charge_updated = float(0)

    if forced_charge is not None:
        monomer_overall_charge = forced_charge
    
    for i, monomer_atom in enumerate(monomer_atoms):
        xmer_atom = xmer_atoms[i]
        xmer_atom_dict = xmer_dict['atom_dict'].get(str(xmer_atom))
        monomer_atom_dict = monomer_dict['atom_dict'].get(str(monomer_atom))

        if not xmer_atom_dict:
            print(f'[adjust_polymer_charges] No xmer atom {xmer_atom} found')
            print('[adjust_polymer_charges] no charge adjustment performed')
            return
        if not monomer_atom_dict:
            print(f'[adjust_polymer_charges] No monomer atom {monomer_atom} found')
            print('[adjust_polymer_charges] no charge adjustment performed')
            return
 
        xmer_atom_charge = xmer_atom_dict['charge']
        monomer_dict['atom_dict'][str(monomer_atom)].update({'charge': xmer_atom_charge})

        if verbose:
            print(f"Updated atom {monomer_atom} charge from {monomer_atom_dict.get('charge')} to {xmer_atom_charge}")

    # getting new overall charge information
    for key, value in monomer_dict['atom_dict'].items():
        overall_charge_updated += float(value['charge'])

    uniform_charge_increment = (monomer_overall_charge - overall_charge_updated)/len(monomer_dict['atom_dict'])

    for key, value in monomer_dict['atom_dict'].items():
        value['charge'] = float(value['charge']) + uniform_charge_increment
   
    return monomer_dict

# *need to add user_match
def add_ff_params(lammps_dict: dict,
                  ff_dict: dict,
                  missing_ff_params: dict,
                  pseudoatoms: bool,
                  approx_match: bool = True,
                  user_match: Optional[dict[str]] = None,
                  verbose: Optional[bool] = True):
    '''
    pylift.builder.add_ff_params

    Adds forcefield parameters from a specified forcefield to 
    pair, bond, angle, dihedral, and improper topology information
    
    Arguments:
        lammps_dict (dict) : contains skeleton topology for molecule
        ff_dict (dict) : contains forcefield to apply
        missing_ff_params (dict) : contains FRCMOD information
        pseudoatoms (bool) : whether the simulation is all-atom (False) or united-atom (True) 
        approx_match (bool) : whether only exact matches to forcefield types should be found 
            (e.g., if  approx_match = True then c3-c3-ca-c3 would first look 
            for parameters for c3-c3-ca-c3 and if not found, then X-c3-ca-X)
        user_match (dict[str]) : contains a dictionary of current atoms (key) to atom to search 
            for instead (value) (e.g., if {'c3', 'cy'} is specified then X-c3-cy-X would first 
            look for parameters for X-c3-cy-X and if not found, then X-c3-c3-X)
    
    Returns:
        dict: lammps_dict like dictionary with forcefield parameters
    '''

    if user_match is not None:
        print('user_match is not yet implemented. Please contact Brandon (bctapia@mit.edu) for faster implementation if needed')

    lammps_dict_user_form = lammps_dict.get('lammps_dict')

    lammps_dict_ff_form = lammps_dict.get('lammps_dict_ff_form')

    atoms_dict_with_h = lammps_dict_user_form.get('mass_dict')

    atoms_dict = lammps_dict_ff_form.get('mass_dict')
    pair_dict = lammps_dict_ff_form.get('pair_dict')
    bond_dict = lammps_dict_ff_form.get('bond_dict')
    angle_dict = lammps_dict_ff_form.get('angle_dict')
    dihedral_dict = lammps_dict_ff_form.get('dihedral_dict')
    improper_dict = lammps_dict_ff_form.get('improper_dict')

    ff_atoms = ff_dict.get('Atoms')
    ff_pairs = ff_dict.get('Pairs')
    ff_bonds = ff_dict.get('Bonds')
    ff_angles = ff_dict.get('Angles')
    ff_dihedrals = ff_dict.get('Dihedrals')
    ff_impropers = ff_dict.get('Impropers')

    if missing_ff_params is not None:
        missing_ff_mass = missing_ff_params.get('mass_dict')
        missing_ff_pairs = missing_ff_params.get('nonbon_dict')
        missing_ff_bonds = missing_ff_params.get('bond_dict')
        missing_ff_angles = missing_ff_params.get('angle_dict')
        missing_ff_dihedrals = missing_ff_params.get('dihedral_dict')
        missing_ff_impropers = missing_ff_params.get('improper_dict')
    else:
        print('[add_ff_params] missing_ff_params dictionary was not supplied')

    if pseudoatoms:
        print('[add_ff_params] pseudoatoms=True, unable to add pair interactions...')

    if verbose:
        print('\n============[add_ff_params] Starting forcefield search for parameters============')

    for key, value in pair_dict.items():
        pair = value['atom']
        found = False
        for gaff_key, gaff_value in ff_pairs.items():
            if pair == gaff_value['Atom_1']:
                if pseudoatoms:
                    eps = 'epsilon_placeholder'
                    sigma = 'sigma_placeholder'
                    found = True
                    break
                else:
                    eps = gaff_value['eps']
                    sigma = gaff_value['sigma']
                    found = True
                    break
        if found is False:
            for frcmod_key, frcmod_value in missing_ff_pairs.items():
                if pair == frcmod_value['atom']:
                    eps = frcmod_value['eps']
                    sigma = frcmod_value['sigma']
                    found = True
                    if verbose:
                        print(f"FRCMOD REPLACEMENT FOUND: Atom {pair} found in frcmod as {frcmod_value['atom']}")
                    found = True
                    break
        if found:
            lammps_dict['lammps_dict_ff_form']['pair_dict'][key].update({'eps': eps,
                                                                        'sigma': sigma})
        elif verbose:
            print(f"FF values for Pair {pair} not found in ff_dict; please populate manually")

    for key, value in atoms_dict.items():
        atom = value['atom']
        found = False
        for gaff_key, gaff_value in ff_atoms.items():
            if atom == gaff_value['Atom']:
                if pseudoatoms:
                    atom_with_h = atoms_dict_with_h[key]['atom']
                    num_h = int(atom_with_h[-1])
                    mass = gaff_value['Mass']+num_h*1.008
                    found = True
                    break
                else: 
                    mass = gaff_value['Mass']
                    found = True
                    break
        if found is False:
            for frcmod_key, frcmod_value in missing_ff_mass.items():
                if atom in frcmod_value['atom']:
                    mass = frcmod_value['mass']
                    found = True
                    if verbose:
                        print(f"FRCMOD REPLACEMENT FOUND: Mass {mass} found in frcmod as {frcmod_value['mass']}")
                    found = True
                    break
        if found:
            lammps_dict['lammps_dict_ff_form']['mass_dict'][key].update({'mass': mass})
        elif verbose:
            print(f"Mass values for Atom {atom} not found in ff_dict, keeping TopoTools mass")

    for key, value in bond_dict.items():
        bond = str(value['atom_1']+'-'+value['atom_2'])
        found = False
        for gaff_key, gaff_value in ff_bonds.items():
            if bond == gaff_value['Bond'] or bond == gaff_value['r_Bond']:
                K = gaff_value['K']
                r = gaff_value['r']
                found = True
                break
        if found is False:
            for frcmod_key, frcmod_value in missing_ff_bonds.items():
                if bond in (frcmod_value['Bond'], frcmod_value['r_Bond']):
                    K = frcmod_value['K']
                    r = frcmod_value['r']
                    found = True
                    if verbose:
                        print(f"FRCMOD REPLACEMENT FOUND: Bond {bond} found in frcmod as {frcmod_value['rep_Bond']}")
                    found = True
                    break
        if found:
            lammps_dict['lammps_dict_ff_form']['bond_dict'][key].update({'K': K,
                                                                         'r': r})
        elif verbose:
            print(f"FF values for Bond {bond} not found in ff_dict; please populate manually")


    for key, value in angle_dict.items():
        angle = str(value['atom_1']+'-'+value['atom_2']+'-'+value['atom_3'])
        found = False
        for gaff_key, gaff_value in ff_angles.items():
            if angle == gaff_value['Angle'] or angle == gaff_value['r_Angle']:
                K = gaff_value['K']
                theta = gaff_value['theta']
                found = True
                #print(f"Angle: {angle}, {gaff_value['Angle']}, {gaff_value['r_Angle']}, {K}, {r}")
                break
        if found is False:
            for frcmod_key, frcmod_value in missing_ff_angles.items():
                if angle in (frcmod_value['Angle'], frcmod_value['r_Angle']):
                    K = frcmod_value['K']
                    theta = frcmod_value['theta']
                    found = True
                    if verbose:
                        print(f"FRCMOD REPLACEMENT FOUND: Angle {angle} found in frcmod as {frcmod_value['rep_Angle']}")
                    found = True
                    break
        if found:
            lammps_dict['lammps_dict_ff_form']['angle_dict'][key].update({'K': K,
                                                                          'theta': theta})
        elif verbose:
            print(f"FF values for Angle {angle} not found in ff_dict; please populate manually")

    # adding dihdedral forcefield data
    for key, value in dihedral_dict.items():
        dihedral = str(value['atom_1']+'-'+value['atom_2']+'-'+value['atom_3']+'-'+value['atom_4'])
        if approx_match:
            approx_dihedral = str('X-'+value['atom_2']+'-'+value['atom_3']+'-X')
        else:
            approx_dihedral = 'NOT A VALUE TO LOOK UP'
        found = False
        for gaff_key, gaff_value in ff_dihedrals.items():
            if dihedral in (gaff_value['Dihedral'], gaff_value['r_Dihedral']) or \
                        approx_dihedral in (gaff_value['Dihedral'], gaff_value['r_Dihedral']):
                m = gaff_value['m']
                K = gaff_value['K']
                n = gaff_value['n']
                d = gaff_value['d']
                found = True
                #  print(f"Dihedral: {dihedral}, {gaff_value['Dihedral']}, {gaff_value['r_Dihedral']}, {m}, {K}, {n}, {d}")
                break
        if found is False:
            for frcmod_key, frcmod_value in missing_ff_dihedrals.items():
                if (dihedral in (frcmod_value['Dihedral'], frcmod_value['r_Dihedral']) or \
                    approx_dihedral in (frcmod_value['Dihedral'], frcmod_value['r_Dihedral'])):
                    m = frcmod_value['m']
                    K = frcmod_value['K']
                    n = frcmod_value['n']
                    d = frcmod_value['d']
                    found = True
                    if verbose:
                        print(f"FRCMOD REPLACEMENT FOUND: Dihedral {dihedral} found in frcmod as {frcmod_value['rep_Dihedral']}")
                    #  print('found replacement:')
                    #  print(f"Dihedral: {dihedral}, {frcmod_value['Dihedral']}, {frcmod_value['rep_Dihedral']}, {m}, {K}, {n}, {d}")
                    found = True
                    break
        if found:
            lammps_dict['lammps_dict_ff_form']['dihedral_dict'][key].update({'m': m,
                                                                             'K': K,
                                                                             'n': int(n),
                                                                             'd': d})
        elif verbose:
            print(f"FF values for Dihedral {dihedral} not found in ff_dict")

    for key, value in improper_dict.items():
        improper = str(value['atom_1']+'-'+value['atom_2']+'-'+value['atom_3']+'-'+value['atom_4'])

        if approx_match:
            approx_improper = str('X-'+value['atom_2']+'-'+value['atom_3']+'-X')
        else:
            approx_improper = 'NOT A VALUE TO LOOK UP'

        found = False

        for gaff_key, gaff_value in ff_impropers.items():
            if improper in (gaff_value['Improper'], gaff_value['r_Improper']) or \
                        approx_improper in (gaff_value['Improper'], gaff_value['r_Improper']):
                K = gaff_value['K']
                n = gaff_value['n']
                d = gaff_value['d']
                found = True
                # print(f"Improper: {improper}, {gaff_value['Improper']}, {gaff_value['r_Improper']}, {K}, {n}, {d}")
                break

        if found is False:
            for frcmod_key, frcmod_value in missing_ff_impropers.items():
                if (improper in (frcmod_value['Improper'], frcmod_value['r_Improper']) or \
                        approx_improper in (frcmod_value['Improper'], frcmod_value['r_Improper'])):
                    K = frcmod_value['K']
                    n = frcmod_value['n']
                    d = frcmod_value['d']
                    found = True
                    if verbose:
                        if frcmod_value['rep_Improper']:
                            print(f"FRCMOD REPLACEMENT FOUND: Improper {improper} found in frcmod as {frcmod_value['rep_Improper']}")
                        else:
                            print(f"FRCMOD REPLACEMENT FOUND: Improper {improper} found in frcmod as default value")
                    # print('found replacement:')
                    # print(f"Improper: {improper}, {frcmod_value['Improper']}, {frcmod_value['rep_Improper']}, {K}, {n}, {d}")
                    found = True
                    break
        if found:
            lammps_dict['lammps_dict_ff_form']['improper_dict'][key].update({'K': K,
                                                                          'n': n,
                                                                          'd': d})
        elif verbose:
            print(f"FF values for Improper {improper} not found in ff_dict")

    if verbose:
        print('=============[add_ff_params] Ending forcefield search for parameters=============\n')


    return lammps_dict
