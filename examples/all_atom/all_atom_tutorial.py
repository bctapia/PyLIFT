from pylift import amber, reader, writer, vmd, utilities, builder

def pipeline():

    amber.antechamber('PIM-1_monomer.mol2',
                    'PIM-1_monomer_Amber.mol2',
                    forcefield = 'gaff2',
                    charge_method = None,
                    missing_search = 'parmchk2')

    molecule = reader.read_mol2('PIM-1_monomer_Amber.mol2')

    molecule = builder.remove_h(molecule,
                            specific_atoms=[4,5,34,35],
                            num_delete=[1,1,1,1],
                            charge_distribution='uniform',
                            h_identifiers = ['h','H'])

    molecule = builder.assign_linkers(molecule,
                                linker_atoms=[4,35],
                                linker_identifier='L')

    molecule = builder.types_to_names(molecule)

    writer.write_mol2(molecule, 'PIM-1_forTopo.tmp.mol2')

    vmd.topo_write('PIM-1_forTopo.tmp.mol2',
                'PIM-1_fromTopo.tmp.lmps', 
                bonds=True,
                angles=True,
                dihedrals=True,
                impropers=True)

    molecule = reader.read_topo(topo_in = 'PIM-1_fromTopo.tmp.lmps',
                                pseudoatoms = False)

    #gaff2 = utilities.read_json('gaff2.json')
    gaff2 = reader.read_gaff2(out_json='gaff2.json')

    missing_params = reader.read_frcmod('missing_ff_params.frcmod')

    molecule = builder.add_ff_params(molecule,
                                    gaff2,
                                    missing_params,
                                    pseudoatoms = False,
                                    approx_match = True,
                                    user_match = None)
    
    writer.write_lammps(molecule, 'lammps_file.lmps', comment_style=',')

    utilities.cleanup_pylift()

if __name__ == '__main__':
    pipeline()
