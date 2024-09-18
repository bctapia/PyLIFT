from pylift import amber, reader, writer, vmd, utilities, builder
import numpy as np

def pipeline():

    amber.antechamber('PIM-1_monomer.mol2',
                  'PIM-1_monomer_Amber.mol2',
                   forcefield = 'gaff2',
                   charge_method = None,
                   missing_search = 'parmchk2')

    molecule = reader.read_mol2('PIM-1_monomer_Amber.mol2')

    molecule = builder.convert_to_pseudo(molecule,
                             h_identifiers = ['h', 'H'])

    molecule = builder.remove_h(molecule,
                    charge_distribution='heavy',
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

    molecule = reader.read_topo('PIM-1_fromTopo.tmp.lmps',
                                pseudoatoms=True)

    gaff2 = utilities.read_json('gaff2.json')
    #gaff2 = reader.read_gaff2()

    missing_params = reader.read_frcmod('missing_ff_params.frcmod')

    molecule = builder.add_ff_params(molecule,
                                    gaff2,
                                    missing_params,
                                    pseudoatoms=True,
                                    approx_match = True,
                                    user_match = None)

    writer.write_lammps(molecule, 'PIM-1.lmps', comment_style=',')

    #utilities.cleanup_pylift()

if __name__ == '__main__':
    pipeline()
