from pylift import vmd
from pylift import reader
from pylift import utilities
from pylift import builder
from pylift import writer

#from pysimm import forcefield
#from pysimm import system
#from pysimm import amber

def pipeline():

    # forcefield to use
    #ff = forcefield.Gaff2()
    # mol2 output from quantum chemistry program (verified with Gaussian) or PyRed
    #polymer_AA = system.read_mol2('CANAL-Me-Me2F_H.mol2')

    #amber.get_forcefield_types(polymer_AA, types='gaff2', f=ff, fo_type='mol2',
                    #mol2_name='CANAL-Me-Me2F_H_gaff2.mol2') # finding the correct gaff2 file types
    #amber.get_missing_ff_params(types='gaff2') # finding any missing forcefield parameters

    molecule = reader.read_mol2(input_file='CANAL-Me-Me2F_H.mol2')
    molecule = builder.remove_h(mol2_dict=molecule, h_identifiers=['h', 'H'], specific_atoms=[1, 37], num_delete=[1,1])
    molecule = builder.types_to_names(mol2_dict=molecule)
    molecule = builder.assign_linkers(linker_atoms=[1, 37], mol2_dict=molecule)