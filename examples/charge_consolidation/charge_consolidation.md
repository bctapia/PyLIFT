# Understanding Uniform vs Heavy Charge Consolidation
PyLIFT provides two functions, ```builder.remove_h()``` and ```adjust_polymer_charges()``` that are able to adjust the partial charges within the molecule and/or full molecule charge. This tutorial walks through the different options available for how to adjust charges when needed.
## ```builder.remove_h()```
```builder.remove_h()``` as the name suggests, provides functionality to remove hydrogen atoms from a molecule. By running:

```
from pylift import builder

# ...previous code to load molecule here...


molecule = builder.remove_h(molecule,
                            charge_distribution='heavy', 
                            h_identifiers=['h', 'H']):

# or something like 

molecule = builder.remove_h(molecule,
                            specific_atoms: [1, 35],
                            num_delete: [1, 1],             
                            charge_distribution='heavy')
```
I am either removing all or specific atoms (see ```all_atom``` or ```united_atom``` tutorials for an in-depth understanding of ```remove_h()```).

Here, we are focusing on what ```charge_distribution``` means. If we specify ```charge_distribution='heavy'```, the charge of each hydrogen that is removed is bundled into the respective heavy atom. Alternatively, we could specify ```charge_distribution='uniform'``` where the charge is spread across all remaining atoms. Regardless, the overall charge of the output molecule will be identical to the overall charge of the input molecule. Lastly, we could specify ```charge_distribution=None``` where no charge consolidation will occur and the outlet molecule will keep the charges of remaining atoms unchanged. This will affect the overall charge of the molecule.

### When to use 'Heavy', 'Uniform', or 'None'?
The choice of charge distribution method is best considered based upon where in the code it shows up. A ```heavy``` distribution minimizes unwanted charge distortion by affecting a minimal number of atoms, however, the charge on these atoms can be shifted significantly. A ```uniform``` distribution affects many atoms, however, the charge distortion is minimal. A ```none``` distribution affects no atoms, however, the overall molecule charge will be shifted. 

While it is up to the user to determine the best methods for themselves, we recommend the following:

1. If attempting to make a united-atom molecule from an all-atom equivalent, perform a ```heavy``` distribution as the heavy atom is now acting as a pseudoatom and should take into account the (now) implicit hydrogen charges.

2. If removing specific hydrogens (e.g., from linker atoms), perform a ```heavy``` distribution if you plan on further addressing partial charges later (such as with ```pylift.builder.adust_charges()```). We recommend this because it limits unecessary distortion to additional atoms as the heavy distribution becomes irrelevant when ```adjust_charges()``` overwrites the charge and automatically performs a ```uniform``` distribution later on.

3. If removing specific atoms (e.g., from linker atoms), perform a ```uniform``` distribution if you do not plan on further addressing partial charges later. We recommend this because there is no physical reason to bundle charge into the heavy atom as it is not attemping to become a pseudoatom, rather, we are simply attempting to remediate an unwanted charge descrepency. By using a ```uniform``` distribution, each atom is minimally affected.

4. We never recommend using ```none``` as a distribution because the molecule's charge will change. Of course, this can be remedied later (e.g., we could have specified ```none``` instead of ```heavy``` in bulletpoint 2), however, there are many potential problems that arise due to charge issues that can easily be avoided by avoiding using ```none```. 

## builder.adjust_charges()
```builder.adjust_charges()``` provides a method to refine charges on any atoms. For example, when determining charges, a monomer is typically used as it can be a time and resource itensive step, especially if using QM methods. However, the charges that are outputted may be significantly different than what is expected. This is especially noticeable on linker atoms where the monomer cannot predict the charge affect of polymerization.

By also knowing the partial charges within an xmer (e.g., a dimer), the charges on the linker atoms are better approximated. By reading in the monomer and the xmer via ```pylift.reader_read_mol2()``` it is then possible to transfer the linker atom charges (and any other atom charges) from the xmer to the monomer via ```adjust_charges()```

```
def adjust_charges(monomer_dict: monomer,
                    xmer_dict: xmer,
                    monomer_atoms: [1, 34],
                    xmer_atoms: [4, 56])
```
Because the adjustment of atom charges will affect the molecule charge, ```adjust_charges()``` will uniformly distribute the charge discrepancy across the entire molecule. We can also specify ```force_charge``` to uniformly distribute a charge to uniformly reach the desired ```forced_charge```. 