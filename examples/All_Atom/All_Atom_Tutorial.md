# All Atom Molecule Build
This tutorial covers how to use PyLIFT to build an all-atom LAMMPS input script for the polymer of intrinisic microporisity, PIM-1:

 ![image](../../images/Tutorials/PIM-1_monomer.png)
## Prerequisites
Before using PyLIFT we have to perform a couple of GUI-based tasks:
### 1. ChemDraw
 First, we need to generate a 2D molecule structure of the molecule we want to create. If this is a new molecule, you can build it in [ChemDraw](https://revvitysignals.com/products/research/chemdraw) (or [similar](https://alternativeto.net/software/chemdraw/)) and save it as ```PIM-1.cdxml```

### 2. Avogadro (or similar)
We can now open ```PIM-1.cdxml``` in [Avogadro](https://avogadro.cc/). Because ```PIM-1.cdxml``` is a 2D file, it may be in an unrealistic conformation. We can quickly adjust this via the toolbar: 
```
Extensions > Optimize geometry
```
to generate a physically-realistic 3D PIM-1 visualization:
 ![image](../../images/Tutorials/PIM-1_monomer_Avogadro.png)

We can save this file as both a MOL2 and PDB files: ```PIM-1_monomer.mol2``` and ```PIM-1_monomer.pdb```.

### 3. PyRed
To generate accurate partial charges, we recommend submitting the PDB file you created to [PyRed](https://upjv.q4md-forcefieldtools.org/REDServer-Development/). This process might take a couple of days to generate the final MOL2 file. If you want you can either perform this step, or simply continue with the MOL2 file you generated in Avogadro and we will add less-accurate, but usable, charges later.

## PyLIFT
Now, we are ready to use PyLIFT. First, we open a new Python file in our favorite text editor (I like [Visual Studio Code](https://code.visualstudio.com/)) and import the required modules from PyLIFT
```
from pylift import amber, reader, writer, vmd, utilities, builder
```
Next, we create a function so that all the commands are contained and can be performed with a single call:
```
def pipeline():
```
We need to generate the proper forcefield information for the molecule. We want to use the GAFF2 forcefield so we will use the command:
```
amber.antechamber('PIM-1_monomer.mol2',
                  'PIM-1_monomer_Amber.mol2',
                   forcefield = 'gaff2',
                   charge_method = None,
                   missing_search = 'parmchk2')
```
This command reads in the MOL2 file either from PyRed (if used), or that was created in Avogadro. Amber will generate a MOL2 file called ```CANAL-Me-Me2F_monomer_Amber.mol2``` containing the necessary forcefield information. If we used PyRed, then ```charge_method=None``` because charges were already determined. If we didn't use PyRed then we should specify ```charge_method=abcg2``` or another method available in Antechamber (See Ch. 15 of the [Amber Manual](https://ambermd.org/doc12/Amber24.pdf)).
It is possible that the GAFF2 forcefield will not have every parameter necessaery. By specifying ```missing_search='parmchk2'```, the ```parmchk2``` function available in Antechamber will try to find acceptable replacements and will be available in the  ```missing_ff_params.frcmod``` file.
 
To allow PyLIFT to understand ```CANAL-Me-Me2F_monomer_Amber.mol2```, we need to read it into PyLIFT as 
```
molecule = reader.read_mol2('CANAL-Me-Me2F_Amber.mol2')
```
Now, we have a dictionary named ```molecule``` which contains all of the molecules information. Because PIM-1 is a polymer, it will need to eventually undergo polymerization. pysimm provides Polymatic as a simulated-polymerization program. To use Polymatic later, we have to remove one hydrogen from each linking atom so that there is not an extra hydrogen attached after polymerization.

From inspecting the molecule in Avogadro, we know that the linker atoms are atoms 4, 5, 34, and 35. By removing one hydrogen from each linker, we are also inadvertantly adjusting the overall molecule charge. To maintain a net-neutral molecule, we can tell PyLIFT to uniformly redistribute the charge from each removed hydrogen back into the remaining molecule's atoms. We recommended reading the ```charge_consolidation``` tutorial for a further explanation of ```charge_distribution``` and what other types of distributions are possible! To do all this we run:
```
molecule = builder.remove_h(molecule, 
                        specific_atoms=[4,5,34,35], 
                        num_delete=[1,1,1,1],          
                        charge_distribution='uniform',
                        h_identifiers = ['h','H'])
```
Hydrogens are identified as any atom that starts with an 'h' or 'H'. ```num_delete[i]``` atoms are removed from ```specific_atoms[i]```.

Because we plan on polymerizing later, we need to set which atoms are linkers by running:
```
molecule = builder.assign_linkers(molecule, 
                                linker_atoms=[4,35], 
                                linker_identifier='L')
``` 
We have now identified atoms 4 and 35 as the atoms we want to link together by adding an L in front of each atom type. Note: we do not need to assign '5' and '34' as linkers because Polymatic treats these as "secondary" linkers. Regardless, if we wanted to assign them as linkers, PyLIFT would have no problem doing that.

In MOL2 files, atom types and atom names are different entities. Amber has added the forcefield information as atom types, however, We also want to make sure the atom names include the forcefield information to avoid confusion. Therefore, we run
```
molecule = builder.types_to_names(molecule)
```
To create a skeleton LAMMPS data file using TopoTools we need a MOL2 file containing all of our previous work which we generate by running
```
writer.write_mol2(molecule, 'PIM-1_forTopo.tmp.mol2')
```
We can now read the ```PIM-1_forTopo.tmp.mol2``` file using VMD and TopoTools.
```
vmd.topo_write('PIM-1_forTopo.tmp.mol2',
            'PIM-1_fromTopo.tmp.lmps', 
            bonds=True,
            angles=True,
            dihedrals=True,
            impropers=True)
```
This takes ```PIM-1_forTopo.tmp.mol2```, creates bonds, angles, dihedrals, and improper information and writes all of this to ```PIM-1_fromTopo.tmp.lmps```.

Now that We have a LAMMPS file, we need to read this new information back into a PyLIFT dictionary:
```
molecule = reader.read_topo(topo_in = 'PIM-1_fromTopo.tmp.mol2',
                            pseudoatoms = False)
```
PyLIFT needs to understand the forcefield, which is available as a JSON file in PyLIFT in the ff_data folder
```
gaff2 = utilities.read_json('gaff2.json')
```
If we want to use the ```gaff2.dat``` file in Amber instead, see the ff_update tutorial. 
PyLIFT also needs the missing paramter FRCMOD file from Antechamber
```
missing_params = reader.read_frcmod('missing_ff_params.frcmod')
```
We are now ready to add the forcefield information into the LAMMPS file:
```
molecule = builder.add_ff_params(molecule,
                                gaff2,
                                missing_params,
                                pseudoatoms = False,
                                approx_match = True,
                                user_match = None)
```
```molecule``` is our system, ```gaff2``` is the forcefield, and ```missing_params``` are the missing forcefield parameters. We set ```pseudoatoms=False``` because this is an all-atom simulation. ```approx_match=True``` means that if we have some type of parameter (e.g., ```c3-c3-ca-c3```), PyLIFT will first look for a GAFF2 parameter that matches ```c3-c3-ca-c3``` (or ```c3-ca-c3-c3```) but if that doesn't exist, PyLIFT will look for a gaff2 parameter that matches ```X-c3-ca-X``` (or ```X-ca-c3-X```) as ```X``` denotes a wildcard in Amber. ```user_match``` allows a user to specify any other approximate matches to search for. For example for the parameter ```c3-c3-cy-c3```,  ```X-c3-cy-X``` does not exist, but ```cy``` is just a ```c3``` atom in a square system, thefore if we specify ```user_match = {'cy':'c3'}```, the system will also search for ```X-c3-c3-X```.

We can now ask PyLIFT to write the LAMMPS file:
```
writer.write_lammps(molecule, 'PIM-1.lmps',           
                            comment_style=',')
```
Finally, outside of the function we write
```
if __name__ == '__main__':
    pipeline()
```
which means pipeline() will only run when we run the file from the command line. 

We can save the file as ```all_atom_tutorial.py``` and run the file in the terminal with ```python3 all_atom_tutorial.py```.