# Updating Forcefield File from Amber Source
PyLIFT has the GAFF and GAFF2 forcefields from Amber available in the ```PyLIFT/pylift/ff_data``` directory in JSON files. This data is structured as dictionaries, allowing for current PyLIFT functionality while also being readily extensible if needed. To access this data, we can run
```
from pylift import utilities

gaff2 = utilities.read_json('gaff2.json')
```
which now stores ```gaff2.json``` as a dictionary in the ```gaff2``` variable. Note that ```utilities.read_json``` is configured to always look for data in the ```ff_data``` directory. However, we might worry that ```gaff2.json``` is not the most up-to-date forcefield packaged in Amber. To account for this, PyLIFT provides the ```read_gaff2``` function to read in the classic ```gaff2.dat``` file in Amber and convert it into the dictionary that PyLIFT can understand which we can access via
```
from pylift import reader

gaff2 = reader.read_gaff2('gaff2.dat', default_loc=True)
```
which will read the gaff2.dat file straight from Amber. This of course, requires that Amber be installed on your system. By specifying ```default_loc=True```, we are telling PyLIFT that we expect ```gaff2.dat``` to be in the ```~/amber24/dat/leap/parm``` directory as this is the default Amber location. If this is not true for you, then simply specify ```reader.read_gaff2('path/to/your/gaff2.dat', default_loc=False)``` Provided that the ```gaff2.json``` file in PyLIFT is up-to-date, the following code is interchangeable:

```
gaff2 = reader.read_gaff2('gaff2.dat', default_loc=True)
gaff2 = utilities.read_json('gaff2.json')
```
If you would like, it is also possible to create an updated JSON by running 

```
gaff2 = reader.read_gaff2('gaff2.dat',
                            default_loc=True, 
                            out_json='gaff2_updated.json')
```
Then, moving ```gaff2_updated.json``` into ```ff_data```, we can read it into the system as 
```
gaff2 = utilities.read_json('gaff2_updated.json')
```