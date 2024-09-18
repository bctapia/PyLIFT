import inspect
import shutil
import os
from pathlib import Path
import pylift
from pylift import reader

def pipeline():
    '''
    Updates the GAFF2 file within PyLIFT from the one within your AMBER distribution
    '''

    pylift_loc = Path(inspect.getfile(pylift)).parent
    amber_loc = Path(shutil.which('antechamber')).parent.parent

    gaff2_loc = amber_loc / 'dat' / 'leap' / 'parm' / 'gaff2.dat'

    if gaff2_loc.exists():
        print('Updating gaff2.json in pylift')
        reader.read_gaff2(gaff2_loc, default_loc=False, out_json=pylift_loc / 'ff_data' / 'gaff2.json')
        print('Update successful')
    else:
        print('ANTECHAMBER not found on path, searching default location...')
        gaff2_loc = Path(os.path.expanduser("~")) / 'amber24' / 'dat' / 'leap' / 'parm' / 'gaff2.dat'
        
        if gaff2_loc.exists():
            print('Updating gaff2.json in pylift')
            reader.read_gaff2('gaff2.dat', default_loc=True, out_json=pylift_loc / 'ff_data' / 'gaff2.json')
            print('Update successful')
        else:
            print(f'[ERROR] GAFF2 file not found in expected location: {gaff2_loc}')

if __name__ == '__main__':
    pipeline()
