import inspect
import shutil
import os
from pathlib import Path
import pylift
from pylift import reader


def pipeline():
    """
    Updates the GAFF file within PyLIFT from the one within your AMBER distribution
    """

    pylift_loc = Path(inspect.getfile(pylift)).parent
    amber_loc = Path(shutil.which("antechamber")).parent.parent

    gaff_loc = amber_loc / "dat" / "leap" / "parm" / "gaff.dat"

    if gaff_loc.exists():
        reader.read_gaff2(
            gaff_loc, default_loc=False, out_json=pylift_loc / "ff_data" / "gaff.json"
        )
    else:
        print("ANTECHAMBER not found on path, searching default location...")
        gaff_loc = (
            Path(os.path.expanduser("~"))
            / "amber24"
            / "dat"
            / "leap"
            / "parm"
            / "gaff.dat"
        )

        if gaff_loc.exists():
            reader.read_gaff2(
                "gaff.dat",
                default_loc=True,
                out_json=pylift_loc / "ff_data" / "gaff.json",
            )
        else:
            print(f"[ERROR] GAFF file not found in expected location: {gaff_loc}")


if __name__ == "__main__":
    pipeline()
