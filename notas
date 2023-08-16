#calcular DVH
#esto son notas de codigo random
#pydvh is a package for generating and reading DVH exported from VARIAN Eclipse treatment planning system with Python.

#Presently, pyDVH parses VARIAN Eclipse DVH files. It has been tested with Eclipse v15.6.

#The package is available on pypi . To install the package, just run:
pip install pydvh
# to upgrade from previous version, run:
pip install --upgrade pydvh
from pydvh import DVHType, DVHData, DVHFile
import matplotlib.pyplot as plt

def test_cdvh_plot(self):
        dvhfile = DVHFile.from_file_eclipse("P1_Summed_Plan_DDVH.txt")
        structures = dvhfile.structure_names()
        ddvh = dvhfile.get_dvh_by_name(structures[0])
        cdvh = ddvh.to_cumulative()
        # ddvh = cdvh.to_differential()

        dose_array, volume_array = zip(*cdvh.curve())
        plt.xlabel(f"Physical dose ({cdvh.dose_unit})")
        plt.ylabel(f"Volume ({cdvh.volume_unit})")
        plt.grid()
        plt.plot(dose_array, volume_array)
        plt.savefig("CDVH_bladder.png")
        plt.show()


#https://pyplati.github.io/platipy/_examples/dvh_analysis.html

#https://github.com/glucee/DVH

  
#https://simpleitk.org/