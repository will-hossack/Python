"""
   Test for materials 
"""

import optics.material as mat
import optics.wavelength as w
import matplotlib.pyplot as plt

def main() :

    index = mat.MaterialData().getIndex()

    index.draw(w.RefractiveIndexColour(index))
    plt.legend()
    plt.show()
    




main()
