"""
   Test for materials 
"""

import optics.wavelength as w
import sys
import matplotlib.pyplot as plt

def main() :

    index = w.MaterialIndex("SF11")     # Get material
    if not index:
        print("No material")
    else:
        print(repr(index))
        index.draw(w.RefractiveIndexColour(index))
        plt.legend()
        plt.show()
    




main()
