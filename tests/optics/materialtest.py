"""
   Test for materials 
"""

import optics.wavelength as w
import sys
import matplotlib.pyplot as plt

def main() :

    index = w.MaterialIndex("SF11")
    print(str(index.valid))
    if not index:
        print("No material")
    else:

        index.draw(w.RefractiveIndexColour(index))
        plt.legend()
        plt.show()
    




main()
