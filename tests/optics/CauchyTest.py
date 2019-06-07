"""
   Test for materials 
"""

import optics.wavelength as w
import sys
import matplotlib.pyplot as plt

def main() :

    index = w.CauchyIndex(564745)
    if not index:
        print("No material")
    else:
        print(repr(index))
        index.draw(w.RefractiveIndexColour(index))
        plt.legend()
        plt.title(repr(index))
        plt.show()
    




main()
