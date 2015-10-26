"""
   Test for materials 
"""

import optics.material as mat
import matplotlib.pyplot as plt

def main() :

    index = mat.MaterialData().getIndex()

    index.draw()
    plt.legend()
    plt.show()
    




main()
