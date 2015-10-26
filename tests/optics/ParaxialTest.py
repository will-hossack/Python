""""
    Test of Paraxial matrix classes
"""

import optics.matrix as mat
import matplotlib.pyplot as plt

def main():

    obj = mat.ParaxialGroup(10,mat.ThickLensMatrix(0.002,1.5,10,-0.015),20)
    print(repr(obj))
    eye = mat.ParaxialGroup(200,mat.ThinLensMatrix(20),5)

    system = mat.ParaxialSystem(obj,eye)

    system.draw()
    


    plt.show()

main()
