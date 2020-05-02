"""
Example code to take a system of two thick singlet lenses the first on art 50 mm
and plot the  focal length of as the second lens is moved from 60 to 100 mm.

"""
import optics.lens as l
import optics.ray as ray
import matplotlib.pyplot as plt
import numpy as np
import math

def main():

    lens = mat.ParaxialThickLens(30,0.025,1.61,6.0,-0.035,5.0)
    lens.setFocalLength(50)
    f = lens.backFocalLength()
    print("New focal length : " + str(f))
    print(repr(lens))

    distance = np.logspace(2.5,5,30)
    zpos = np.zeros(distance.size)
    ima = mat.ParaxialPlane(60,10)

    for i,d in enumerate(distance):
        obj = mat.ParaxialPlane(-d,50)
        mag = lens.setWithPlanes(obj,ima)
        zpos[i] = lens.inputPlane()


    plt.plot(distance,zpos)
    plt.xscale("log")
    plt.show()

    
main()
