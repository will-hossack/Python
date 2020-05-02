"""
Example code to take a single thick lens, scale ot 50 mm focal length
and plot out it optimal location for a range on ojbect locations.

"""
import optics.matrix as mat
import matplotlib.pyplot as plt
import numpy as np

def main():

    lens = mat.ParaxialThickLens(30,0.025,1.61,6.0,-0.035,5.0)
    lens.setFocalLength(50)
    f = lens.backFocalLength()
    print("New focal length : " + str(f))
    print(repr(lens))
    
    # for distance in log space for object distance
    distance = np.logspace(2.5,5,30)
    zpos = np.zeros(distance.size)
    
    #     Fix image plane at 60 mm
    ima = mat.ParaxialPlane(60,10) 

    for i,d in enumerate(distance):
        obj = mat.ParaxialPlane(-d,50)   # Set object plane
        lens.setWithPlanes(obj,ima)      # Set lens position
        zpos[i] = lens.inputPlane()      # record lens position

    # Plot out grap with distance on a log scale
    plt.plot(distance,zpos)
    plt.xscale("log")
    plt.show()

    
main()
