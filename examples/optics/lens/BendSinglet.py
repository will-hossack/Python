"""
   Programme to explore the bending of a single lens

Author: Will Hossack, The University of Edinburgh
"""

import optics.lens as ln
import tio as t
import matplotlib.pyplot as plt
import numpy as np

def main():
    #      Get the lens paramters
    c = t.getFloat("Curvature",1.0e-2)
    th = t.getFloat("Thickness",5.0)
    r = t.getFloat("Radius",10.0,0.0)
    
    #      Form the lens bfrom BK7 glass
    lens = ln.Singlet(0.0,c,th,-c,r,"BK7")
    t.tprint("Focal length : ",lens.backFocalLength()," bend is : ",lens.getBend())

    #        Create np arrays to hold the data
    bendData = np.linspace(-1.0,1.0,50)
    focalData = np.zeros(bendData.size)
    
    
    for i,b in enumerate(bendData):    # For each bend
        lens.setBend(b,False)          # Bend the lens, but do not rests focal length
        focalData[i] = lens.backFocalLength()  # Record the focal length

    #      Plot the graph
    plt.plot(bendData,focalData)
    plt.title("Focal length against bend for a thick singlet")
    plt.xlabel("Bend")
    plt.ylabel("Focal Length")
    plt.show()

    

main()
