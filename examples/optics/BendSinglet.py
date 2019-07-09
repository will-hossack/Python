"""
   Programme to explore the bending of a single lens

Author: Will Hossack, The University of Edinburgh
"""

import optics.lens as ln
import tio as t
import matplotlib.pyplot as plt

def main():
    c = t.getFloat("Curvature",1.0e-2)
    th = t.getFloat("Thickness",5.0)
    r = t.getFloat("Radius",10.0,0.0)

    lens = ln.Singlet(0.0,c,th,-c,r)
    t.tprint("Focal length : ",lens.backFocalLength()," bend is : ",lens.getBend())

    bData = []
    fData = []
    
    
    for i in range(-100,101):
        b =i/100.0
        bData.append(b)
        lens.setBend(b,False)
        fData.append(lens.backFocalLength())

    plt.plot(bData,fData)
    plt.show()

    

main()
