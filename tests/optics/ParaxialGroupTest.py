"""
   Test of ParaxialGroup class
"""
import optics.matrix as m
import vector as v
import tio as t
import  matplotlib.pyplot as plt

def main():
    tl = m.ParaxialThickLens(20,0.015,1.74,15,-0.01,5)
    t.tprint(tl)

    t.tprint("Back focal length : ",tl.backFocalLength())

    ml = m.ParaxialMirror(100,-0.03,10)

    plt.figure(1)
    plt.subplot(2,1,1,axisbg="y")
    tl.draw(legend = True)

    plt.subplot(2,1,2)
    ml.draw()

    plt.show()

    

main()
