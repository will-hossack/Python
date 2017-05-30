"""
  Test paraxial rays
"""
import optics.ray as r
import optics.matrix as m
import matplotlib.pyplot as plt
import math
import tio

def main():
    

    lens = m.ParaxialThickLens(50,0.01,1.74,10.0,-0.002,10)
    pencil = r.RayPencil().addCollimatedParaxialBeam(lens,math.radians(5.0))
    pencil.addMonitor(r.RayPath())

    pencil *= lens
    focal = lens.backFocalPlane()
    pencil *= focal

    lens.draw(legend=True)
    pencil.draw()
    plt.show()
    
main()
    
