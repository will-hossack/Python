"""
  Test paraxial rays
"""
import optics.ray as r
import optics.matrix as m
import matplotlib.pyplot as plt
import math
import tio

def main():
    ray = r.ParaxialRay(0,math.radians(5),10)
    ray.addMonitor(r.RayPath())
    
    #tio.tprint(repr(ray))
    lens = m.ParaxialThickLens(160,0.01,1.7,20.0,-0.015,radius=20)
    #tio.tprint(lens)
    tio.tprint("Focal length is : " + str(lens.backFocalLength()))
   
    ray *= lens
    #tio.tprint(ray)
    ray += 160
    #tio.tprint(ray)

    lens.draw(legend=True)
    ray.draw()
    plt.show()
    
main()
    
