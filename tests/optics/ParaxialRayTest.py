"""
  Test paraxial rays
"""
import optics.ray as r
import optics.matrix as m
import matplotlib.pyplot as plt
import math
import tio

def main():
    ray = r.ParaxialRay(0,math.radians(5),100)
    ray.addMonitor(r.RayPath())
    

    ray += 160
    tio.tprint(repr(ray))
    lens = m.ThickLensMatrix(0.01,1.5,20.0,-0.015)
    pg = m.ParaxialGroup(160,lens)
    tio.tprint(lens)
    tio.tprint("Focal length is : " + str(lens.backFocalLength()))
   
    ray *= lens
    tio.tprint(ray)
    ray += 160
    tio.tprint(ray)


    ray.draw()
    plt.show()
    
main()
    
