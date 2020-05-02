"""
   Programme to explore the bending of a single lens

Author: Will Hossack, The University of Edinburgh
"""

import optics.lens as ln
import optics.psf as p
import vector as v
import optics.ray as r
import tio as t
import matplotlib.pyplot as plt
import numpy as np

def main():
    c = t.getFloat("Curvature",1.0e-2)
    th = t.getFloat("Thickness",5.0)
    radius = t.getFloat("Radius",10.0,0.0)

    lens = ln.Singlet(0.0,c,th,-c,radius,"BK7")
    t.tprint("Focal length : ",lens.backFocalLength()," bend is : ",lens.getBend())
    
    #           Get angle of beam
    angle = t.getFloat("Angle in degrees",0.0,0.0,15.0)
    u = v.Unit3d(v.Angle().setDegrees(angle))

    bendData = np.linspace(-1.0,2.0,50)
    areaData = np.zeros(bendData.size)
    minarea = float("inf")
    
    for i,b in enumerate(bendData):
        lens.setBend(b)             # Bend lens and retain focal length
    
        #    Make a ray pencil
        pencil = r.RayPencil().addBeam(lens,u,"array")
        bf = lens.backFocalPlane()

        #            Propagate through lens to back focal plane
        pencil *= lens
        pencil *= bf
        
        # Form a psf in the back focal plane
        psf = p.Psf().setWithRays(pencil,bf)
        # Record the area of the psf
        area = psf.area()
        areaData[i] = area
        if area < minarea:
            minarea = area
            minbend = b
        
        
    t.tprint("Minimin area as bend : ","{0:5.3f}".format(minbend))
    plt.plot(bendData,areaData)
    plt.title("Plot of Psf area against bend for a singlet")
    plt.xlabel("Bend")
    plt.ylabel("Psf area")
    plt.show()
    

main()
