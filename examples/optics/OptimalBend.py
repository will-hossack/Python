"""
   Programme to explore the bending of a single lens

Author: Will Hossack, The University of Edinburgh
"""

import optics.lens as ln
import optics.analysis as a
import vector as v
import optics.ray as r
import tio as t
import matplotlib.pyplot as plt

def main():
    c = t.getFloat("Curvature",1.0e-2)
    th = t.getFloat("Thickness",5.0)
    radius = t.getFloat("Radius",10.0,0.0)

    lens = ln.Singlet(0.0,c,th,-c,radius)
    t.tprint("Focal length : ",lens.focalLength()," bend is : ",lens.getBend())
    
    #           Get angle of beam
    angle = t.getFloat("Angle in degrees",0.0,0.0,15.0)
    u = v.Unit3d(v.Angle().setDegrees(angle))

    bData = []
    aData = []

    for i in range(-100,101):
        b = i/100.0
        bData.append(b)

        lens.setBend(b)
    
        #    Make a ray pencil
        pencil = r.RayPencil().addCollimatedBeam(lens,u,"array")
        bf = lens.backFocalPlane()

        #            Propagate through lens to back focal plane
        pencil *= lens
        pencil *= bf

        psf = a.Psf().optimalArea(pencil,bf)

        aData.append(psf.area())

        
    plt.plot(bData,aData)
    plt.show()
    

main()
