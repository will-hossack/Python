"""
   Example Programme to for a Spot Diagram


"""

import optics.lens as ln
import optics.psf as p
import vector as v
import optics.ray as r
from optics.wavelength import Default
from optics.analysis import SpotAnalysis
import tio as t
import matplotlib.pyplot as plt

def main():

    #      Get lens from database 
    lens = ln.DataBaseLens()       
    
    #           Get angle of beam and wavelnegth 
    angle = t.getFloat("Angle in degrees",0.0,0.0,15.0)
    u = v.Unit3d(v.Angle().setDegrees(angle))     # Angle as unit vectr
    w = t.getFloat("Wavelength",Default)

    #    Make two ray pencils, one for spot diagram and one for display (vertical only)
    vpencil = r.RayPencil().addBeam(lens,u,"vl",wave=w).addMonitor(r.RayPath())

    #            Propagate through lens to back focal plane
    bf = lens.backFocalPlane()
    vpencil *=lens
    vpencil *=bf

    #            Get optimal area psf and create a SpotDiagram 
    sd = SpotAnalysis(lens,u,0)

    #             Go round loop plotting the sopt diagram as various zplane positions

    while True:
        zp = t.getFloat("Delta",0.0)
        plt.subplot(2,1,1)
        lens.draw()
        vpencil.draw()
        plt.axis("equal")
        plt.title("Lens " + lens.title)

        plt.subplot(2,1,2)
        sd.draw(zp, True)
        pos = sd.plane.getPoint().z 
        print("Plane pos " + str(pos))
        plt.title("Spot diagram")
        plt.show(block = True)
    

main()
