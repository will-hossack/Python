"""
   Programme to for a Spot Diagram

Author: Will Hossack, The University of Edinburgh
"""

import optics.lens as ln
import optics.psf as p
import vector as v
import optics.ray as r
from optics.surface import OpticalPlane
from optics.wavelength import Default
import tio as t
import matplotlib.pyplot as plt

def main():

    lens = ln.DataBaseLens()
    
    #           Get angle of beam
    angle = t.getFloat("Angle in degrees",0.0,0.0,15.0)
    w = t.getFloat("Wavelength",Default)
    u = v.Unit3d(v.Angle().setDegrees(angle))

    
    #    Make two ray pencils, one for spot diagram and one for display (vertical only)
    pencil = r.RayPencil().addCollimatedBeam(lens,u,"array",wave=w)
    vpencil = r.RayPencil().addCollimatedBeam(lens,u,"vl",wave=w).addMonitor(r.RayPath())
    bf = lens.backFocalPlane()

    #            Propagate through lens to back focal plane
    pencil *= lens
    pencil *= bf
    vpencil *=lens
    vpencil *=bf

    #            Get optimal area psf and create a SpotDiagram 
    psf = p.Psf().optimalArea(pencil,bf)
    sd = p.SpotDiagram(pencil)

    #             Go round loop plotting the sopt diagram as various zplane positions

    while True:
        zp = t.getFloat("Zplane",psf.z)
        plane = OpticalPlane(zp)
        plt.subplot(2,1,1)
        lens.draw()
        vpencil.draw()
        plt.axis("equal")
        plt.title("Lens " + lens.title)

        plt.subplot(2,1,2)
        sd.draw(plane)
        plt.title("Spot diagram")
        plt.show(block = True)
    

main()
