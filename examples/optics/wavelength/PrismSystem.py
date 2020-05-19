"""
         Example program create a prism with three colout input beam 
         with the output beam straightened and put through a vieweing lens

"""
import math
from optics.lens import Prism
from optics.wavelength import Red,Green,Blue
from optics.ray import Disc,RayPencil,RayPath
from optics.lens import DataBaseLens
import matplotlib.pyplot as plt


def main():
    
    #      Make a 40m high prism of F4 glass at the origin
    prism = Prism(0,height = 40, index = "F4")
    #      Calcuate min deviation angle a 
    minDev = prism.minDeviation(Green)
    print("Min deviation at Green: {0:6.4f}".format(math.degrees(minDev)))
    
    # Form a pencil and add in three beams of R,G,B along the optical axis
    # Start of beeams given by 10mm diameter disc at -50mm before prism
    disc = Disc(-50.0 ,10.0)    
    pencil = RayPencil()               # Make blank pencil
    pencil.addBeam(disc,0.0,wave = Green)  # Add the three coloured beams.
    pencil.addBeam(disc,0.0,wave = Blue)   # default of 21 rays in vertical line
    pencil.addBeam(disc,0.0,wave = Red)
    angle = minDev/2                       # Angle of beam at half min deviation
    pencil.rotateAboutX(angle,prism.getInputPoint()) # retate beam about inpout point
    print("Numnber of rays : " + str(len(pencil)))
    pencil.addMonitor(RayPath())           # Add a monitor to allow plotting
    
    pencil *= prism                        # propagate beam through prism
    
    # Beam will be a dirction -angle with green rays centered at output point of prsim
    # so rotate beam so green pencil is parallel to optical axis
    pencil.rotateAboutX(angle,prism.getOutputPoint())
    
    #       Get a imagin lens, (a 140mm optimised doublet)
    lens = DataBaseLens("$LENS/Linos140Doublet")
    lens.setPoint(50.0)        # Set it 50mm beyond the centre of the prism
    bf = lens.backFocalPlane(Green)  # Get back focal plane where rays wil be images
    
    #      propagate pencil theorugh lens to back focal plane
    pencil *=lens
    pencil *= bf
    
    #       Generate the diagram
    prism.draw()
    lens.draw(True,True)
    bf.draw()
    pencil.draw()
    plt.axis("equal")
    plt.grid()
    plt.xlabel("mm along optical axis")
    plt.show()
    
    
    
main()
