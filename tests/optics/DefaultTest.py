
"""
    Test of the default /degign wavelength
"""

from optics.wavelength import Default,getDefaultWavelength,setDefaultWavelength
from optics.ray import RayPencil,Disc

def getW(w = None):
    if w == None:
        w = getDefaultWavelength()
        
    return w

def main():
    
    w = getDefaultWavelength()
    print("Current is : " + str(w))
    print("Method access" + str(getW()))

    setDefaultWavelength(0.7)

    w = getDefaultWavelength()
    print("New is : " + str(w))
    print("Method access" + str(getW()))
    
    
    disc = Disc(100,10)
    pencil = RayPencil().addBeam(disc,0.0,wavelength = 0.65, intensity = 10)
    #   print("Ray is " + repr(pencil))
    
main()