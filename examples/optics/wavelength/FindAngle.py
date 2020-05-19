
"""
Created on Sat May 16 09:28:14 2020



@author: wjh
"""

from optics.lens import Prism
from optics.ray import IntensityRay
from vector import Unit3d,Angle
import tio as t
import math

def getWavelength(prism, inAngle, outAngle, wavelengths = [0.25,1.0]):
    
    #        Get prism point and angle of input at Unit3d
    #
    pt = prism.getInputPoint()
    u = Unit3d(Angle(inAngle))
    
    #         Guess at initial wavelngth
    wave = (wavelengths[1] - wavelengths[0])/2
    #         Make input ray at guess wavelength
    ray = IntensityRay(pt,u,wave)
    
    #       Parameters for seaerch
    delta = 0.1
    forward = True
    na = float("inf")   # New angle
    
    while abs(na - outAngle) > 1.0e-9/abs(outAngle) :
        nray = ray*prism       #      New Ray through prism
        na = nray.getAngle()
        na = na.theta*math.cos(na.psi)    # In radians
        if na < outAngle:                       # Less that target
            wave += delta
            forward = True
        else:   
            if forward:                   # Half step
                delta *= 0.5
            forward = False
            wave -= delta
        if wave < wavelengths[0] or wave > wavelengths[1]:
            print("Out of wavelength range :")
            return float("nan")
        
        ray.wavelength = wave             # Update wavelength of ray
        
    return ray.getWavelength()                            # End of loop, so success


def main():
    
    prism = Prism()          # Default prism
    
    
    
    ia = math.radians(t.getFloat("Input angle"))
    oa = math.radians(t.getFloat("Output Angle"))
    
    w = prism.getWavelength(ia,oa)
    print("Wavelength from function " + str(w))
    
main()