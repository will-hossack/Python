"""
   Set of test methods for the Paraxial Matrix class
"""

import optics.matrix as m
import tio as t

def main():


    
    
    l = m.ParaxialThinLens(10,0.01,1.6,-0.015,10)
    t.tprint(l)
    t.tprint("Front Focal Plane : " + str(l.frontFocalPlane()))
    t.tprint("Back focal length : " + str(l.backFocalLength()))
    t.tprint("Front focal length : " + str(l.frontFocalLength()))



    

main()
