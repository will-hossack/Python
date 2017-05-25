"""
   Set of test methods for the Paraxial Matrix class
"""

import optics.matrix as m
import tio as t

def main():

    first = t.getFloat("First lens",80.0)
    separ = t.getFloat("Seperation",30.0)
    second = t.getFloat("Second lens",50.0)
    mat = m.ThinLensMatrix(first)
    mat += separ
    mat *= m.ThinLensMatrix(second)
    t.tprint("Focal length of system is : " + str(mat.backFocalLength()))
    
    p = 1.0/first + 1.0/second - separ/(first*second)
    f = 1.0/p
    t.tprint("Direct calcualtion gives : " + str(f))
    

main()
