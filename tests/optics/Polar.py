"""          Test for polarsiation and Jones classes
"""
import math
import optics.jones as j
import tio as t

def main():
    a = j.LeftCircularPolarisedBeam(3.0)
    t.tprint(repr(a))
    t.tprint("Intensity is ","{0:5.3f}".format(a.getIntensity()))
    t.tprint("Phase is ","{0:5.3f}".format(a.getPhase()))
    


main()
