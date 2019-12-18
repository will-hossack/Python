"""
Simple lens explorer program
"""

from optics.lens import DataBaseLens
from optics.wavelength import Default
import tio as t
import sys

opts = ("exit","wavelength","iris")


def main():

    lens = DataBaseLens()
    t.tprint(repr(lens))
    design = Default
    wavelength = Default

    
    while True:
        opt,nopt = t.getOption("Opt",opts)
        if opt == 0:
            sys.exit(0)
        elif opt == 1:
            wavelength = t.getFloat("Wavelength",Default)
            t.tprint("Wavelength set to : " + str(wavelength))
        elif opt == 2 :
            iris = t.getFloat("Iris ratio",1.0,0.0,1.0)
            lens.setIris(iris)
            t.tprint("Iris set to : " + str(iris))
    
        else:
            t.tprint("Unknown option")

main()
