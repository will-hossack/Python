"""
Example programme to read in lens, trace 2d ray pencil and calcualte optical
Zernike expansion of the wavefont in the exit pupil wrt to geometric PSF in the
paraxial image plane.

"""

import tio as t
import matplotlib.pyplot as plt
import optics.lens as l
import math
import optics.wavefront as a



def main():
    
    #        Get lens from database, get iris ration and set it
    lens = l.DataBaseLens()
    ratio = t.getFloat("Aperture Ratio",1.0)
    lens.setIris(ratio)
    
    #         Get angle and wavelength
    angle = math.radians(t.getFloat("Angle",0.0))
    wave = t.getFloat("Wavelength",0.55)
    design = 0.55         # Hard code design wavelength
    
    # set up wavefront analysis
    wa = a.WaveFrontAnalysis(lens,design)
    
    # do a 4th order Zernike fit with Collated lens 
    ze = wa.fitZernike(angle,wave,4,0)
    t.tprint("Reference point is : ",wa.refpt)
    
    t.tprint(repr(ze))

    #    Plot zernike as interfometer plot with tilt of 2 fringes
    inter = a.Interferometer(ze)
    inter.draw()
    plt.show()

main()
