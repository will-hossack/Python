"""
Example programme to read in lens, trace 2d ray pencil and calcualte optical
Zernike expansion of the wavefont in the exit pupil wrt to geometric PSF in the
paraxial image plane.

"""

import tio as t
import matplotlib.pyplot as plt
import optics.lens as l
import optics.ray as ray
import math
import optics.psf as p
import optics.wavefront as a



def main():
    lens = l.DataBaseLens()
    angle = math.radians(t.getFloat("Angle",0.0))
    ratio = t.getFloat("Aperture Ratio",1.0)
    lens.setIris(ratio)
    wave = t.getFloat("Wavelength",0.55)
    design = 0.55         # Hard code design wavelength
    
    wa = a.WaveFrontAnalysis(lens,design)
    ze = wa.fitZernike(angle,wave,4,0)
    t.tprint("Reference point is : ",wa.refpt)
    
    t.tprint(repr(ze))

    #    Plot zernike as interfometer plot with tilt of 2 fringes
    ze.plotImage(ytilt = 2.0)
    plt.show()
    

main()
