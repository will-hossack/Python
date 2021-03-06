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
    se = wa.fitSeidel(angle,wave,0)
    t.tprint("Reference point is : ",wa.refpt)
    
    t.tprint(repr(se))

    #    Plot zernike as interfometer plot with tilt of 2 fringes


    plt.subplot(2,1,1)
    se.plotImage(xtilt = 2.0)
    plt.subplot(2,1,2)
    #se.plotOTF(128,"b")
    plt.show()
    

main()
