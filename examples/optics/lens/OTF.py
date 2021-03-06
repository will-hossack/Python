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
    lens = l.DataBaseLens()
    angle = math.radians(t.getFloat("Angle",0.0))
    ratio = t.getFloat("Aperture Ratio",1.0)
    lens.setIris(ratio)
    wave = t.getFloat("Wavelength",0.55)
    design = 0.55         # Hard code design wavelength
    
    wa = a.WaveFrontAnalysis(lens,design)
    ze = wa.fitZernike(angle,wave,4,1)
    t.tprint("Reference point is : ",wa.refpt)
    
    t.tprint(repr(ze))

    #    Plot zernike as interfometer plot with tilt of 2 fringes


    plt.subplot(2,1,1)
    inter = a.Interferometer(ze,xtilt=2.0)
    inter.draw()
    plt.subplot(2,1,2)
    ze.plotOTF(128,"b")
    plt.show()
    

main()
