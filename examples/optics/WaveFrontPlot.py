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


    wave = a.WaveFront().readFromFile()
    wave.plotImage(ytilt = 2.0)
    
    plt.show()
    

main()
