"""
Test of wavefronts
"""

import matplotlib.pyplot as plt
import optics.lens as lens
import optics.wavelength as w
import optics.ray as ray
import math
import optics.psf as p
import optics.analysis as a



def main():
    ln = lens.DataBaseLens("$LENS/Tessar-F4.5")
    pen = ray.RayPencil().addCollimatedBeam(ln,math.radians(0.0),"array",path=True)
    
    ip = ln.backFocalPlane()
    wp = ln.exitAperture()

    pen *= ln
    pen *= ip

    psf = p.Psf().optimalArea(pen,ip)
    print("PSF at : " + repr(psf))

    wf = a.WavePointSet().setWithPencil(pen,wp,psf)
    ze = wf.fitZernike()
    print(repr(ze))
    ze.draw()
    plt.show()
    

main()
