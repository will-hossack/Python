"""
Test of wavefronts
"""

import tio as t
import matplotlib.pyplot as plt
import optics.lens as lens
import optics.wavelength as w
import optics.ray as ray
import math
import optics.psf as p
import optics.analysis as a



def main():
    ln = lens.DataBaseLens("$LENS/Singlet-F6.3")
    angle = math.radians(t.getFloat("Angle",0.0))
    pen = ray.RayPencil().addCollimatedBeam(ln,angle,"array",path=True)
    
    ip = ln.backFocalPlane()
    wp = ln.exitAperture()

    pen *= ln
    pen *= ip

    psf = p.Psf().setWithRays(pen,ip)
    print("PSF at : " + repr(psf))

    wf = a.WavePointSet().setWithPencil(pen,wp,psf)
    ze = wf.fitZernike()
    print(repr(ze))


    psf = ze.getPSF(log=True)
    plt.imshow(psf,cmap=plt.cm.gray,extent=(-1.0,1.0,-1.0,1.0))


    
    #ze.draw(xtilt = 2.0, ytilt = 2.0)
    plt.show()
    

main()
