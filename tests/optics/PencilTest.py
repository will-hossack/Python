"""
Test of pencils
"""

import matplotlib.pyplot as plt
import optics.lens as lens
import optics.wavelength as w
import optics.ray as ray
import math

def main():
    ln = lens.DataBaseLens("$LENS/Tessar-F4.5")
    pen = ray.RayPencil()
    sp = w.PlanckSpectrum(5000)
    pen.addCollimatedBeam(ln,math.radians(10.0),wave=0.6,intensity=sp)
    pen.addMonitor(ray.RayPath())
    ip = ln.backFocalPlane(0.6)

    pen *= ln
    pen *= ip


    fig,ax = plt.subplots()
    ln.draw()
    pen.draw()
    ax.set_aspect(1.0)

    plt.show()


main()
