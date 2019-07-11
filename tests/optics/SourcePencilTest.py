"""
Test of source pencils
"""

import matplotlib.pyplot as plt
import optics.lens as lens
import optics.wavelength as w
import optics.ray as ray
import vector as v
import math

def main():
    ln = lens.DataBaseLens("$LENS/Tessar-F4.5")
    ln.setIris(0.5)
    obj,ima = ln.planePair(-0.2,300,300)

    onaxis = obj.getSourcePoint(0,0)
    offaxis = obj.getSourcePoint(0,100)
    pen = ray.RayPencil()
    pen.addSourceBeam(ln,onaxis)
    pen.addSourceBeam(ln,offaxis)
    pen.addMonitor(ray.RayPath())

    pen *= ln
    pen *= ima
    
    fig,ax = plt.subplots()
    obj.draw()
    ima.draw()
    pen.draw()
    ln.draw()
    ax.set_aspect(1.0)

    plt.show()


main()
