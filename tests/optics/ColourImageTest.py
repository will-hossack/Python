"""      Tests of Optical Image
"""

import matplotlib.pyplot as plt
import optics.analysis as ana
import optics.surface as s
import optics.lens as l

def main():

    oi = ana.ColourImage(0,300,300,300,300)
    oi.addTestGrid(40,40,intensity = [1,1,1])

    oi.draw()
    plt.show()
    lens = l.SimpleSinglet(200,80,10,index="SF11")
    #l.DataBaseLens("$LENS/Tessar-F4.5")
    #lens.setFocalLength(80)

    #obj,ima = lens.planePair(-0.2,oi.xsize,oi.ysize)

    #print("Object plane : " + repr(obj))
    #print("Image plane : "+ repr(ima))

    #oi.setPoint(obj.point)


    im = oi.getSystemImage(lens,-0.2)
    
    print(repr(oi))
    print(repr(im))

    #im = oi.getImage(lens,ima)

    im.draw()

    plt.show()


main()
