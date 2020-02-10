import tio as t
import optics.analysis as anal
import matplotlib.pyplot as plt
import optics.lens as l
import vector as v
import optics.ray as r
import optics.psf as p


def main():


    lens = l.SimpleSinglet(200,80,10)
    obj,im = lens.planePair(-0.2,200.0)
    print("Object plane : " + repr(obj))
    print("Image plane : " + repr(im))
    

    objectTarget = anal.TargetPlane(obj,wave=0.6)
    imageTarget = anal.TargetPlane(im)

    objectTarget.addGrid(5,5)

    fig = plt.figure()
    panel = fig.add_subplot(1,1,1)
    panel.axis('equal')

    imageTarget.draw()

    targetPencil = objectTarget.getPencils(lens, wave = 0.45)

    for pencil in targetPencil:
        pencil *= lens
        pencil *=imageTarget
        psf = p.Psf().setWithRays(pencil,imageTarget)
        psf.draw()

    plt.show()
    
main()
