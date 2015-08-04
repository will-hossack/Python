import optics.lens as len
import optics.analysis as anal
import optics.surface as sur
import matplotlib.pyplot as plt
import vector
import optics.ray as ray
import math
import tio

def main():

    lens = len.DataBaseLens()
    #SimpleSinglet(0.0,100.0,12.5,"planoconvex")

    #cp = lens.cardinalPoints()
    #fl = lens.focalLength()

    #knife = sur.KnifeEdgeAperture(cp[1],10.0,-0.01)
    #output = sur.ImagePlane(cp[1] + vector.Vector3d(0,0,fl),30,30)

    #xsize  = 3.0*lens.entranceAperture().maxRadius

    angle = tio.getFloat("Angle in degrees")

    output = anal.knifeEdgeTest(lens,math.radians(angle),0.0)
    #u = vector.Unit3d(vector.Angle(math.radians(angle)))
    #output = anal.OpticalImage(cp[3].propagate(2*fl,u),xsize,xsize)

    #pencil = ray.RayPencil().addCollimatedBeam(lens,u,"array",50)
    #.addMonitor(ray.RayPath())

    #pencil *= lens
    #psf = anal.Psf().optimalArea(pencil,cp[1].z)
    #knife = sur.KnifeEdgeAperture(psf,10.0,0.0,math.radians(0))
    #pencil *= knife
    #pencil *= output

    #lens.draw()
    #knife.draw()
    #pencil.draw()
    output.draw()
    plt.show()

main()
    
