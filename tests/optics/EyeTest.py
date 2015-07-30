import optics.lens as len
import optics.ray as ray
import optics.wavelength as wl
import matplotlib.pyplot as plt
import math
from vector import Vector3d

def main():


    eye = len.Eye(0.0)
    print("Original Focalength : " + str(eye.focalLength()))
    #eye.draw()
    eye.setIris(1.0)
    eye.setNearPoint(1000.0)


    #print(repr(eye[4]))

    
    theta = 0.0
    u = ray.Director(ray.Angle(math.radians(theta)))
    pt = eye.frontNodalPoint() - Vector3d(0,20,1000)
    rpencil = ray.RayPencil().addSourceBeam(eye,pt,"vl",10,wl.Red).\
             addMonitor(ray.RayPath())
    rpencil *= eye
    gpencil = ray.RayPencil().addSourceBeam(eye,pt,"vl",10,wl.Green).\
             addMonitor(ray.RayPath())
    gpencil *= eye
    bpencil = ray.RayPencil().addSourceBeam(eye,pt,"vl",10,wl.Blue).\
             addMonitor(ray.RayPath())
    bpencil *= eye

    eye.draw()
    rpencil.draw()
    gpencil.draw()
    bpencil.draw()
    plt.grid()
    plt.show()

main()
