import optics.lens as len
import optics.ray as ray
import optics.wavelength as wl
import matplotlib.pyplot as plt
import math
from vector import Vector3d,Unit3d,Angle

def main():


    eye = len.Eye(0.0)
    print("Original Focalength : " + str(eye.backFocalLength(wl.PhotopicPeak)))
    eye.setIris(0.5)
    
    eye.setNearPoint(300)
    print("Modified Focalength : " + str(eye.backFocalLength(wl.PhotopicPeak)))
    
    theta = 10
    u = Unit3d(Angle(math.radians(theta)))
    rpencil = ray.RayPencil().addBeam(eye,u,"vl",10,wl.Red).\
             addMonitor(ray.RayPath())
    rpencil *= eye
    gpencil = ray.RayPencil().addBeam(eye,u,"vl",10,wl.Green).\
             addMonitor(ray.RayPath())
    gpencil *= eye
    bpencil = ray.RayPencil().addBeam(eye,u,"vl",10,wl.Blue).\
             addMonitor(ray.RayPath())
    bpencil *= eye

    eye.draw()
    #rpencil.draw()
    gpencil.draw()
    #bpencil.draw()
    plt.axis("equal")
    plt.grid()
    plt.show()
    
    

main()
