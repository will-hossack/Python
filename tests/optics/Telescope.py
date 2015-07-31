import optics.lens as lens
import optics.ray as ray
import vector as v
import tio
import matplotlib.pyplot as plt


def main():

    front = lens.Singlet(0.0,7.7639e-3,2.6,0.0,12.7,"N-BK7")
    bfp = front.paraxialGroup().backFocalPlane()
    print("Back focal is at : " + str(bfp))

    back = lens.Singlet(0.0,1.95695e-2,2.6,-1.95695e-2,6.35,"N-BK7")
    ffp = back.paraxialGroup().frontFocalPlane()

    print("Front focal is at : " + str(ffp))

    back.setPoint(bfp - ffp)

    eye = lens.Eye(320.0)
    eye.setIris(0.7)

    pencil = ray.RayPencil().addCollimatedBeam(front,v.Unit3d(0,0,1)).addMonitor(ray.RayPath())

    pencil *= front
    pencil *= back
    pencil *= eye

    front.draw()
    back.draw()
    eye.draw()
    pencil.draw()
    plt.grid()
    plt.show()

main()
