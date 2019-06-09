"""
         Test of Paraxial Rays
"""
import tio
import optics.ray as ray
import optics.matrix as mat
import matplotlib.pyplot as plt

def main():

    lens = mat.DataBaseMatrix("$LENS/Tessar-100")
    mag = -10.0
    h = 5.0
    op,ip  = lens.planePair(h,mag)
    ps = mat.ParaxialSystem(op,lens,ip)

    pencil = ray.RayPencil().addSourceParaxialBeam(lens,-5.0, op)
    pencil.addMonitor(ray.RayPath())

    pencil *= ps
    fig,ax = plt.subplots()
    ps.draw()
    pencil.draw()
    #ax.set_aspect(1.0)
    plt.show()

main()
