"""
         Test of Paraxial Rays
"""
import tio
import optics.ray as ray
import optics.matrix as mat
import matplotlib.pyplot as plt

def main():

    lens = mat.DataBaseMatrix("$LENS/Tessar-100")
    mag = -2.0
    h = 5.0
    op,ip  = lens.planePair(h,mag)
    ps = mat.ParaxialSystem(op,lens,ip)

    rt = ray.ParaxialRay(2.0,0.02,op.inputPlane())
    pathtop = ray.RayPath(rt)
    rl = ray.ParaxialRay(2.0,-0.03,op.inputPlane())
    pathlower = ray.RayPath(rl)
    tio.tprint(repr(rt))
    rt *= ps
    rl *=ps
    tio.tprint(repr(rt))
    tio.tprint(pathtop.getInfo())

    ps.draw()
    pathtop.draw()
    pathlower.draw()
    plt.show()

main()
