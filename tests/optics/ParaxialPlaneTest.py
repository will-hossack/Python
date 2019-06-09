"""
     Test of Paraxial Plane 
"""
import optics.matrix as m
import matplotlib.pyplot as plt
import tio as t

def main():

    lens = m.DataBaseMatrix("$LENS/Tessar-100")
    lens.setFocalLength(50)
    t.tprint(lens.getInfo())


    mag = -2.0
    h = 5.0
    op,ip  = lens.planePair(h,mag)

    t.tprint(op.getInfo())
    t.tprint(ip.getInfo())

    ps = m.ParaxialSystem(op,lens,ip)

    fig,ax = plt.subplots()
    ps.draw()
    ax.set_aspect(1.0)
    
    plt.show()


main()

                        
