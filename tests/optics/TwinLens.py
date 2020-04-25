import matplotlib.pyplot as plt
import optics.lens as l
import optics.ray as r
import optics.surface as s
import optics.matrix as m
import tio


def main():

    front = l.SimpleSinglet(0,100,10)
    back = l.SimpleSinglet(15,-300,10)
   
    system = l.OpticalSystem("Twin",front,back)
    
    # system = l.SimpleSinglet(0,80)


    print("Fl : " + str(system.backFocalLength()))
    obj,im = system.planePair(-0.2,10.0)
    source = obj.getSourcePoint(0.0,3.0)
    pencil = r.RayPencil().addBeam(system,source).addMonitor(r.RayPath())
    pencil *= system
    pencil *= im
    
    system.draw(False)
    system.paraxialGroup().draw(True)
    obj.draw()
    im.draw()
    pencil.draw()
    plt.show()


main()
    
