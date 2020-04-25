import matplotlib.pyplot as plt
import optics.lens as l
import optics.ray as r
import optics.surface as s
import tio


def main():

    objective = l.DataBaseLens("$LENS/Tessar-F6.3")
    objective.setFocalLength(300.0)
    block = s.CircularAperture(150,20)
   

    eyelens = l.SimpleSinglet(330,30,10)

    system = l.OpticalSystem("Telescope",objective,block,eyelens)
    
    system.setIris(0.5)
    system.movePoint(-25.0,1)
    
    pencil = r.RayPencil().addBeam(system,0.0).addMonitor(r.RayPath())
    pencil *= system

    fm = system[0].paraxialGroup()
    
    system.draw(False)
    fm.draw()
    pencil.draw()
    plt.show()


main()
    
