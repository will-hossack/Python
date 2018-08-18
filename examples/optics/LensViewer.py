"""
View an display lens in MatPlotLib with a typical ray pencil.

Author: Will Hossack, The University of Edinburgh
"""
import optics.lens as len
import optics.ray as ray
import vector as v
import tio as t
import matplotlib.pyplot as plt
    

def main():
    
    #
    #      Read lens in from database
    #
    lens = len.DataBaseLens()
    lens.setIris(0.7)           # Set iris to 0.7 of max
    #
    #       Make default pencil and add ray monitor to each ray
    angle = 2.0
    u = v.Unit3d(v.Angle().setDegrees(angle))       # Unit3d of ray pencil.
    pencil = ray.RayPencil().addCollimatedBeam(lens,u,"vl").addMonitor(ray.RayPath())
    #
    t.tprint("Focal length is : ",lens.focalLength())
    t.tprint("Petzal sum is : ",lens.petzvalSum())
    #
    #        Set the output plane (being the back focal plane)
    op = lens.backFocalPlane()

    #         Propagate pencil through lens and one to back plane
    pencil *= lens       # Through lens
    pencil *= op         # To plane
    #
    #                    Draw the diagram
    fig = plt.figure()
    panel = fig.add_subplot(1,1,1)
    panel.axis('equal')
    lens.draw()
    op.draw()
    pencil.draw()
    #                    Add decorations.
    plt.grid()
    plt.xlabel("Optical Axis")
    plt.ylabel("Height")
    plt.title("Diagram of lens " + lens.title)
    plt.show()

main()




    





