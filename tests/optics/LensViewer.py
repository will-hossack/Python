"""
View an display lens in MatPlotLib
"""
import optics.lens as len
import optics.ray as ray
import optics.surface as sur
import vector as v
import tio
import matplotlib.pyplot as plt
import math

def main():

    #
    #      Read lens in from database
    #
    lens = len.DataBaseLens()
    angle = 2.0
    lens.setIris(0.7)
    #
    #       Make default pancil and add ray monitor to each ray
    u = v.Unit3d(v.Angle(math.radians(angle)))
    pencil = ray.RayPencil().addCollimatedBeam(lens,u,"vl").addMonitor(ray.RayPath())
    #
    print("Focal length is : {0:7.5f}".format(lens.focalLength()))
    print("Petzal sum is : {0:7.5f}".format(lens.petzvalSum()))
    #
    #        Set the output plane
    op = sur.ImagePlane(lens.imagePoint(u).z)   

    pencil *= lens       # Propagare through lens
    pencil *= op         # Then to output plane
    #
    #                    Draw the diagram
    lens.draw()
    op.draw()
    pencil.draw()
    plt.grid()
    plt.xlabel("Optical Axis")
    plt.ylabel("Height")
    plt.title("Diagram of lens " + lens.title)
    plt.show()

main()




    





