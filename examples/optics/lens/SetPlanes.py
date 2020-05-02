"""
Example code to take a use a DataBaseLens and form an optical system with 
object plane at -300mm and image plane at +60mm. The lens is positioned 
at the paxaxial image position. 

The full ray traced pencil is then produced at a wavelnegth of 0.65 (red)

"""

import optics.lens as l
import optics.ray as ray
import optics.surface as sur
import matplotlib.pyplot as plt



def main():

    lens = l.DataBaseLens("$LENS/Tessar-F4.5") # Read in lens
    lens.setFocalLength(50)                    # set focal length by scaling
    print(repr(lens))


    #    Form the object and image planes at -300 and + 60 respectively
    obj = sur.ImagePlane(-300,100)
    ima = sur.ImagePlane(60)

    mag = lens.setWithPlanes(obj,ima)          # Set lens at paraxial locatiom
    print("Magnifications is : " + str(mag))
    print("Object Plane : " + str(repr(obj)))
    print("New Image Plane : " + str(repr(ima)))
    print(repr(lens))
    
    #     Form a ray pencil from a point in the object plane at (0,10) mm
    pencil = ray.RayPencil().addBeam(lens,obj.getSourcePoint(0.0,10.0),wave=0.65)
    pencil.addMonitor(ray.RayPath())      # Add a monitor to record positions

    #      Propagate through lens and to the image plane
    pencil *= lens
    pencil *= ima

    
    #     Make the drawing of lens, planes and pencil
    lens.draw(True,True)   # Plot with paraxial planes and legend
    obj.draw()
    ima.draw()
    pencil.draw()
    plt.title(lens.title)
    plt.axis("equal")
    plt.show()
main()
