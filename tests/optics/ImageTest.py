import optics.lens as len
import optics.analysis as an
import optics.ray as ray
import matplotlib.pyplot as plt
import tio
import time

def main():

    lens = len.DataBaseLens()
    print("Focal length is " + str(lens.focalLength()))

    mag = tio.getFloat("Magnification",-0.1)
    osize = tio.getFloat("Input plane size in mm",200)

    obj = an.OpticalImage(0.0,osize,osize,200,200).addTestGrid()

    print(repr(obj))

    startTime = time.clock()          # Start timer

    image = obj.getImage(lens,mag)


    print(repr(image))

    deltaTime = time.clock() - startTime # fine CPU time
    
    print("Time was  {0:6.2f} CPU seconds".format(deltaTime))

    image.draw()
    plt.grid()
    plt.show()
    
main()
