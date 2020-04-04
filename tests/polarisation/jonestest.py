import optics.jones as j
import matplotlib.pyplot as plt
import tio
import math

def main():

    theta = math.radians(tio.getFloat("Angle"))
    
    beam = j.RightCircularPolarisedBeam()
    print(repr(beam))

    polar = j.LinearPolariser(theta)

    beam *= polar
    print(repr(beam))
    

    beam.polarPlot()
    
    plt.show()

main()
