"""
 

    Example program to similate the output of of a prism spectrometer for a range
    of wavelengths read from a csv file
"""
import math
from optics.lens import Prism
from optics.wavelength import MaterialIndex,Mercury_e
from optics.ray import RayPencil
import numpy as np
from scipy.special import j1
import tio as t
import matplotlib.pyplot as plt
from csvfile import readCSV

def jinc(x):
    if x == 0.0:
        return 1.0
    else:
        return 2*j1(x)/x


def main():
    
    #      First set up prism
    n = MaterialIndex()                   # get materail, angle and height
    prismAngle = t.getFloat("Prism angle in degrees",60.0)
    prismHeight = t.getFloat("Height of prism",100.0)
    prism = Prism(0.0,prismAngle,prismHeight,n)   # Make prism at the origon
    t.tprint(repr(prism))                         # print descrciption
    
    #     Get the bean and alignment wavelength.
    radius = t.getFloat("Beam Radius",10.0)
    setupWavelength = t.getFloat("Set up wavelength",Mercury_e)
    deviation = prism.minDeviation(setupWavelength)
    t.tprint("Min deviation : {0:6.4f} at : {1:6.4f}".format(math.degrees(deviation),setupWavelength))
    t.tprint("Max resolutions is : {0:7.3f}".format(prism.maxResolution(setupWavelength)))
    t.tprint("Resolution with specified beam : {0:7.3f}".format(prism.resolution(radius,setupWavelength)))
    
    #          Get the wavelengths and intensities of the spcetrum
    fileName = t.getFilename("Wavelength file","csv")
    wavelengths,intensities = readCSV(fileName)
    
    #          Set up ray paramteers
    pt = prism.getInputPoint()         # Input point on prism
    ang = deviation/2                 # Angle of imput beam
    pencil = RayPencil().addRays(pt,ang,wavelengths,intensities) # Make pencil of rays
    pencil *= prism                    # Put through prism
    
    #          Extarct the information about the rays into numpy arrays.
    angle = np.zeros(len(pencil))
    peaks = np.zeros(len(pencil))
    widths = np.zeros(len(pencil))
    
    for i,ray in enumerate(pencil):
        a = ray.getAngle()                    # Get the angle in theta,psi format
        angle[i] = a.theta*math.cos(a.psi)    # Extract -ve angle 
        peaks[i] = ray.getIntensity()
        #       Note width is lambda/(pi*diameter) but 1,000 to convert nn -> um
        widths[i] = ray.getWavelength()/(2000*math.pi*radius)

    #       Set up the plot paramters
    minField = np.min(angle) - 10*np.max(widths)
    maxField = np.max(angle) + 10*np.max(widths)
    npoints = max(int((maxField - minField)/np.min(widths)),400)
    
    #            Plot out the results
    fieldAngle = np.linspace(minField,maxField,npoints)
    spectralOutput = np.zeros(fieldAngle.size)
    
    for a,p,w in zip(angle,peaks,widths):
    
        for i,af in enumerate(fieldAngle):
            da = af - a
            s = jinc(da/w)
            spectralOutput[i] += p*s**2
        
    #     Do the actual plot        
    plt.plot(np.degrees(fieldAngle),spectralOutput)
    plt.grid()
    plt.title("Spectral plot")
    plt.xlabel("Angle in degrees")
    plt.ylabel("Intensty")
    plt.show()



    
    
main()
