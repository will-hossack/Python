"""
 

    Example program to similate the output of of a prism spectrometer for a range
    of wavelengths read from a csv file

"""

from optics.lens import Prism
from optics.wavelength import MaterialIndex,Mercury_e,Sodium_D
import tio as t
from  vector import Angle 
from optics.ray import RayPencil
import math
import numpy as np
from scipy.special import j1
import matplotlib.pyplot as plt
from csvfile import readCSV

def jinc(x):
    if x == 0.0:
        return 1.0
    else:
        return 2*j1(x)/x


def main():
    
    #      First set up spctrometer
    
    
    n = MaterialIndex()
    prismAngle = t.getFloat("Prism angle in degrees",60.0)
    prismHeight = t.getFloat("Height of prism",100.0)
    prism = Prism(0.0,prismAngle,prismHeight,n)
    
    t.tprint(repr(prism))
    
    radius = t.getFloat("Beam Radius",10.0)
    setupWavelength = t.getFloat("Set up wavelength",Mercury_e)
    deviation = prism.minDeviation(setupWavelength)
    t.tprint("Min deviation : ", math.degrees(deviation), " at : ",setupWavelength)
    t.tprint("Max resolutions is ",prism.maxResolution(setupWavelength))
    t.tprint("Resolution with specified beam : ", prism.resolution(radius,setupWavelength))
    
    
    #          Get the wavelengths and intensities of the spcetrum
    fileName = t.getFilename("Wavelength file","csv")
    wavelengths,intensities = readCSV(fileName)
    
    #          Set to ray paramteers
    pt = prism.getInputPoint()         # Input point on prism
    u = Angle(deviation/2)             # Angle of input beam            
    pencil = RayPencil().addRays(pt,u,wavelengths,intensities) # Make pencil of rays
    pencil *= prism                    # Put through prism
    
    #          Extarct the information about the rays
    angle = np.zeros(len(pencil))
    peaks = np.zeros(len(pencil))
    widths = np.zeros(len(pencil))
    
    for i,ray in enumerate(pencil):
        a = ray.getAngle()
        angle[i] = a.theta*math.cos(a.psi)
        peaks[i] = ray.getIntensity()
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
