#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example program to simulate a simple prism spectrometer. The system is set to min
deviation at Mercury e line, then output angle are calcualted between Mercury i line
and Helium r line by ray tracing


"""

import optics.ray as r
import optics.lens as l
import optics.wavelength as w
import matplotlib.pyplot as plt
import vector as v
import math
import tio as t
import numpy as np

def main():
    
    #      Get the materail type and make a prism of default angle, size and location
    n = w.MaterialIndex()
    prism = l.Prism(index = n)
    t.tprint(repr(prism))
    
    
    #      Get input point on prism and min deviation at Mercury_i line
    pt = prism.getInputPoint()
    dev = prism.minDeviation(w.Mercury_e)
    t.tprint("Min deviation : ", math.degrees(dev), " at : ",w.Mercury_e)
    t.tprint("Max resolutions is ",prism.maxResolution(w.Mercury_e))
    t.tprint("Resolution with 20 mm diameter beam : ", prism.resolution(10,w.Mercury_e))
    
    u = v.Angle(dev/2)      # Set ray input angle at half min deviation
    
    #      Form np array of wavelength and angle
    wavelengths = np.linspace(w.Mercury_i,w.Helium_r,50)
    angle = np.zeros(wavelengths.size)
    
    #      Go through each wavelength, make a ray and trace it
    for i,wave in enumerate(wavelengths):
        ray = r.IntensityRay(pt,u,wave)
        ray *= prism
        #      Extract angle of ray in degrees
        angle[i] = math.degrees(v.Angle(ray.director).theta)
         
    # Do the plotting
    
    plt.plot(wavelengths,angle)
    plt.title("Spectrometer Output angle")
    plt.xlabel("Wavelength in microns")
    plt.ylabel("Output angle in degrees")
    plt.grid()
    

    plt.show()
    
    
main()