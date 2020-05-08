#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 11:22:57 2020

@author: wjh
"""

import optics.lens as len
import optics.ray as ray
import optics.wavelength as wl
import matplotlib.pyplot as plt
import math
from vector import Vector3d,Unit3d,Angle
import numpy as np
from scipy.optimize import curve_fit


def fit(x,a,b,c,d):
    return a*x**3 + b*x**2 + c*x + d

def main():


    eye = len.Eye(0.0)
    print("Original Focalength : " + str(eye.backFocalLength()))
    eye.accommodation(1.5)
    print("Modifiec Focalength : " + str(eye.backFocalLength()))
    
    
    acc = np.linspace(1.0,1.5,50)
    focal = np.zeros(acc.size)
    
    for i,a in enumerate(acc):
        eye.accommodation(a)
        focal[i] = eye.backFocalLength(wl.PhotopicPeak)
        
    popt,pvar = curve_fit(fit,focal,acc)
    print(str(popt))
    
        
    plt.plot(focal,acc)
    plt.plot(focal,fit(focal,*popt))
    plt.show()
    
    
    
main()