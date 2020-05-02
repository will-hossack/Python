#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 11:51:41 2020

@author: wjh
"""

import optics.surface as sur
import optics.ray as ray
import matplotlib.pyplot as plt
import vector as v
import optics.wavelength as w
from optics.lens import Prism
import math



def main():
    #w.setFixedAirIndex(True,1.3)
    
    prism = Prism(0,height = 100, index = "F4")
    mindev = prism.minDeviation(w.Green)
    print("Min deviation" + str(math.degrees(mindev)))
    
    # Form a circular aperture 
    ca = sur.CircularAperture([0,-35,-100],30)
    u = v.Unit3d(v.Angle(mindev/2))
    
    #          Build mulipcoloured raypencil
    pencil = ray.RayPencil().addBeam(ca,u,wave=w.Blue).addMonitor(ray.RayPath())
    pencil.addBeam(ca,u,wave=w.Green).addMonitor(ray.RayPath())
    pencil.addBeam(ca,u,wave=w.Red).addMonitor(ray.RayPath())
    print("Numnber of rays : " + str(len(pencil)))
    
    
    pencil *= prism
    pencil += 100
    
    angle = pencil[30].director.getAngle()
    print("Angle is:" + str(angle.getDegrees()))
    delta = 2*angle.theta - mindev
    print("Error is " + str(math.degrees(delta)))
    
    prism.draw()
    pencil.draw()
    ca.draw()
    plt.axis("equal")
    plt.show()
    
    
    
main()
