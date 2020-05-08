#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 09:13:07 2020

Program to explore stop diagrams on the back of the retina in the Eye class
"""
import optics.lens as l
import optics.ray as r
import tio as t
import vector as v
import optics.psf as psf
import matplotlib.pyplot as plt
import math

def main():
    
        lens = l.Eye()
        iris = t.getFloat("Iris",1.0)
        lens.setIris(iris)
        angle = t.getFloat("Angle in degrees",5.0)
        u = v.Unit3d(v.Angle(math.radians(angle)))
        
        vpencil = r.RayPencil().addBeam(lens,u,key="vl").addMonitor(r.RayPath())
        spencil = r.RayPencil().addBeam(lens,u,key="array")
        
        
        vpencil *= lens
        spencil *= lens
        plane = lens.getRetina()
        
        ps = psf.Psf().setWithRays(spencil,plane)
        t.tprint("PSF is",repr(ps))
        
        
        spot = psf.SpotDiagram(spencil)
        spot.draw(plane,True)
        
        plt.show()
        """       
        
        lens.draw()
        vpencil.draw()
        plt.axis("equal")
        plt.show()
        
        """
         
         
main()
