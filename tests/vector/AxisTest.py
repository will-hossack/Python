
"""
vector and axis tests
"""

from vector import Vector3d,Axis3d,Angle,Unit3d
import tio as t
from optics.ray import IntensityRay,RayPath,RayPencil,Disc
from optics.surface import CircularAperture
import math
import matplotlib.pyplot as plt


def main():
    
    a = Vector3d(1,2,3)
    t.tprint(a)    
    
    origin = Vector3d(0,0,0)
    
    axis = Axis3d(origin)
    
    b = axis.transform(a)
    
    t.tprint(b)
    
    
    x = Unit3d(1,0,0)
    y = Unit3d(0,1,0).rotateAboutX(math.cos(math.radians(30)))
    z = Unit3d(0,0,1).rotateAboutX(math.cos(math.radians(30)))
    
    axis = Axis3d(origin,x,y,z)
    
    b = axis.transform(a)
    t.tprint(b)
    
    pt = Vector3d(0,0,50)
    ca = Disc(pt,10)
    u = Unit3d(Angle(math.radians(0)))
    pencil = RayPencil().addBeam(ca,u).rotateAboutX(math.radians(-40),pt)
    
    
    
    
    #for r in pencil:
     #   r.rotateAboutX(math.radians(-40),pt)
        #r.position -= pt
        #r.position.rotateAboutX(math.radians(-40))
        #r.position += pt
        #r.director.rotateAboutX(math.radians(-40))
      
        
      
    pencil.addMonitor(RayPath())    
    for r in pencil:
        t.tprint(50*r.director)
        
    pencil += 50
    
    #ca.draw()
    pencil.draw()
    plt.show()
    
    
    """
    
    pt = Vector3d(0,3,4)
    
    p = axis.transform(pt)
    v = axis.transform(u)
    
    ray = IntensityRay(p,v).addMonitor(RayPath())
    ray += 100
    ray.draw()
    plt.show()
    """
    
    
main()

