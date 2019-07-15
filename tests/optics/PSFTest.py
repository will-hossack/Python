"""
Test of PSF funcctions
"""
import math
import vector as v
import random
import matplotlib.pyplot as plt
import optics.psf as psf


def main():

    xData = []
    yData = []
    vData = []

    rad = 5.0
    for i in range(5000):
        while True:
            x = random.uniform(-rad,rad)
            y = random.uniform(-rad,rad)
            if x*x + y*y <= rad*rad:
                break
        pt = v.Vector2d(x,7*y)
        pt.rotate(math.radians(-30))
        pt += v.Vector2d(20,30)
        xData.append(pt.x)
        yData.append(pt.y)
        vData.append(pt)
        

    mom = psf.FixedMoments(vData)

    c = mom.centroid()
    print("Centriod : " + str(c))
    radius = mom.radius()
    print("Radius : " + str(radius))
    area = math.pi*radius*radius
    print("First area : " + str(area))

    major,minor,alpha = mom.ellipse()

    print("Major is : " + str(major))
    print("Minor is : " + str(minor))
    print("Alpha is : " + str(math.degrees(alpha)))

    area = mom.area()
    print("Second area : " + str(area))
    e = mom.eccentricity()
    print("Ecc is : " + str(e))
    
    fig = plt.figure()
    panel = fig.add_subplot(1,1,1)
    panel.axis('equal')
    plt.scatter(xData,yData)



    n = 20
    dtheta = 2*math.pi/n
    xval = []
    yval = []

    for i in range(0,n + 1):
        theta = i*dtheta
        p = v.Vector2d(major*math.cos(theta),minor*math.sin(theta))
        p.rotate(-alpha)
        p += c
        xval.append(p.x)
        yval.append(p.y)

    plt.plot(xval,yval,"k")

    
    plt.show()


main()
