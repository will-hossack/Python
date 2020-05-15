"""
    Test for scipy intregration
"""

import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt

def circle(x,y,r):
    if x*x + y*y <= r*r:
        return 1.0
    else:
        return 0.0


def shear(x,y,r,a,b):
    one = circle(x,y,r)
    two = circle(x+a,y+b,r)
    return one*two





def fn(x, a, b):
    """
    Function with two arguments
    """

    return  a*x**2 + b


def dfn(x,y):
    return float(x*y)

def ymax(y):
    return 1.0 - 2*y


def main():

    
    

    a = 3.0
    b = 2.0
    lower = 1.0
    upper = 4.0
    
    # val = integrate.quad(fn,lower,upper,args = (a,b))
    
    #print("Integration is : " + str(val))
    

    
    radius = 3.0
    
    sData = np.linspace(0.0,2*radius,20)
    aData = np.zeros(sData.size)
    
    
    xyRange = np.linspace(-radius,radius,50)
    
    
    
    
    
    
    for i,s in enumerate(sData):
        area = 0
        for y in xyRange:
            for x in xyRange:
                area += shear(x,y,radius,s,0.0) 


        #           area = integrate.dblquad(shear,-radius,radius,-radius,radius,args = (radius,x,y))
        aData[i] = area
        
     
    plt.plot(sData,aData)
    plt.show()

main()