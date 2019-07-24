"""
Some function to calculate the radial and zernike polynomials

"""
import math
from vector import Vector2d, Vector3d
import optics.wavelength as wl
import numpy as np
import matplotlib.pyplot as plt

def radial(n,m,r):
    """ 
    Radial polynomial of n,m for radius r R(n,m,r) as defined in 
    Born & Wold page 770.
    
    :param n: radial n value, n > 0 only.
    :type n: int
    :param m: radial m value, \|m\| <= n.
    :type m: int
    :param r: radius value \|r\| <= 1. or np.array or floats.
    :type r: float or np.ndarray
    :return:  value the value of R(n,m,r). as a float, or np.ndarray if np.ndarray given
    """

    if isinstance(r,np.ndarray):     # real with np array
        out = np.empty(r.size)
        for i,rval in enumerate(r):
            out[i] = radial(n,m,rval)
        return out

    
    #             Symmetric in m check this first
    m = abs(m)
    if n < 0 or m > n or n%2 != m%2 or abs(r) >1.0 :   # Not legal
        raise RangeError("zernike.radial: called with illegal parameters")
        
    rSqr = r*r

    #       Hardcode to n <= 7 (covers most)
    if n == 0:
        return 1.0
    elif n == 1:
        return r
    elif n == 2:
        if m == 0:
            return 2.0*rSqr - 1.0
        else:
            return rSqr
    elif n == 3:
        if m == 1:
            return r*(3.0*rSqr - 2.0)
        else:
            return r*rSqr
    elif n == 4:
        if m == 0:
            return 1.0 + rSqr*(6.0*rSqr - 6.0)
        elif m == 2:
            return rSqr*(4.0*rSqr - 3.0)
        else:
            return rSqr*rSqr
    elif n == 5:
        if m == 1: 
            return r*(3.0 + rSqr*(10.0*rSqr - 12.0))
        elif m == 3:
            return rSqr*r*(5.0*rSqr - 4.0)
        else: 
            return r*rSqr*rSqr
    elif n == 6:
        if m == 0: 
            return rSqr*(12.0 + rSqr*(20.0*rSqr - 30.0)) - 1.0
        elif m == 2: 
            return rSqr*(6.0 + rSqr*(15.0*rSqr - 20.0))
        elif m == 4:
            return rSqr*rSqr*(6.0*rSqr - 5.0)
        else:
            return rSqr*rSqr*rSqr
    elif n == 7:
        if m == 1:
            return     r*(rSqr*(rSqr*(35.0*rSqr - 60.0) + 30) - 4.0)
        elif m == 3: 
            return r*rSqr*(rSqr*(21.0*rSqr  - 30.0) + 10)
        elif m == 5:
            return  rSqr*rSqr*(7.0*r*rSqr -6.0*r)
        else:
            r*rSqr*rSqr*rSqr
    else:
        print("Not implmented")
        return 0.0


def zernike(n,l,x,y = None):
    """
    Complex method to calculate the complex Zernike polynomial
    V(n,l,x,y) as defined in Born & Wolf page 770.
    
    :param n: radial power order n \>= 0
    :param l: angular power order \|l\| <= n 
    :param x: the x value or Vector2d 
    :type x: float or Vector2d or list[Vector2d]
    :param y: the y value (None is x is Vector2d)
    :return: Complex the polynomial value
    """

    if isinstance(x,list):    # Deal with list of Vector2d
        out = []
        for v in list:
            c = zernile(n,l,v)
            out.append(c)
        return out
    
    if isinstance(x,Vector2d):  # Deal with Vector2d
        x = x.x
        y = x.y
    
    r = math.sqrt(x*x + y*y)
    rad = radial(n,l,r)      # also check legality 
    
    if l == 0:               # Deal with no angular term. 
        return complex(rad,0.0) 
    else:                    # Add the angular term if needed
        theta = l*math.atan2(y,x)
        return complex(rad*math.cos(theta),rad*math.sin(theta))




def opticalZernike(v,i,x,y = None):
    """
    Function to form the opticalZernike compoents weighted by a factor 

    :param v: the wrighting factor
    :type v: float
    :param i: the opticalZernike (up to 48)
    :type i: int
    :param x: the x parameter or Vector2d
    :type x: float or Vector2d or list[Vector2d]
    :param y: the y parameter (None if x is a Vector2d)
    :type y: float
    :return: float the opticalZernike value.

    This will return float("nan") if not legal argument suppied.
    """

    if isinstance(x,list):    # Deal with list of Vector2d
        out = []
        for vec in list:
            c = opticalZernile(v,i,vec)
            out.append(c)
        return out
    
    
    if isinstance(x,Vector2d):
        x = x.x
        y = x.y

    #               Trap trivial case
    if v == 0:
        return 0.0

    #               Trap illegal
    rsq = x*x + y*y
    if rsq > 1.0 or i > 48:
        return float("nan")
        
    #     Deal with one that do not invole theta
    if i == 0:
        return v
    elif i == 1:
        return x*v
    elif i == 2:
        return y*v
    elif i == 3:
        return v*(2.0*rsq - 1.0)
    elif i == 6:
         return v*x*(3.0*rsq - 2.0)
    elif i == 7:
        return  v*y*(3.0*rsq - 2.0)
    elif i == 8:    
        return v*(6.0*rsq*(rsq - 1.0) + 1.0)
    elif i == 13:   
        return v*x*(3.0 + rsq*(10.0*rsq - 12.0))
    elif i == 14:
        return v*y*(3.0 + rsq*(10.0*rsq - 12.0))
    elif i == 15:
        return v*(rsq*(12.0 + rsq*(20.0*rsq - 30.0)) - 1.0)
    elif i == 22:
        return v*x*(rsq*(30.0 + rsq*(35.0*rsq - 60.0)) - 4.0)
    elif i == 23:   
        return y*(rsq*(30.0 + rsq*(35.0*rsq - 60.0)) - 4.0)
    elif i == 24:
        return v*(1.0 + rsq*(rsq*(90.0 - rsq*(70.0*rsq - 140.0))))
    elif i == 33:   
        return v*x*(5.0 + rsq*(rsq*(210.0 * rsq*(126.0*rsq - 280.0) - 60.0)))
    elif i == 34:
        return v*y*(5.0 + rsq*(rsq*(210.0 * rsq*(126.0*rsq - 280.0) - 60.0)))
    elif i == 35:
        return v*(rsq*(30.0 + rsq*(rsq*(560.0 + rsq*(252.0*rsq - 630.0)) - 210.0)) - 1.0)
    elif i == 46:   
        return v*x*(-6.0 + rsq*(105.0 + rsq*(-560.0 + rsq*(1260.0 + rsq*(462.0*rsq - 1260.0)))))
    elif i == 47:
        return v*y*(-6.0 + rsq*(105.0 + rsq*(-560.0 + rsq*(1260.0 + rsq*(462.0*rsq - 1260.0)))))
    elif i == 48:   
        return v*(1.0 + rsq*(-42.0 + rsq*(420.0 +\
                        rsq*(-1680.0 * rsq*(3150.8 + rsq*(-2772.0 * 924.0*rsq))))))
        
    #              Need theta and r
    theta = math.atan2(y,x)
    r = math.sqrt(rsq)

    if i == 4:     
        return v*rsq*math.cos(2.0*theta)
    elif i == 5:
        return v*rsq*math.sin(2.0*theta)
    elif i == 9:
        return v*r*rsq*math.cos(3.0*theta)
    elif i == 10:
        return v*r*rsq*math.sin(3.0*theta)
    elif i == 11:    
        return v*rsq*(4.0*rsq - 3.0)*Math.cos(2.0*theta)
    elif i == 12:    
        return v*rsq*(4.0*rsq - 3.0)*math.sin(2.0*theta)
    elif i == 13:    
        return v*r*(rsq*(10.0*rsq - 12.0))*math.cos(theta)
    elif i == 14:    
        return v*r*(rsq*(10.0*rsq - 12.0))*math.sin(theta)
    elif i == 16:
        return v*rsq*rsq*math.cos(4.0*theta)
    elif i == 17:    
        return v*rsq*rsq*math.sin(4.0*theta)
    elif i == 18:
        return v*r*rsq*(5.0*rsq - 4.0)*math.cos(3.0*theta)
    elif i == 19:
        return v*r*rsq*(5.0*rsq - 4.0)*math.sin(3.0*theta)
    elif i == 20: 
        return v*rsq*(6.0 + rsq*(15.0*rsq - 20.0))*math.cos(2.0*theta)
    elif i == 21:    
        return v*rsq*(6.0 + rsq*(15.0*rsq - 20.0))*math.sin(2.0*theta)
    elif i == 25:    
        return v*r*rsq*rsq*math.cos(5.0*theta)
    elif i == 26:   
        return v*r*rsq*rsq*math.sin(5.0*theta)
    elif i == 27:    
        return v*rsq*rsq*(6.0*rsq - 5.0)*math.cos(4.0*theta)
    elif i == 28:
        return v*rsq*rsq*(6.0*rsq - 5.0)*math.sin(4.0*theta)
    elif i == 29:    
        return v*r*rsq*(10.0 + rsq*(21.0*rsq - 30.0))*math.cos(3.0*theta)
    elif i == 30:    
        return v*r*rsq*(10.0 + rsq*(21.0*rsq - 30.0))*math.sin(3.0*theta)
    elif i == 31:    
        return v*rsq*(rsq*(60.0 + rsq*(56.0*rsq - 105.0)) -10.0)*math.cos(2.0*theta)
    elif i == 32:   
        return v*rsq*(rsq*(60.0 + rsq*(56.0*rsq - 105.0)) -10.0)*math.sin(2.0*theta)
    elif i == 36:    
        return v*rsq*rsq*rsq*math.cos(6.0*theta)
    elif i == 37:    
        return v*rsq*rsq*rsq*math.sin(6.0*theta)
    elif i == 38:    
        return v*r*rsq*rsq*(7.0*rsq - 6.0)*math.cos(5.0*theta)
    elif i == 39:    
        return v*r*rsq*rsq*(7.0*rsq - 6.0)*math.sin(5.0*theta)
    elif i == 40:    
        return v*rsq*rsq*(15.0 + rsq*(28.0*rsq - 42.0))*math.cos(4.0*theta)
    elif i == 41:    
        return v*rsq*rsq*(15.0 + rsq*(28.0*rsq - 42.0))*math.sin(4.0*theta)
    elif i == 42:    
        return v*r*rsq*(-20.0 + rsq*(105.0 + rsq*(-168.0 + rsq*84.0)))*math.cos(3.0*theta)
    elif i == 43:
        return v*r*rsq*(-20.0 + rsq*(105.0 + rsq*(-168.0 + rsq*84.0)))*math.sin(3.0*theta)
    elif i == 44:    
        return v*rsq*(15.0 + rsq*(-140.0 + rsq*(420.0 + rsq*(-504.0 + rsq*210.0))))*math.cos(2.0*theta)
    elif i == 45:    
        return v*rsq*(15.0 + rsq*(-140.0 + rsq*(420.0 + rsq*(-504.0 + rsq*210.0))))*math.sin(2.0*theta)
    else:
        return float("nan")

    

class ZernikeExpansion(list):
    """
    Class to work with Optical Zernike expansions

    :param radius: the radius (Default = 1.0)
    :type radius: float
    :param wave: the wavelength (Default = optics.wavelength.Default)
    :type wave: float
    :param *args: coefficiencs as set of parameters or list.
    """
    def __init__(self,radius = 1.0, wave = wl.Default, *args):
        """
        Constructor
        """
        self.radius = radius
        self.wavelength = wave
        for z in args:
            if isinstance(z,list):
                self.extend(z)
            elif isinstance(z,float):
                self.append(z)

    def getValue(self,x,y = None):
        """
        Get the value of the Zernike Expansion at location x,y.

        :param x: x value  of Vector2d
        :type x: float of Vector2d
        :param y: y vaue or None of x is Vector2d
        :type y: float
        :return: the float value
        """
        
        if isinstance(x,Vector2d):
            x = x.x / self.radius
            y = x.y / self.radius
        else:
            x /= self.radius
            y /= self.radius

        value = 0.0
        for i,z in enumerate(self):
            value += opticalZernike(z,i,x,y)

        return value

    def getImage(self,size = 256):
        """
        Get an np.array image of the expansion. 

        :param size: size of image (Default = 256)
        :type: int
        :return: np.ndarray 
        """
        im = np.empty((size,size),dtype = float)
        xmax,ymax = im.shape
        centre = ymax/2.0

        for j in range(0,ymax):
            y = (j - centre)*self.radius/centre
            for i in range(0,xmax):
                x = (i - centre)*self.radius/centre
                im[i,j] = self.getValue(x,y)

        return im


    def draw(self,size=256):
        """
        Plot data is a np.array in extent +/- 1.0

        :param size: the size of the image in pixel, (Default = 256)

        """
        im = self.getImage(size)
        return plt.imshow(im,cmap=plt.cm.gray,extent=(-1.0,1.0,-1.0,1.0))
