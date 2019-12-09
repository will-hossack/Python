"""
Some function to calculate the radial and zernike polynomials

"""
import math
import cmath
from vector import Vector2d
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

def radial(n,m,r):
    """ 
    Radial polynomial of n,m for radius r R(n,m,r) as defined in 
    Born & Wold page 770.
    
    :param n: radial n value, n > 0 only.
    :type n: int
    :param m: radial m value, \|m\| <= n.
    :type m: int
    :param r: radius value \|r\| <= 1. or np.array or floats, all less that 1.0.
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
            c = zernike(n,l,v)
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


    
zernikeNames = ("Piston","X-tilt","Y-tilt","Defocus",\
                "X-astigmatism","Y-astigmatism","X-coma","Y-coma","Primary Spherical",\
                "X-trefoil","Y-Trefoil","Secondary X-Astig","Secordary Y-Astig","Secondary X-coma",\
                "Secondary Y-coma","Secondary Spherical",\
                "X-tetrafoil","Y-tertafoil","Secondary X-trefoil","Secondary Y-trefoil","Terniary X-Astig","Terniaty Y-Astrig",\
                "Terniary X-coma","Ternairy Y-coma","Terniary Spherical",\
                "X-pentafoil","Y-pentafoil","Secondary X-tetrafoil","Secondary Y-tertafoil",\
                "Terniary X-trefoil","Terniary Y-trefoil","Quatenary X-Astig","Quantenary Y-Astig",\
                "Quantary X-coma","Quantary y-coma","Quantary Spherical")

def opticalZernikeName(i):
    """
    Function to lookup and retuirn the name of the specified optical zernike component
    """
    if i < len(zernikeNames):
        return zernikeNames[i]
    else:
        return "Zernike component {0:d} not specified".format(i)
    

def opticalZernike(v,i,x,y = None):
    """
    Function to form the opticalZernike components weighted by a factor. These are defined for order 0 to 48. 

    :param v: the weighting factor
    :type v: float
    :param i: the opticalZernike 0 to 48 (49 components)
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
            c = opticalZernike(v,i,vec)
            out.append(c)
        return out
    
    
    if isinstance(x,Vector2d):
        y = x.y
        x = x.x
   
    #               Trap illegal
    rsq = x*x + y*y
    if rsq > 1.0 or i > 48:
        return float("nan")

    if v == 0.0:
        return 0.0
        
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
        return v*rsq*(4.0*rsq - 3.0)*math.cos(2.0*theta)
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
    Class to hold a zernike expansion, being a list of optical zernike components. There are also method to evaluate and
    display the expansion.

    :param radius: the radius (Default = 1.0)
    :type radius: float
    :param \*args: coefficiencs as set of parameters or list, may be blank.

    """
    def __init__(self,radius = 1.0, *args):
        """
        Constructor
        """
        self.radius = radius
        for z in args:
            if isinstance(z,list):
                self.extend(z)
            elif isinstance(z,float):
                self.append(z)


    def __str__(self):
        """
        Print out list inclduing the component names.
        """ 
        s =  "r: {0:6.4f}".format(self.radius)
        for i in range(len(self)):
            s += "\n{0:s} : {1:8.4e}, ".format(opticalZernikeName(i),self[i])
        
        return s

    def __repr__(self):
        """
        Return repr of class, being class name + str(self)
        """
        return "{0:s} ".format(self.__class__.__name__) + str(self)
        
    def getValue(self,x,y = None):
        """
        Get the value of the Zernike Expansion at location x,y, note the x/y values are divided by radius before evaluation.

        :param x: x value  of Vector2d
        :type x: float of Vector2d
        :param y: y vaue or None of x is Vector2d
        :type y: float
        :return: the float value

        Note: if x/y outside range (do outside circle specifed by self.radius) this will return "NaN".
        """
        
        if isinstance(x,Vector2d):
            y = x.y / self.radius
            x = x.x / self.radius
        else:
            x /= self.radius
            y /= self.radius

        value = 0.0
        for i,z in enumerate(self):
            value += opticalZernike(z,i,x,y)

        return value

    def getImage(self,size = 256, xtilt = None, ytilt = None):
        """
        Get an np.array image of the expansion. If both tilts are None, 
        then image set to raw phase value, othwise
        will be simulated interferometer with specified tilt with output in the range 0.0 -> 2.0.

        Note: the pixel elements outside the unit circle will be set to Zero.

        :param size: size of image (Default = 256)
        :type: int
        :param xtilt: Interferometer xtilt, may be None
        :param ytilt: Interferometer ytilt, may be None
        :return: two dimensional np.ndarray  

        """
        im = np.empty((size,size),dtype = float) # Empty array

        if xtilt == None and ytilt == None:      # Sort out the fringe.
            fringe = False
        else:
            fringe = True
            if xtilt == None:
                xtilt = 0.0
            if ytilt == None:
                ytilt = 0.0
            
        xmax,ymax = im.shape
        ycentre = ymax/2.0
        xcentre = xmax/2.0

        for j in range(0,ymax):
            y = (j - ycentre)*self.radius/ycentre       # In range -1.0 to 1.0
            for i in range(0,xmax):
                x = (i - xcentre)*self.radius/xcentre   # In range -1.0 to 1.0
                v = self.getValue(x,y)
                if not math.isnan(v) and fringe:      # Fringe if valid
                    v = 1.0 + math.cos(2.0*math.pi*(x*xtilt + y*ytilt) + v)
                im[i,j] = v                           # Note will be NaN if outside unit circle.

        return im


    def getPSF(self,size = 256, log = True):
        """
        Get the PSF by fourier means.

        :param size: size of PDF array, Default = 256
        :type size: int
        :param log: logical if log of intensity taken, Default = True
        :type log: bool
        """
        im = self.getImage(size)
        r = np.cos(im)     # Real part
        i = np.sin(im)     # Imaginary part
        #        z = np.vectorize(complex)(r,i)
        z = r + 1j*i       # Combine to form complex array

        z = np.nan_to_num(z)
        psf = np.fft.fft2(z)
        psf = np.fft.fftshift(psf)  # Shift to centre
        psf = abs(psf)
        if log :                # Take the log if required.
            psf = np.log(psf + 1.0)

        return psf


    def getOTF(self,size = 256, horizontal = True):
        """
        Get the one-dimenensioal normalsied OFT as np array
        """
        im = self.getImage(size)     # Get the phase image
        # Make the complex image and it complex conjugate
        r = np.cos(im)
        i = np.sin(im)
        z = r + 1j*i
        zc = np.conj(z)

        xsize,ysize = im.shape

        if horizontal:              # Sort out ditection of shift
            shiftSize = xsize
            fullRange = range(0,ysize)
        else:
            shiftSize = ysize
            fullRange = range(0,xsize)
        

        otfData = np.zeros(shiftSize)  # np array to hold the OFT

        #      Loop ovre the shifts
        for shift in range(0,shiftSize): 
            otf = 0.0
            shiftRange = range(shift,shiftSize)

            #
            #       Do 2d sum over overlap area
            for i in shiftRange:
                for j in fullRange:
                    ish = i - shift
                    if horizontal:
                        zr = z[i,j]
                        zl = zc[ish,j]
                    else:
                        zr = z[j,i]
                        zl = zc[j,ish]
                    if not (cmath.isnan(zr) or cmath.isnan(zl)) :
                        otf += (zr * zl).real
            otfData[shift] = otf

        #        Normalise this output
        max = otfData[0]
        otfData /= max
        return otfData    
        



    def plotOTF(self,size = 256, horizontal = True):
        """
        Calcualte and plot OFT with sensible plot paramters.
        """
        otfData = self.getOTF(size,horizontal)   # Get the OTF
        shiftData = np.linspace(0.0,1.0,otfData.size)
        plt.plot(shiftData,otfData)
        plt.xlim(0.0,1.0)
        plt.grid()
        plt.xlabel("Normalised spatial frequency")
        plt.ylabel("OFT")
        plt.title("Plot of OTF")
        return otfData                  # Return otf Data in case it is needed

    def draw(self,size = 256 ,xtilt = None, ytilt = None):
        """
        Plot data is a np.array in extent +/- 1.0

        :param size: the size of the image in pixel, (Default = 256)
        :type size: int
        :param xtilt: X-Interferometer tilt, if None then plot raw phase values
        :type xtilt: float of None
        :param ytilt:  Y-Interferometer tilt
        :type ytilt: float or none

        """
        im = self.getImage(size,xtilt,ytilt)
        plt.imshow(im,cmap=plt.cm.gray,extent=(-1.0,1.0,-1.0,1.0))
