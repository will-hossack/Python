"""
Set of classes to handle wavefronts, analysis and interferoneters.

Author: Will Hossack, The University of Edinburgh
"""
from vector import Vector2d,Vector3d,Unit3d
import math
import optics.wavelength as w
from optics.zernike import opticalZernike
import numpy as np
import matplotlib.pyplot as plt
import tio


class WaveFront(object):
    """
    Basic wavefront class, abstract at this point
    """
    def __init__(self,radius = 1.0, wavelength = w.Default):
        """
        Constuctor to form WaveFront object, only wavelength and maxradius is set 
        here.
        """
        self.wavelength = float(wavelength)
        self.maxradius = float(radius)

    def getValue(self,x,y):
        """
        Defaults getValue, to we supplied in extending classes.
        """
        return float("nan")

    def readFromFile(self,fn = None):
        """
        Read a wavefront from a file
        """
        if fn == None:
            wfile = tio.openFile("Wavefront file","r","wf")
        else:
            fn = tio.getExpandedFilename(fn)   # Sort out logicals
            if not fn.endswith("wf"):        # Append ".lens" if not given
                fn += ".wf"
            wfile= open(fn,"r")             # open file

        #          read file and process one line at a time
        #
        coef = []                          # Local coefficients
        wave = self.wavelength
        rad = self.maxradius
        type ="zzz"
        theta = 0.0
        for line in wfile.readlines():
            line = line.strip()
            if not line.startswith("#") and len(line) > 0:   # Kill comments and blanks
                token = line.split()

                if token[0].startswith("type"):
                    type = str(token[1])
                elif token[0].startswith("radius"):
                    rad = float(token[1])
                elif token[0].startswith("wave"):
                    wave = float(token[1])
                elif token[0].startswith("field"):   # Only for Seidel
                    theta = float(token[1])      
                else: # Assume a coefficient
                    c = float(token[0])
                    coef.append(c)

        #          Work cout what we have
        if type.startswith("seid"):
            return SeidelWaveFront(coef,theta,rad,wave)
        elif type.startswith("zernike"):
            return ZernikeWaveFront(coef,rad,wave)
        elif type.startswith("poly"):
            return PolynomialWaveFront(coef,rad,wave)
        else:
            print("WaveFront.readFromFile: unknown type : " + str(type))
            return None
                

class SeidelWaveFront(WaveFront):
    """
    Class to calcualte the basic Siedel aberrations
    """

    def __init__(self,coeff,theta = 0.0,radius = 1.0,wavelength = w.Default):
        """
        Form the siedel class with the coefficients, note the coefficients are in microns.
        """
        WaveFront.__init__(self,radius,wavelength)
        self.coef = coeff
        self.theta = float(theta)

        self.name = ("Defocus","Spherical Aberration","Coma","Astigmatism","Field Curvature","Distortion")

    def __str__(self):
        """
        The str to print out values on single line.
        """
        s = "Seidel: S0: {0:5.3f} S1: {1:5.3f} S2: {2:5.3f} S3: {3:5.3f} S4 : {4:5.3f} S5: {5:5.3f} Theta: {6:5.3f}". \
            format(self.coef[0],self.coef[1],self.coef[2],self.coef[3],self.coef[4],self.coef[5],self.theta)
        return s

    def __repr__(self):
        """   Detailed of the Seidel in tabular form.
        """
        s = "Seidel Aberrations \n"
        for i in range(0,len(self.coef)):
            s += "{0:<24}: {1:7.5f}\n".format(self.name[i],self.coef[i])
        s += "{0:<24}: {1:7.5f}".format("Field Angle",self.theta)
        return s
        

    def getValue(self,x,y):
        """
        Get the value of phase at specified poistion and field angle.
        """
        x /= self.maxradius              # Normalise
        y /= self.maxradius

        rSqr = x*x + y*y
        if rSqr > 1.0:
            return float("nan")

        phi = 0.5*self.coef[0]*rSqr + 0.125*self.coef[1]*rSqr*rSqr
        if self.theta != 0.0:
            phi += 0.5*self.coef[2]*y*rSqr*self.theta + \
                   0.5*self.coef[3]*y*y*self.theta*self.theta + \
                   0.25*(self.coef[3] + self.coef[4])*rSqr*self.theta*self.theta + \
                   0.5*self.coef[5]*y*self.theta**3
        return 2.0*math.pi*phi/self.wavelength

class ZernikeWaveFront(WaveFront):
    """
    Zernike expansion of a wavefrom in terms of the opticalZernike compoents.
    """
    def __init__(self,coeff,radius = 1.0,wavelength = w.Default):
        """
        Form the Zernike class with the coefficients, coefficeints are in microns.
        """
        WaveFront.__init__(self,radius,wavelength)
        self.coef = coeff

        self.name = ("bais","x-tilt","y-tilt","defocus","x-astigmatism","y-astigmatism",\
           "x-coma","y-coma","primary spherical","x-trefoil","y-trefoil",\
           "secondary x-Astigmatism","secondary y-Astigmatism",\
           "secondary x-coma", "secondary y-coma","secondary spherical",\
           "x-tetrafoil","y-tretafoil","secondary x-trefoil","secondary y-trefoil",\
           "tertiary x-astigmatism","tertiary y-astigmatism",\
           "tertiary x-coma","tertiary y-coma","terniary spherical",\
           "x-pentafoil","y-pentafoil","secondary x-tetrafoil","secondary y-tetrafoil",\
           "tertiary x-trefoil","tertiary x-trefoil","quatenary x-astigmatism",\
           "quatenary y-astigmatism","quatenary x-coma","quatenary y-coma",\
           "quaternary spherical")

    def __repr__(self):
        """   Detailed of the Zernike in tabular form.
        """
        s = "Zernike Expansion\n"
        for i in range(0,len(self.coef)):
            s += "{0:<24}: {1:7.5f}\n".format(self.name[i],self.coef[i])
        return s

    def getValue(self,x,y):
        """
        Get the value as specified poistion and field angle
        """
        x /= self.maxradius              # Normalise
        y /= self.maxradius

        rSqr = x*x + y*y
        if rSqr > 1.0:
            return float("nan")

        phi = 0.0
        for i in range(0,len(self.coef)):
            phi += opticalZernike(self.coef[i],i,x,y)
            
        return 2.0*math.pi*phi/self.wavelength
    

class PolynomialWaveFront(WaveFront):
    """
    Polynomial expansion of a wavefront in terms of x/y
    """
    def __init__(self,coeff,radius = 1.0,wavelength = w.Default):
        """
        Form the Polynomial class with the coefficients, coefficeints are in microns.
        """
        WaveFront.__init__(self,radius,wavelength)
        self.coef = coeff

        self.name = ("bias","xtilt","ytilt","x^2","xy","y^2","x^3","x^2 y","x y^2","y^3"\
                     "x^4","x^3 y","x^2 y^2","x y^3","y^4")

    def __repr__(self):
        """   Detailed of the Zernike in tabular form.
        """
        s = "Polynomial Expansion\n"
        for i in range(0,len(self.coef)):
            s += "{0:<24}: {1:7.5f}\n".format(self.name[i],self.coef[i])
        return s

    def polynomial(self,v,i,x,y):
        """
        Internal function to form one polynomial complonent
        """    

        # Trap trivial case
        if v == 0.0:
            return 0.0

        if x*x + y*y > 1.0:
            return float("nan")

        #         hard code the basic set
        if i == 0:
            return v
        elif i == 1:      # First order
            return x*v
        elif i == 2:
            return y*v
        elif i == 3:      # Second order
            return v*x*x
        elif i == 4:
            return v*x*y
        elif i == 5:
            return v*y*y
        elif i == 6:       # Third order
            return v*x**3
        elif i == 7:
            return v*x*x*y
        elif i == 8:
            return v*x*y*y
        elif i == 9:
            return v*y**3
        elif i == 10:     # Fourth order
            return v*x**4
        elif i == 11:
            return v*x**3*y
        elif i == 12:
            return v*x**2*y**2
        elif i == 13:
            return v*x*y**3
        elif i ==14:
            return v*y**4
        else:
            return float("nan")
        
        
    def getValue(self,x,y):
        """
        Get the value as specified poistion and field angle
        """
        x /= self.maxradius              # Normalise
        y /= self.maxradius

        rSqr = x*x + y*y
        if rSqr > 1.0:
            return float("nan")

        phi = 0.0
        for i in range(0,len(self.coef)):
            phi += self.polynomial(self.coef[i],i,x,y)
            
        return 2.0*math.pi*phi/self.wavelength


#class WaveFrontPoint(object):
    """
    Class to hold the value of a wavefront at specified point 
    """

    

class Interferometer(object):
    """
    Defeine an interferometer to view wavefront aberrations
    """

    def __init__(self,size = 200, type = "Twyman"):
        """
        Create an interferometer to view fringes for a WaveFront 
        """

        self.size = int(size)
        self.type = str(type)

        #         Create the image as an np array of zeros
        self.image = np.zeros((self.size,self.size),dtype = float)

    def addWaveFront(self,wf,xtilt = 4.0, ytilt = 0.0, draw = True):
        """
        Add the wavefront and render if requested
        """
        self.wavefront = wf
        if draw:
            self.setTilt(xtilt,ytilt)
            self.draw()
       

    def setTilt(self, xtilt = 4.0, ytilt = 0.0):
        """
        Render the interferoneter with specified tilt
        """
        xtilt *= 2.0*math.pi/self.wavefront.maxradius
        ytilt *= 2.0*math.pi/self.wavefront.maxradius

        xr = range(0,self.size)

        for j in xr:
            y = 2.0*self.wavefront.maxradius*(j - self.size/2)/self.size
            for i in xr:
                 x = 2.0*self.wavefront.maxradius*(i - self.size/2)/self.size
                 phi =self.wavefront.getValue(x,y) + x*xtilt + y*ytilt
                 self.image[i,j] = math.cos(phi) + 1.0
                 
        

    def draw(self):
        """
        Display the image in pyplot
        """
        plt.title(self.type)
        plt.imshow(self.image,cmap=plt.cm.gray, \
                   extent=(-self.wavefront.maxradius,self.wavefront.maxradius,\
                           -self.wavefront.maxradius,self.wavefront.maxradius))


class ScalarPSF(object):
    """
    Class to form the Scalar PSF from a wavefront.
    """

    def __init__(self,size = 256):
        """
        Creat the object of specified size
        """
        self.size = int(size)
        self.image = np.zeros((self.size,self.size),dtype = complex)

    def addWaveFront(self,wf):
        """
        Add the wavefront, this will calsulate the ScarPSF and return the image,
        but will also hold in interallally so that it can be rendered
        """
        xr = range(0,self.size)

        for j in xr:
            y = 2.0*wf.maxradius*(j - self.size/2)/self.size
            for i in xr:
                 x = 2.0*wf.maxradius*(i - self.size/2)/self.size
                 if x*x + y*y <= wf.maxradius*wf.maxradius:
                     phi =wf.getValue(x,y) 
                     self.image[i,j] = complex(math.cos(phi),math.sin(phi))

        
        self.image = np.fft.fft2(self.image)
        self.image = np.fft.fftshift(self.image)
        self.image = np.absolute(self.image)
        return self.image

    def draw(self):
        """
        Render the ScalaSPF using imshow with grayscale colour map
        """
        plt.imshow(self.image,cmap=plt.cm.gray)



