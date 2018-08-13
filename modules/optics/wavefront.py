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


class Seidel(WaveFront):
    """
    Class to calcualte the basic Siedel aberrations
    """

    def __init__(self,coeff,theta = 0.0,radius = 1.0,wavelength = w.Default):
        """
        Form the siedel class with the coefficients 
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
        Get the value as specified poistion and field angle
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

class Zernike(WaveFront):
    """
    Zernike expansion of a wavefrom in terms of the opticalZernike compoents.
    """
    def __init__(self,coeff,radius = 1.0,wavelength = w.Default):
        """
        Form the Zernike class with the coefficients 
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
        xtilt *= math.tau/self.wavefront.maxradius
        ytilt *= math.tau/self.wavefront.maxradius

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
