"""
Set of classes to implement Polarsiatisaton analysis via Jones methods
"""
import math
import cmath
from optics.wavelength import getDefaultWavelength,getDesignWavelength,WavelengthColour
from  optics.ray import Ray
from matplotlib.pyplot import polar,show
import numpy as np


def parseAngle(theta = 0.0):
    """
    Method to parse angle, for float / int assume radians, 
    """
    if isinstance(theta,float) or isinstance(theta,int):
        return float(theta)              
    elif(theta,str):
        if theta.startswith("h"):  # Horizontal
            return 0.0
        if theta.startswith("v"):  # Vertical
            return math.pi/2

        return math.radians(float(theta)) # Try and convert from degress
    else:
        print("Angle with unknown paramteer type")
        return 0.0


class JonesVector(Ray):
    """
    JonesVector to hold polarised state of light

    :param x: x complex component of jones vector (defaults to 1.0)
    :type x: complex
    :param y: component of jones vector (defults to 0.0)
    :type y: complex
    :param wavelength: Wavelength Jones vector (defaults to Default)

    Will also accept single argument of current JonesVector
    """
    def __init__(self,x = 1.0 ,y = 0.0 ,wavelength = getDefaultWavelength()):
        """
        Constructor to create a JonesVector

        """
        if isinstance(x,JonesVector):
            self.__init__(x.x,x.y,x.wavelength)
        else:
            Ray.__init__(self,wavelength)       # Set underlying wavelength
            self.x = complex(x)
            self.y = complex(y)

    
    def __str__(self):
        """
        Implment the str() method to print out x/y components and wavelength.
        """
        return "({0:s} , {1:s} , {2:7.5f} )".format(str(self.x),\
                                str(self.y),self.wavelength)
    #
    #
    def __repr__(self):
        """
        The repr function with the class name.
        """
        return "{0:s} ".format(self.__class__.__name__) + str(self)


    #       
    def copy(self):
        """
        Method to make copy and a new JonesVector being a copy of the current.
        
        :return: copy of current JonesVector

        """
        return JonesVector(self)    


    def __iadd__(self,b):
        """
        implment self += b, if b is  JonesVector then the components add, otherwise b is assumed constant
        and complex(b) is added to each component.
        """
        if isinstance(b,JonesVector):
            self.x += b.x
            self.y += b.y
        else:
            self.x += complex(b)
            self.y += complex(b)
        return self

    def __add__(self,b):
        """
        Implement the c = self + b where b is JonesVector then componenst add, else complex(b) 
        added to both componnets
        """
        c = self.copy()
        c += b
        return c

    def __radd__(self,b):
        """
        Implement the c = b + self where b is JonesVector then componnets add, else complex(b) is 
        added to both components.
        """
        c = self.copy()
        c += b
        return c


    def __isub__(self,b):
        """
        implment self -= b If b is a JonesVector then the componets are substracted, else complex(b) is subracted from each.
        """
        if isinstance(b,JonesVector):
            self.x -= b.x
            self.y -= b.y
        else:
            self.x -= complex(b)
            self.y -= complex(b)
        return self
        

    def __sub__(self,b):
        """
        Implement the c = self - b where b is JonesVector then componenst add, else complex(b) added to both componnets
        """
        c = self.copy()
        c -= b
        return c

    def __rsub__(self,b):
        """
        Implement the c = b - self where b is JonesVector then componnets add, else complex(b) is added to both components.
        """
        if isinstance(b,JonesVector):
            x = b.x - self.x
            y = b.y - self.y
        else:
            z = complex(b)
            x = z - self.x
            y = z - self.y
        return JonesVector(x,y,self.wavelength)



    def getIntensity(self):
        """
        Method to get the intensity of the current JonesVector.
        returns intensity as a float.

        :return: intensity as a float.
        """
        return  self.x.real**2 + self.x.imag**2 + self.y.real**2 + self.y.imag**2
    #
    #     
    def getPhase(self):
        """
        Method to get the phase difference between the two componets of 
        the current JonesVector

        :return: phase difference between the two componts in the range -pi ->p pi.
        """
        delta = cmath.phase(self.y) - cmath.phase(self.x)
        if delta > math.pi :
            delta -= 2*math.pi
        if delta < -math.pi:
            delta += 2*math.pi
        return delta


    def getAngle(self):
        """
        Method to get the angle of the polarsiation ellipse of the
        current JonesVector in the range given by atan2
        """
        xm = abs(self.x)
        ym = abs(self.y)
        pa = self.getPhase()
        return 0.5*math.atan2(2.0*xm*ym*math.cos(pa) , (xm*xm - ym*ym))

    def getEllipicity(self):
        """
        Method to get the ellipticity of the polarisation ellipse
        of the current JonesVector 

        :return: Ellipicity as a float on range -pi/2 -> pi/2
        """
        xm = abs(self.x)
        ym = abs(self.y)
        pa = self.getPhase()
        return 0.5*math.asin(2*xm*ym*math.sin(pa)/(xm*xm + ym*ym))

    
    def __imul__(self,m):
        """
        Implment a *= m  to pre-multiply JonesSystemMatrix a single JoneMatrix, othwise each component multiplied by complex(m)
        """

        if isinstance(m,JonesMatrixSystem) :     # Deal with system
            m = m.getMatrix()
            
        if isinstance(m,JonesMatrix):
            p = self.x*m.A + self.y*m.B
            q = self.x*m.C + self.y*m.D
            self.x = p
            self.y = q
        else:
            c = complex(m)
            self.x *= c
            self.y *= c
        return self

       
    def __mul__(self,m):
        """
        Implment the c = self * m tio pre-multiply by a single JoneMatrix, or complex(m)
        """
        c = self.copy()
        c *= m
        return c

    def __rmul__(self,b):
        """
        Implemnts c = b * self which multiplies each component by complex(b), is NOT valid for JonesVector
        """
        c = self.copy()
        d = complex(b)
        c.x *= d
        c.y *= d
        return c
        
            

    def throughPolariser(self,theta):

        """
        Method to put the current JonesVector through an idea  Linear Polarsier at specfied angle and
        get the output intensity.

        :param theta: angle of polarised from y-axis
        :type theta: float

        The current JoneVector is not altered.
        """
        pol = LinearPolariser(theta)
        r = self * pol
        return r.getIntensity()


   
    def polarPlot(self,points = 200):
        """
        Method to generate polar plot through an indeal rotated linear polarsied polariser.

        :param key: plot key to plt.plot, (Default = "r")
        :type key: str
        :param points: number of points (Default = 200)
        :type points: int
        """
        theta = np.linspace(0.0,2*math.pi,num=points,endpoint=True)
        intensity = np.empty(theta.size)

        pol = LinearPolariser()
        
        for i,angle in enumerate(theta):
            pol.setAngle(angle)
            b = self*pol
            intensity[i] = b.getIntensity()

        col = WavelengthColour(self.wavelength)
        polar(theta,intensity,col)
        


class LinearPolarisedBeam(JonesVector):
    """ 
    Class Linear Polarsied Beam

    :param theta: angle wrt x axis in radians (default  = 0.0)
    :type theta: float
    :param intensity: the intensity of the beam (default = 1.0)
    :type intensity: float
    :param wavelength: the wavelnegth (defaults to wl.Default)
    :type wavelength: float
     
    """
    def __init__(self,theta = 0.0, intensity = 1.0, wavelength = getDefaultWavelength()):
        """
        Linear polarised beam 
        """
        amp = math.sqrt(intensity)
        theta = parseAngle(theta)
        JonesVector.__init__(self,amp*math.cos(theta),amp*math.sin(theta),wavelength)


class RightCircularPolarisedBeam(JonesVector):
    """     
    Right Circular Polarsied Beam
    
    :param intensity: the intensity of the beam (defaults to 1.0)
    :type intensity: float
    :param wavelength: the weavelength (defaults to wl.Default)
    :type wavelength: float

    """
    def __init__(self,intensity = 1.0, wavelength = getDefaultWavelength()):
        """   
        Right circular polarsied beam
        
        """
        amp = math.sqrt(intensity/2.0)
        JonesVector.__init__(self,amp,complex(0.0,-amp),wavelength)

class LeftCircularPolarisedBeam(JonesVector):
    """    
    Left Circular Polarsied Beam

    :param intensity: the intensity of the beam (defaults to 1.0)
    :type intensity: float
    :param wavelength: the weavelength (defaults to wl.Default)
    :type wavelength: float

    """
    def __init__(self,intensity = 1.0, wavelength = getDefaultWavelength()):
        """   Left circular polarsied beam
        """
        amp = math.sqrt(intensity/2.0)
        JonesVector.__init__(self,amp,complex(0.0,amp),wavelength)


class JonesMatrix(object):
    """
    JonesMatrix class to implement a component matrix with four complex components, defaults to unit matrix.

    :param a: the A component (Default = 1.0)
    :type a: complex or Jonesmatrix
    :param b: The B component (Default = 0.0)
    :type b: complex
    :param c: The C compoent (Default = 0.0)
    :type c: complex
    :param d: the D component (Default = 1.0)

    """
    
    def __init__(self, a = 1.0, b = 0.0, c = 0.0, d = 1.0):
        """
        Set the values
        """
        if isinstance(a,JonesMatrix):
            self.set(a.A,a.B,a.C,a.D)
        else:
            self.set(a,b,c,d)


    def set(self,a , b, c, d):
        """
        Method to set (or reset)  the actual matrix values.
        
        :param a: The A component
        :type a: complex (parsed by complex(a))
        :param b: The B component
        :type b: complex (parsed by complex(b))
        :param c: The C component
        :type c: complex (parsed by complex(c))
        :param d: The D component
        :type d: complex (parsed by complex(d))
        
        """
        self.A = complex(a)
        self.B = complex(b)
        self.C = complex(c)
        self.D = complex(d)
        return self


    def setWavelength(self,w = None):
        """
        Method to set wavelength, abstract here, but overloaded in retarders
        """

    def setAngle(self,theta):
        """
        Method to set angle of component, abstarct here but overloaded where it make sense
        """
    
    def __str__(self):
        """
         The str method
         """
        return "[ " + str(self.A) + "  " + str(self.B) + " " + \
            str(self.C) + "  " + str(self.D) + "]"

    
    def __repr__(self):
        """
        The repr function with the class name.
        """
        return "{0:s} ".format(self.__class__.__name__) + str(self)
     
    def determinant(self) :
        """
        Form the determinant of Matrix 

        :return: the complex determinant
        """
        
        det = self.A*self.D - self.B*self.C
        return det
    
   
    def trace(self):
        """
        Form the compex trace of the matrix

        :return: the complex trace
        """
        tr = self.A + self.D
        return tr

    #      
    def copy(self):
        """
         Method to copy the current JonesMatrix

        :return: copy of current JonesMatrix
        """
        return JonesMatrix(self)


    def __imul__(self,m):
        """
        Implement the self *= m operators, will do complex mulgtiply of matrix if m is JonesMatrix, if not
        will multiply all elements by complex(m)
        """
        if isinstance(m,JonesMatrix):
            a = m.A*self.A + m.B*self.C
            b = m.A*self.B + m.B*self.D
            c = m.C*self.A + m.D*self.C
            d = m.C*self.B + m.D*self.D
            self.set(a,b,c,d)
        else:
            m = complex(m)
            self.A *= m
            self.B *= m
            self.C *= m
            self.D *= m
        return self

    def __mul__(self,m):
        """
        Implement the n = self * m operator
        """
        n = self.copy()
        n *= m
        return n
            


    def rotate(self,angle):
        """
        Method to rotate a general JonesMatix by angle using the pre and post rotation matrices and 
        returns a new JonesMatrix
        
        :param angle: Rotation angle in radians
        :type angle: float
        :return: a new JonesMatrix
        """
        cos = math.cos(angle)
        sin = math.sin(angle)
        #           First rotation matrix
        rp = JonesMatrix(cos,-sin,sin,cos)
        a = self * rp
        ra = JonesMatrix(cos,sin,-sin,cos)
        return ra * a

    
    def rotateBy(self,angle):
        """
        Method to rotate the current JonesMatrix by specified angle in place.

        :param angle: Rotation angle in radians
        :type angle: float
        :return: the modified Jonesmatrix
        """
        cos = math.cos(angle)
        sin = math.sin(angle)
        #           First rotation matrix
        rp = JonesMatrix(cos,-sin,sin,cos)
        self *= rp
        ra = JonesMatrix(cos,sin,-sin,cos)
        ra *= self
        self.set(ra.A,ra.B,ra.C,ra.D)        # set self with rotated values.

       

class LinearPolariser(JonesMatrix):
    """ 
    Class to form a LinearPolariser 

    :param theta: Polarsistaion axis in radians (Default = 0.0)
    :type theta: float
    :param transmission: Intensity transmission, (Default = 1.0)
    :type transmission: float
    """
    def __init__(self,theta = 0.0,transmission = 1.0):
        self.transmission = transmission
        self.setAngle(theta)

    def __str__(self):
        """
        The str to include transmission
        """
        return JonesMatrix.__str__(self) + "a: {0:6.4f} t: {1:6.4f}".format(self.angle,self.transmission)
        
    
    def copy(self):
        """
        Method to return a copy

        :return: copy of currennt LinearPolariser
        """
        return LinearPolariser(self.angle,self.transmission)

    
    def setAngle(self,theta = 0.0):
        """
        Method to set the angle of the current Linear Polarsier.

        :param theta: angle of the polarsier wrt to x-axis in radians (Default = 0.0)
        :type theta: float
        """
        self.angle = parseAngle(theta)      #   Record the angle
        cos = math.cos(self.angle)
        sin = math.sin(self.angle)
        amp = math.sqrt(self.transmission)
        self.set(amp*cos*cos,amp*cos*sin,amp*cos*sin,amp*sin*sin)

    def rotate(self,delta):
        """
        Method to return a new LinearPolarised with its angle incremended by specified angle 
        (overload of general method in JonesMatrix)

        :param delta: Angle increment in radians
        :type delta: float
        """
        return LinearPolariser(self.angle + delta, self.transmission)

       
    def rotateBy(self,delta):
        """
        Method to increment  angle of a polariser, (overload of general method in JonesMatrix)

        :param delta: Angle increment in radians
        :type delta: float
        """
        self.setAngle(self.angle + delta)

    


class Retarder(JonesMatrix):
    """
    Class to implemate a general retarder secified by phase shift and orientation of slow axis

    :param phase: Phase differebnce beween fast and slow axis in radians 
    :type phase: float
    :param theta: angle of slow axis
    :type theta: float or str
    :param transmission: Intensity transmission (Default = 1.0)
    :type transmission: float
    """

    def __init__(self, phase, theta = 0.0, transmission = 1.0):
        self.transmission = float(transmission)
        self.phase = float(phase)
        self.setAngle(theta)

    def copy(self):
        """
        Method to make copy of the current retander

        :return: copy of current
        """
        return Retarder(self.phase,self.angle,self.transmission)

                                  
    def setAngle(self,theta = None):
        """
        Methods to set the angle of the retarder

        :param theta: Angle of fast axis
        :type float:
        """
        if theta != None:
            self.angle = parseAngle(theta)
        amp = math.sqrt(self.transmission)
        val = cmath.rect(amp,0.5*self.phase)  # amp*exp(-i phase/2)

        a = val                               # X compenent
        b = complex(0,0)
        c = complex(0,0)
        d = val.conjugate()                   # Y component
        self.set(a,b,c,d)
        self.rotateBy(self.angle) # Rotate to specified angle


class WavePlate(Retarder):
    """
    Class to hold a waveplate 

    :param delta: phase thickness of the plate is wavelengths (Default = 0.5)
    :param order: order of the plate (Default = 0)
    :param theta: angle of the slow axis
    :param design: the design wavelength (Default = w.Design)
    """
    def __init__(self,delta = 0.5, order = 0, theta = 0.0, design = getDesignWavelength(),transmission = 1.0):
        """
        """
        self.delta = float(delta)
        self.order = int(order)
        self.design = float(design)
        self.transmission = float(transmission)
        self.angle = parseAngle(theta)
        self.setWavelength(getDefaultWavelength())
        

    def setWavelength(self,wave = getDefaultWavelength()):
        """
        Set the wavelength
        """
        if isinstance(wave,JonesVector):
            self.wavelength = wave.wavelength
        else:
            self.wavelength = float(wave)
        self.phase = 2*math.pi*(self.order + self.delta)*self.design/self.wavelength
        self.setAngle(self.angle)


class HalfWavePlate(WavePlate):
    """
    Class for a half wave retarded
    
    :param theta: angle of fast axis wrt x in radians, (Default = 0.0)
    :type theta: float
    :param order: order of halfwave plate, (Defaault = 0)
    :type order: int
    :param design: design wavellength, Defaults = wl.Design)
    :type design: float
    :param transmission: Transmission (Default =  1.0)
    :type transmittance: float
    
    """
    def __init__(self,theta = 0.0,order = 0, design = getDesignWavelength(),transmission = 1.0):
        WavePlate.__init__(self,0.5,theta,order,design,transmission)


    
        

class QuarterWavePlate(WavePlate):
    """
    Class for a quarter wave retarded
    
    :param theta: angle of fast axis wrt x in radians, (Default = 0.0)
    :type theta: float
    :param order: order of halfwave plate, (Defaault = 0)
    :type order: int
    :param design: design wavellength, Defaults = wl.Design)
    :param transmission: Transmission (Default =  1.0)
    :type transmission: float
    
    """
    def __init__(self,theta =  0.0,order = 0, design = getDesignWavelength(), transmission = 1.0):
        WavePlate.__init__(self,0.25,theta,order,design,transmission)




class JonesMatrixSystem(list):
    """
    Class to hold a list of JonesMatrix compoents
    
    :param ms: matrix or list of matrices
    :type ms: JonesMatrix or list of JonesMatrices
    """
    def __init__(self,ms=None,*args):
        list.__init__(self)
        if isinstance(ms,list) :        # Supplied a list
            for m in ms:
                self.append(d)          # Appled each

        if ms != None:                # Append single supplied Group
            self.append(ms)

        for a in args:                # Append additional groups (if any)
            self.append(a)

   
    def rotateBy(self,index,delta):
        """
        Method to rotate specifed component mnatrix by a specified angle.

        :param index: the component to rotate
        :type index: int
        :param delta: the rotation angle in radians
        :type delta: float
        """
        self[index].rotateBy(delta)

    def setWavelength(self,wave = getDefaultWavelength()):
        """
        Set the wavelength for all components
        """
        for m in self:
            m.setWavelength(wave)

    
    def getMatrix(self):
        """
        Method to form JonesMatrix by multiplying the current compoents together.

        :return: a single JoneMatrix
        """
        s = JonesMatrix()       # Start with unit matrix
        for m in self:
            s *= m

        return s

   
    def polarPlot(self, beam, index, points = 200):
        """
        Method to generate polar my rotating the one specided specfied component.
        The plot will be the beam intensity agaist angle
        
        :param beam: the input beam (this is not modified)
        :type beam: JonesVector
        :param index: the component to be rotated
        :type index: int
        :param key: plot key passes to plt.polar (Default = 'r')
        :type key: str
        :param points: Number of point on polar plot (Default = 200)
        :type points: int
        """
        if isinstance(beam,list):   # deal with a list
            for b in beam:
                self.polarPlot(b,index,points)
                
        else:
            self.setWavelength(beam)      # Set the wavelength once for whole rotation.
            theta = np.linspace(0.0,2*math.pi,num=points,endpoint=True)
            intensity = np.empty(theta.size)
        
            for i,angle in enumerate(theta):
                self[index].setAngle(angle)
                outputBeam = beam * self
                intensity[i] = outputBeam.getIntensity()
                

            col = WavelengthColour(beam.wavelength) # Get the colour.
            polar(theta,intensity,col)
