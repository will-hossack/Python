"""
Set of classes to implement Polarsiatisaton analysis via Jones methods
"""
import math
import cmath
import wavelength as wl
from opticalray import Ray
from matplotlib.pyplot import polar,show

#
#            
#
class JonesVector(Ray):
    """
    JonesVector to hold polarised state of light
    """
    #
    #
    def __init__(self,x = 1.0 ,y = 0.0 ,wavelength = wl.Default):
        """
        Constructorcutor to create a JonesVector
        param x complex component of jones vector (defaults to 1.0)
        param y component of jones vector (defults to 0.0)
        wavelength of Jones vector (defaults to Default)
        also accept single argument of current JonesVector
        """
        if isinstance(x,JonesVector):
            self.__init__(x.x,x.y,x.wavelength)
        else:
            Ray.__init__(self,wavelength)       # Set underlying wavelength
            self.x = complex(x)
            self.y = complex(y)

    #
    #
    def __str__(self):
        """
        Implment the str() method
        """
        return "({0:s} , {1:s} , {2:7.5f} )".format(str(self.x),\
                                str(self.y),self.wavelength)
    #
    #
    def __repr__(self):
        """
        Implement the repr() method
        """
        return "jones.JonesVector" + str(self)

    #       
    def copy(self):
        """
        Method to make copy and a new JonesVector being a copy of the current.
        """
        return JonesVector(self)

    #
    #
    def __add__(self,b):
        """
        Implement the c = self + b where b and c must be both JonesVectors
        """
        if isinstance(b,JonesVector):        # Only sensible add
            return JonesVector(self.x + b.x, self.y + b.y, self.wavelength)
        else:
            raise TypeError("JonesVector _add_ incorrect argument type")

    ##
    #
    def __iadd__(self,b):
        """
        implment self += b where b must be a JonesVector
        """
        if isinstance(p,JonesVector):
            self.x += b.x
            self.y += b.y
        else:
            raise TypeError("JonesVector __iadd__ : incorrect argument type")
        return self
        
    #
    #
    def __sub__(self,b):
        """
        Implement the c = self + b where b and c must be both JonesVectors
        """
        if isinstance(b,JonesVector):        # Only sensible add
            return JonesVector(self.x - b.x, self.y - b.y, self.wavelength)
        else:
            raise TypeError("JonesVector __sub__ : incorrect argument type")

    #
    #
    def __isub__(self,b):
        """
        implment self -= b where b must be a JonesVector
        """
        if isinstance(p,JonesVector):
            self.x -= b.x
            self.y -= b.y
        else:
            raise TypeError("JonesVector __iadd__ : incorrect argument type")
        return self

    #    
    def getIntensity(self):
        """
        Method to get the intensity of the current JonesVector.
        returns intensity as a float.
        """
        v = self.x.real**2 + self.x.imag**2 + self.y.real**2 + self.y.imag**2
        return float(v)
    #
    #     
    def getPhase(self):
        """
        Method to get the phase difference between the two componets of 
        the current JonesVector as a float in the range -pi -> pi
        returns phase as a float.
        """
        delta = cmath.phase(self.y) - cmath.phase(self.x)
        if delta > math.pi :
            delta -= 2*math.pi
        if delta < -math.pi:
            delta += 2*math.pi
        return delta

    #    
    def getAngle(self):
        """
        Method to get the angle of the polarsiation ellipse of the
        current JonesVector in the range given by atan2
        """
        xm = abs(self.x)
        ym = abs(self.y)
        pa = self.getPhase()
        return 0.5*math.atan2(2.0*xm*ym*math.cos(pa) , (xm*xm - ym*ym))

    #
    #
    def getEllipicity(self):
        """
         Method to get the ellipticity of the polarisation ellipse
        of the current JonesVector as a float on range -pi/2 -> pi/2
        """
        xm = abs(self.x)
        ym = abs(self.y)
        pa = self.getPhase()
        return 0.5*math.asin(2*xm*ym*math.sin(pa)/(xm*xm + ym*ym))

    #    Method to multiply by the current JonesVector by a JonesMatrix 
    #    and return new JonesVector
    def mult(self,ms):
        r = self.copy()
        if isinstance(ms,list):
            for m in ms:
                r.multBy(m)
        else:
           r.multBy(ms)
        return r

    #   
    def __mul__(self,m):
        """
        Implment the __mul__ tio pre-multiply by a JoneMatrix
        """
        if isinstance(m,JonesMatrix):
            p = self.x*m.A + self.y*m.B
            q = self.x*m.C + self.y*m.D
            return JonesVector(p,q,self.wavelength)
        else:
            raise TypeError("JonesVector.multBy: argument not a JonesMatrix")
            
         
    def __imul__(self,m):
        """
        Implment the __imul__ to pre-multiply by a JoneMatrix
        """
        if isinstance(m,JonesMatrix):
            p = self.x*m.A + self.y*m.B
            q = self.x*m.C + self.y*m.D
            self.x = p
            self.y = q
            return self
        else:
            raise TypeError("JonesVector.multBy: argument not a JonesMatrix")
            

    #    Method to multiply the current JonesVector by a JonesMatrix
    #    in place.
    def multBy(self,ms):
        if isinstance(ms,list) :
            for m in ms:
                self.multBy(m)
        elif isinstance(ms,JonesMatrix) :
            p = self.x*ms.A + self.y*ms.B
            q = self.x*ms.C + self.y*ms.D
            self.x = complex(p)
            self.y = complex(q)
        else:
            raise TypeError("JonesVector.multBy: argument not a JonesMatrix")

    #     Method to put the current JonesVector through an idea
    #     Linear Polarsier as specfied angle and get the output intensity
    #     The current JoneVector is not altered.
    def throughPolariser(self,theta):
        pol = LinearPolariser(theta)
        r = self * pol
        return r.getIntensity()


    #     Method to generate polar plot through a rotated polariser
    def polarPlot(self,key='r',points = 200):
        theta = [0.0]*(points + 1)
        intensity = [0.0]*(points + 1)
        delta = 360.0/points

        for i in range(points + 1):
            angle = math.radians(i*delta)
            theta[i] = angle
            intensity[i] = self.throughPolariser(angle)

        return polar(theta,intensity,key)
        

#
#       Class to make JonesVector for linear polarsied light
class LinearPolarisedBeam(JonesVector):

    #      Define constuctor
    def __init__(self,theta = 0.0, intensity = 1.0, wavelength = wl.Default):
        amp = math.sqrt(intensity)
        JonesVector.__init__(self,amp*math.cos(theta),amp*math.sin(theta), \
                             wavelength)

#
#        Class to make a JoneVector for Right Circular Polarsied light
class RightCircularPolarisedBeam(JonesVector):

    #
    #     Define constructor
    def __init__(self,intensity = 1.0, wavelength = wl.Default):
        amp = math.sqrt(intensity/2.0)
        JonesVector.__init__(self,amp,complex(0.0,-amp),wavelength)

#        Class to make a JoneVector for Right Circular Polarsied light
class LeftCircularPolarisedBeam(JonesVector):

    #
    #     Define constructor
    def __init__(self,intensity = 1.0, wavelength = wl.Default):
        amp = math.sqrt(intensity/2.0)
        JonesVector.__init__(self,amp,complex(0.0,amp),wavelength)

#
#          Class to definbe a JonesMatrix to represent a polarisation
#          component.
class JonesMatrix(object):

    #       Define constructor with all 4 complex components
    #       defaults to identity matrix, with optional angle
    def __init__(self,a=1.0,b=0.0,c=0.0,d=1.0):
        self.set(a,b,c,d)       # Set the compenents
    #
    #       Internal set method to set 4 components, all complex
    #       a_or_jm componet
    #       b component
    #       c component
    #       d component
    #       also if called with a_or_jm = JoneMatrix 
    def set(self,a_or_jm,b=None,c=None,d=None):
        if isinstance(a_or_jm,JonesMatrix):
            self.A = a_or_jm.A
            self.B = a_or_jm.B
            self.C = a_or_jm.C
            self.D = a_or_jm.D
        else:
            self.A = complex(a_or_jm)
            self.B = complex(b)
            self.C = complex(c)
            self.D = complex(d)

    #     Form the determinant of Matrix 
    #     return complex, the complex determinant
    def determinant(self) :
        det = self.A*self.D - self.B*self.C
        return det
    
    #       Form the compex trace of the matrix
    #       return the complex trace
    def trace(self):
        tr = self.A + self.D
        return tr

    #       Method to copy the current JonesMatrix
    def copy(self):
        return JonesMatrix(self.A,self.B,self.C,self.D)

    #      Method to multiply the current JonesMatrix and return new JonesMatrix
    #      m JonesMatrix to multiply by
    #      returns new JonesMatrix and leave current unchanged.
    def mult(self,m):
        n = Jonesmatrix.copy(self)
        n.multBy(m)
        return n

    #     Overload * method to multiply two JonesMatricces and rerturn 
    #     new JonesMatrix
    def __mul__(self,m):
        return self.mult(m)

    #      Method to multiply current JonesMatrix in place
    #      m JonesMatrix to muluiply by
    #      this overwrited the current values
    def multBy(self,m):
        if isinstance(m,JonesMatrix):
            a = m.A*self.A + m.B*self.C
            b = m.A*self.B + m.B*self.D
            c = m.C*self.A + m.D*self.C
            d = m.C*self.B + m.D*self.D
            self.set(a,b,c,d)  # Undate current
        else:
            raise TypeError("JonesMatrix.multBy: argument not a JonesMatrix.")

    #      Method to rotate a general JonesMatix by angle and return new matrix
    #      using the pre and post rotation matrices.
    def rotate(self,angle):
        cos = math.cos(angle)
        sin = math.sin(angle)
        #           First rotation matrix
        rp = JonesMatrix(cos,-sin,sin,cos)
        a = self.mult(rp)
        ra = JonesMatrix(cos,sin,-sin,cos)
        return ra.mult(a)

    #     Method to rotate the current JonesMatrix by specified angle in place.
    def rotateBy(self,angle):
        cos = math.cos(angle)
        sin = math.sin(angle)
        #           First rotation matrix
        rp = JonesMatrix(cos,-sin,sin,cos)
        self.multBy(rp)
        ra = JonesMatrix(cos,sin,-sin,cos)
        ra.multBy(self)
        self.set(ra)        # set self with rotated values.

    #      String method
    def __str__(self):
        return "JonesMatrix: [ " + str(self.A) + "  " + str(self.B) + "\n" + \
            str(self.C) + "  " + str(self.D) + "]"

    

#         Class to make a linear polariser 
class LinearPolariser(JonesMatrix):

    #      Constructor for a liner polariser
    #      theta angle of polarsier wrt x axis (default 0.0)
    #      transmission intensity tramission of polarsier (default to 1.0)
    #
    #      Note angle is in radians
    def __init__(self,theta = 0.0,transmission = 1.0):
        self.transmission = transmission
        self.setAngle(theta)

    #      Method to return a copy
    def copy(self):
        return LinearPolariser(self.angle,self.transmission)

    #     Method to set the angle of the current Linear Polarsier
    #     theta angle of the polarsier wrt to x-axis
    def setAngle(self,theta):
        self.angle = theta      #   Record the angle
        cos = math.cos(theta)
        sin = math.sin(theta)
        amp = math.sqrt(self.transmission)
        self.set(amp*cos*cos,amp*cos*sin,amp*cos*sin,amp*sin*sin)

    #    Method to increment the the angle of a polariser, (overload)
    #    of method in JonesMatrix for effiency.
    #    delta angle to be incremented
    def rotateBy(self,delta):
        self.setAngle(self.angle + delta)

#       Class to make a general retarder
class Retarder(JonesMatrix):

    #    Full constructor
    def __init__(self, phase, theta, transmission = 1.0):
        self.transmission = transmission
        self.phase = phase
        self.setAngle(theta)

    #    Method to take a copy
    def copy(self):
        return Retarder(self.phase,self.angle,self.transmission)
    #
    #    Method to set the angle of retarder
    def setAngle(self,theta):
        self.angle = theta
        
        amp = math.sqrt(self.transmission)
        val = cmath.rect(amp,0.5*self.phase)  # amp*exp(-i phase/2)

        a = val                               # X compenent
        b = complex(0,0)
        c = complex(0,0)
        d = val.conjugate()                   # Y component
        self.set(a,b,c,d)
        JonesMatrix.rotateBy(self,self.angle) # Rotate to specified angle
        

    #    Method to increment the the angle the retarder , (overload)
    def rotateBy(self,delta):
        theta = self.angle + delta           # keep angle updated
        self.setAngle(theta)


#       Class for HalfWaveRetarder
class HalfWavePlate(Retarder):
    
    #           General constructor for a half wave retarded
    #           theta angle of fast axis wrt x, (defaults to zero)
    #           order order of halfwave plate, (defaults to zero)
    #           wavelength wavelength of halfwave plate (defaults to Default)
    #           transmittance (defaults to 1.0)
    #           design wavellength, defaults to Default
    def __init__(self,theta = 0.0,order = 0.0, wavelength = wl.Default,\
                 transmittance = 1.0, design = wl.Default):
        phase = 2.0*math.pi*(order + 0.5)*design/wavelength
        Retarder.__init__(self,phase,theta,transmittance)


#                Ideal Half wave plate for testing
class IdealHalfWavePlate(Retarder):
    #
    #          Constructor that just take an angle
    def __init__(self,theta = 0.0):
        self.setAngle(theta)

    #          method to set the plate at angle using exact formula
    #
    def setAngle(self,theta):
        self.angle = theta
        cos = math.cos(2.0*theta)
        sin = math.sin(2.0*theta)
        a = complex(cos)
        b = complex(sin)
        c = complex(sin)
        d = complex(-cos)
        self.set(a,b,c,d)


#       Class for QuanterWaveRetarder
class QuarterWavePlate(Retarder):
    
    #           General constrcuor for a half wave retarded
    def __init__(self,theta =  0.0,order = 0.0, wavelength = wl.Default,\
                 transmittance = 1.0, design = wl.Default):
        phase = 2.0*math.pi*(order + 0.25)*design/wavelength
        Retarder.__init__(self,phase,theta,transmittance)



#        Class to hold a list of JonesMatrix compoents
class JonesMatrixSystem(list):
    
    #
    #        Constructor to take a sereies of JonesMatrices
    #        either as a set of argumens or a list.
    #
    def __init__(self,ms=None,*args):
        list.__init__(self)
        if isinstance(ms,list) :        # Supplied a list
            for m in ms:
                self.append(d)          # Appled each

        if ms != None:                # Append single supplied Group
            self.append(ms)

        for a in args:                # Append additional groups (if any)
            self.append(a)

    #     Method to rotate the specified Matrix 
    def rotateBy(self,index,delta):
        self[index].rotateBy(delta)


    #     Method to form a single matrix from the list
    def getMatrix(self):
        s = JonesMatrix()
        for m in self:
            s.multBy(m)

        return s

    #
    #     Method to generate polar plot you rotating the
    #     specfied component
    def polarPlot(self, beam, index, key='r', points = 200):

        
        delta = math.radians(360.0/(points - 1))
        thetas = [0.0]*points
        intensity = [0.0]*points

        for i in range(points):
            theta = i*delta
            thetas[i] = theta
            self.rotateBy(index,delta)
            outputBeam = beam.mult(self)
            intensity[i] = outputBeam.getIntensity()

        return polar(thetas,intensity,key)
