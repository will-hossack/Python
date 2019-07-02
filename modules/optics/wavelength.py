"""
Set of classes to deal with optical wavelengths and functions associated with 
wavelengths including refractive index and spectra.  It aslo handles the default wavelength of the package.

"""
import math
from os import getenv
from matplotlib.pyplot import plot
import numpy as np
from vector import Vector2d,Vector3d
from optics.material import MaterialData, Material


Blue = 0.460        #: Blue in microns
Green = 0.55        #: Green in microns
Red = 0.65          #: Red in microns
BlueLimit = 0.35    #: Visual Blue limit
RedLimit = 0.7      #: Visual Read Limit
BlueColourMatch = 0.425 #: Colour matching Blue in microns
GreenColourMatch = 0.53 #: Colour matching Green in microns
RedColourMatch = 0.65   #: Colour matching Blue in microns
ScotopicPeak = 0.502819   #: Scotopic peak in microns
ScotopicWidth = 0.0555317 #: Scotopic width in microns 
PhotopicPeak = 0.5559087  #: Photopic peak in microns
PhotopicWidth = 0.0599179 #: Photopic with in microns
Mercury_i = 0.36501    #: Mercury i line in microns
Mercury_h = 0.4046561  #: Mercury h line in microns
Mercury_g = 0.4358343  #: Mercury g line in microns
Cadmium_F = 0.4799914  #: Cadium F line in microns.
Hydrogen_F = 0.4861327 #: Hydrodgen F line in microns
Mercury_e = 0.546073   #: Mercury e line in microns
Helium_d = 0.5875618   #: Helium d line
Sodium_D2 = 0.5889953  #: Sodium D2 line
Sodium_D = 0.5892938   #: Sodium D (average) line
Sodium_D1 = 0.5895923  #: Sodium D1 line
Cadmium_C = 0.6438469  #: Cadmium C line
Hydrogen_C = 0.6562725 #: Hydeogen C line
Helium_r = 0.7065188   #: Helium r line
Potassium_A = 0.7682   #: Potatium A line
Caesium_s = 0.85211    #: Caesium s line
Mercury_t = 1.01398    #: Mercury t line
ArgonBlue = 0.4880   #: Argon Ion Blue
ArgonGreen = 0.5145  #: Argon Ion Green
KryptonRed = 0.6471  #: Krypton Ion Red
HeNeGreen = 0.5435   #: Helium Neon Green
HeNeYellow = 0.594   #: Helium Neon Yellow
HeNeRed = 0.6328     #: Helium Neon Red
HeCdBlue = 0.441563  #: Helium Cadmium Blue
HeCdUV = 0.325       #: Helium Cadmium UV
NdYagIR = 1.064      #: Neodinium primary (IR)
NdYagGreen = 0.532   #: Neodinium doubled (Green)
Ruby = 0.6943        #: Ruby red.
DiodeRed = 0.68      #: Laser diode red
DiodeNearIR = 0.785  #: Laser diode near IR
DiodeMidIR = 0.86    #: Laser diode mid IR
DiodeLongIR = 1.5    #: Laser diode long IR


def getDefaultWavelength():
    """
    Function to extarct the default wavelength from environmentalvariable DEFAULTWAVELENGTH.
    It can be set to numerical value or any the wavelengths specified as Globals above.
    If environemntal variale not set, then defaults of Green = 0.55 is used

    This is called automatically once on startup and set the global variable Default.

    """
    val = getenv("DEFAULTWAVELENGTH")
    if val == None:
        return Green
    else:
        return float(eval(val))
#
#          Package default wavelength held in global.

Default = getDefaultWavelength()       #     Package default wevelength

def setDefaultWavelength(w):
    """
    Function to set the default wavelength from within the package my resetting the global variable Default.

    :param w: The new default wavelength in microns
    :type w: float

    """
    global Default
    Default = float(w)


FixedAirIndex = False
FixedAirIndexValue = 1.0

def setFixedAirIndex(type = True ,value = 1.0):
    """
    Set the Refrective index for the package, the default is variale.

    :param type: Sets fixed index for package, (defaults to True)
    :type type: bool
    :param value: the value of the fixed index (default to 1.0)
    :type value: float

    """
    global FixedAirIndex
    global FixedAirIndexValue
    FixedAirIndex = type
    FixedAirIndexValue = value
    

#
class WaveLength(object):
    """
    Define Abstract class to deal with functions of wavelength.  There is a default
    constructor to set defaults, typically called by extending classes only there are no parameters, 
    it just initialises variables. 

    There are three local variable to control plotting and information

    - self.dynamic = False (forces new calcculation for every call)
    - self.minWavelength = BlueLimit
    - self.maxWavelenth = Redlimit
    - self.plotPoints = 100
    - self.title = None


    """          
    def __init__(self):
    
        self.valid = True                         # Defaults to True
        self.currentWavelength = float("Nan")     # Default to illegal
        self.currentValue = float("Nan")
        self.dynamic = False                      # Flag to force dynamics calls
        self.minWavelength = BlueLimit            # Plot paramters
        self.maxWavelength = RedLimit
        self.plotPoints = 100
        self.title = None                         # User title

    def __str__(self):
        """
        The str functoon, will be overwritten
        """
        return "Abstract"

    def __repr__(self):
        """
        The repr function, may get overwritten
        """
        return "{0:s} ".format(self.__class__.__name__) + str(self)

    def __bool__(self):
        """ Validity flag  
        """
        return self.valid

    #         
    def getValue(self,wave = Default ):
        """
        Method to get the current value as the specified wavelength
        param wave the wavelength, this is assumes be a float OR any object
        that has .wavelength as a float variable. 

        :param wave: wavelength (Defaults to package default, usually 0.55 microns)
        :type wave: float or object.wavelength (in microns)
        :return: the value as a float.

        This will cache the value from the last call so if called again with the same wavelength it will not recalculate. 
        This is controlled by self.dynamic, which is set to True will calcualte a new value for every call.

        This is the normal call for all classes.
        """
        if isinstance(wave,float) or isinstance(wave,int):
            w = wave
        elif isinstance(wave.wavelength,float):
            w = wave.wavelength
        else:
            raise TypeError("Wavelength.getValue(): call with unknown type {0:s}".format(str(wave)))
        #
        #                           If value at this wavelength is know, use it
        #       
        if self.dynamic or w != self.currentWavelength:
            self.currentWavelength = w
            self.currentValue = self.__getNewValue__(self.currentWavelength)
            
        return self.currentValue


    def getDerivative(self,wave = Default):
        """
        Get the derivative dn/dl numerically using 4 point approximation with delta = wave / 2000 (which will have
        negligible errors for a smooth function)

        :param wave: wavelength (Defaults to package default, usually 0.55 microns)
        :type wave: float or object.wavelength (in microns)
        :return: the derivative as a float.

        """
        if isinstance(wave,float) or isinstance(wave,int):
            w = wave
        elif isinstance(wave.wavelength,float):
            w = wave.wavelength
        else:
            raise TypeError("Wavelength.getDerivative(): call with unknown type {0:s}".format(str(wave)))
        
        delta = w  / 2000.0
        
        return (self.getValue(w - 2.0*delta) - 8.0*self.getValue(w - delta) + 8.0*self.getValue(w + delta) - \
           self.getValue(w + 2.0*delta))/(12.0*delta)
           

    def __getNewValue__(self,wave):
        """
        Abstract iunternal method to get the value at a new wavelength (needs to be defined)

        :param wave: the wavelength 
        :type wave: float
        :return: the value at this wavelength.

        This is not normally called direclty but is called via .getValue()
        
        """
        raise NotImplementedError("wavelength.Wavelength.getNewValue not implemnted, Class is abstract.")

    def getArrayValues(self,wave_array):
        """
        Get an array of values for a NumPy array

        :param wave_array: a Numpy array of wavelengths in microms
        :type wave_array: np.array
        :return: values as a Numpy.array

        """
        out_array = np.empty(wave_array.size)
        for i,wave in enumerate(wave_array):
            out_array[i] = self.getValue(wave)

        return out_array

    def getArrayDerivatives(self,wave_array):
        """
        Get array of derivatives values for a NumPy array

        :param wave_array: a Numpy array of wavelengths in microms
        :type wave_array: Numpy.array
        :return: values as a Numpy.array

        """
        out_array = np.empty(wave_array.size)
        for i,wave in enumerate(wave_array):
            out_array[i] = self.getDerivative(wave)

        return out_array

    #         
    #
    def draw(self,colour='r', derivative = False):
        """
        Plot to the current Matplotlib axis.

        - self.minWavelenth (Default 0.35)
        - self.maxWavelength (Default 0.65)
        
        :param colour: colour passed to matplotlib (defaults to 'r')
        :type colour: Color or str
        :param derivative: if True will plot the derivative, (default is False)
        :type derivatice: Bool
        :return: None
        

        """
        x = np.linspace(self.minWavelength,self.maxWavelength,self.plotPoints)
        if derivative:
            y = self.getArrayDerivatives(x)
        else:
            y = self.getArrayValues(x)

        if self.title == None:
            plot(x,y,c=colour)
        else:
            plot(x,y,c=colour,label=self.title)


#
#
class RefractiveIndex(WaveLength):
    """
    Class RefrativeIndex which extends WaveLength to handle different types of 
    Refrative Index. Class is Abstract, need to be extended to be useful.
    """
    #
    #
    def __init__(self):
        """
        Default constuctor, same as WaveLength, does nothing but sets defaults.
        """
        WaveLength.__init__(self)
    #
    #         
    def getNd(self):
        """
        Method to get refrative index at the Helium_d line.

        :return: float value of referative index Nd

        """
        return self.getValue(Helium_d)
    #
    #          
    def getNe(self):
        """
        Method to get the refrative index at the Mercury_e line.

        :return:float calue of refratcive index Ne

        """
        return self.getValue(Mercury_e)
    #          
    #
    def getVd(self):
        """
        Method to get the Abbe or Vd number, calculated at the Helium_d, Hydroden_F and Hydroden_C lines.
        Will return zero if non-dispersive.

        :return: float the Vd number.

        """
        nd = self.getNd()
        nf = self.getValue(Hydrogen_F)
        nc = self.getValue(Hydrogen_C)
        if nf != nc:
            return (nd - 1.0)/(nf - nc)
        else:
            return 0.0               # Non dispersive
    #    
    #
    def getVe(self):
        """
        Method to get the Mercury Abbe or Ve number, cacaulted at the Mercury_e, Cadmium_F and Cadmium_C lines.
        Will return zero if non-dispersive. 

        :return: float  the Ve numbers.

        """
        ne = self.getNe()
        nf = self.getValue(Cadmium_F)
        nc = self.getValue(Cadmium_C)
        if nf != nc:
            return (ne - 1.0)/(nf - nc)
        else:
            return 0.0
    #
    #
    def getType(self):
        """
        Method to get the type number on nnnVVV for format using calculated from the Nd and Vd numbers.
        
        :return: int XXXYYY when nd = 1.XXX and Vd = Y.YY

        """
        nd = self.getNd()
        vd = self.getVd()
        nd = int(round((nd - 1.0)*1000))
        vd = int(round(vd*10))
        return nd*1000 + vd
        
#
#
class InfoIndex(RefractiveIndex):
    """
    Implment a RefractiveIndex in the format suppled by RefractiveIndex.info website,
    all have common calls.
    
    :param formula: formula type, (1,2,3,5,6 and 7 implementated)
    :type formula: int
    :param wrange: list of two float giving validity range
    :type wrange: list
    :param coef: list of floats holding the coefficents.
    :type coef: list
    :param name:  name or key (defaults to InfoIndex)
    :type name: str

    This class is not typically called by the used, the normal intreface is MaterialIndex
    
    """
    def __init__(self,formula,wrange,coef,name = "InfoIndex"):
        """
        Constuctor
        
        
        Format is same as in RefrativeIndex.info database.
        """
        RefractiveIndex.__init__(self)
        self.formula = int(formula)
        if self.formula == 0:
            self.valid = False
        self.R = list(wrange)
        self.C = list(coef)
        self.title = str(name)
    #
    #
    def __str__(self):
        """
        Implement str() to give full inpormation.
        """
        return "({0:d}, {1:s}, {2:s}, {3:s})".format(self.formula,str(self.R),str(self.C),self.title)
    #
    #
    def copy(self):
        """
        Define copy method to give deep copy
        """
        return InfoIndex(self.formula,list(self.R),list(self.C),self.title)
    #
    #          
    
    def __getNewValue__(self,wave):
        """
        Method to get the new value as specified wavelength
        param wave float the wavelength
        """
        if wave < self.R[0] or wave > self.R[1]:
            print("InfoIndex: range warning: {0:7.5f} called but range {1:7.5f} to {2:7.5f}".\
                  format(wave,self.R[0],self.R[1]))        

        #        Implment the various formula
        #
        if self.formula == 1:               # Sellmeir with Sqr of lower compoent
            n = 1.0 + self.C[0]             # Put in first C (nornmally zero)
            lSqrInv = 1.0/(wave*wave) 
            for i in range(1,len(self.C),2):
                c = self.C[i + 1]
                n += self.C[i]/(1.0 - lSqrInv*c*c)
            
            return math.sqrt(n)

        elif self.formula == 2:             # Sellmeir with Linear lower component.
            n = 1.0 + self.C[0]             # Put in first C (nornmally zero)
            lSqrInv = 1.0/(wave*wave) 
            for i in range(1,len(self.C),2):
                n += self.C[i]/(1.0 - lSqrInv*self.C[i + 1])
            
            return math.sqrt(n)

        elif self.formula == 3:            # Polynomial
            n = self.C[0]                  # Put in first C 
        
            for i in range(1,len(self.C),2):
                n += self.C[i]*pow(wave,self.C[i+1])
               
            return math.sqrt(n)


        elif self.formula == 5:             # Full Cauchy
            n = self.C[0]                   # Put in first C 
            for i in range(1,len(self.C),2):
                n += self.C[i]*pow(wave,self.C[i+1])
               
            return n

        elif self.formula == 6:            # Gas index
            n = 1.0 + self.C[0]            # First parmeter, usually 0
            overlSqr = 1.0/(wave*wave)
            for i in range(1,len(self.C),2):
                n += self.C[i]/(self.C[i + 1] - overlSqr)

            return n

        elif self.formula == 7:            # Herzberger
            lSqr = wave*wave
            n = self.C[0] + self.C[1]/(lSqr - 0.028) + self.C[2]*(1/(lSqr - 0.028)**2)
            for i in range(3,len(self.C)):
                lSqr *= (wave*wave)
                n += self.C[i]*lSqr
            return n
                                                            
            
    
        else:
            raise NotImplementedError("wavelength.InfoIndex: frommula {0:d} not impmented".\
                                      format(self.formula))        
    
class MaterialIndex(InfoIndex):
    """
    General class for Material Index with indices looked up in ethe internal package dattbase (is normal user interface)

    :param key: the index key, eg "BK7" 
    :type key: str
    :param database: the index database, (Default = None which use package standard)
    :type database: str

    
    """
    def __init__(self,key = None, database = None):
        """ Created a material refratcive index by key
        param key the Index Key (name), if None, you will be prompted
        parameter database Material database, of None, used package default
        """
        if  key != None and key.lower().strip() == "air":     # Trap air as special case
            AirIndex.__init__(self)
        else:
            material = MaterialData(database).getMaterial(key)     # Get material
            InfoIndex.__init__(self, material.formula,material.wrange,material.coef,material.name)


    def __str__(self):
        """           str to give name and range only
        """
        return " n: {0:s} r : ({1:7.5f} , {2:7.5f})".format(self.title, self.R[0], self.R[1])
            
        
class FixedIndex(RefractiveIndex):
    """
    Class to implement a simple fixed index.

    :param n: the refrative index
    :type n: float

    """
    def __init__(self,n):
        """
        Implement a simple fixed index
        param n the fixed index
        """
        RefractiveIndex.__init__(self)
        self.value = float(n)
        self.title = "Fixed Index"
        
    def __str__(self):
        """
        Implment str()
        """
        return " n: ({0:7.5f})".format(self.value)
    
    def __getNewValue__(self,w):
        """
        Get the NewValue, does not depend of wavelength.
        param w the wavelngth, which is ignored.
        """
        return self.value
    
        
class AirIndex(InfoIndex):
    """
    Class for AirIndex, this is either fixed or a special cals of InfoIndex with fixed paramers. Controlled by Global variiable
    self.FixedAirIndex which defaults to False.

    :param none: No parameters


    """
   
    def __init__(self):
        """
        No parameter conctructor.
        """
        InfoIndex.__init__(self,6,[0.23,1.69],[0, 0.05792105, 238.0185, 0.00167917, 57.362],"air")

    #     Set __str__ 
    def __str__(self):
        if FixedAirIndex:
            return " fixed: {0:7.5f}".format(FixedAirIndexValue)
        else:
            return " variable "

    def copy(self):
        """
        Make a copy
        """
        return AirIndex()
    
    def __getNewValue__(self,wave):
        """
        Get a new value as specified wavelenght taking into accouth the Fixed index flag
        """
        if FixedAirIndex:
            return FixedAirIndexValue
        else:
            return InfoIndex.getNewValue(self,wave)   # Do the full calculation

class CauchyIndex(InfoIndex):
    """
    Class to implement the simple a + b/lambda^2 + c / lambda^4 Cauchy index, with either a,b,c or  Nd,Vd or
    XXXYYY int whhere Nd = 1.XXX and Vd = Y.YY

    :param a_or_nd:  a vaule in the Cauchy formula if 3 parmeters, the Nd value if two, or type integer if one.
    :type a_or_nd: float OR int
    :param b_or_vd: b value in Cauchy formula if 3 parameters, else Vd id two given
    :type b_or_vd: float
    :param c: c value in Caucby formula (may be None)
    :type c: float

    """
    def __init__(self, a_or_nd, b_or_vd = None, c = None):
        """
              Implement a simple Cauchy index with variable parameter types
             
        """
        if c != None:                      # all three parameters given
            a = float(a_or_nd)
            b = float(b_or_vd)
            c = float(c)
            InfoIndex.__init__(self,5,[0.35,1.0][a,b,-2,c,-4],"Cauchy")    # Use InforIndex of type 5 with 3 paramters

        elif b_or_vd != None:             # Two paramters
            nd = float(a_or_nd)
            vd = float(b_or_vd)
            lf = Hydrogen_F
            lc = Hydrogen_C
            ld = Helium_d

            b = lf*lf*lc*lc*(nd - 1)/(vd*(lc*lc - lf*lf))
            a = nd - b/(ld*ld)
            InfoIndex.__init__(self,5,[0.35,1.0],[a,b,-2],"Cauchy")       # Use InfoIndex of type 5 with 2 parameters
        elif isinstance(a_or_nd,int):                            # int index/Vd type, unpack and re-call with two params
            ip = a_or_nd//1000
            ap = a_or_nd - 1000*ip
            nd = 1.0 + ip/1000.0
            vd = ap/10.0
            CauchyIndex.__init__(self,nd,vd)
        else:
            raise TypeError("wavelnegth.SimpleCauchyIndex(): called with unknown types")

    def __str__(self):
        """
        Implements str to give nd and Vd values.
        """
        return " nd: {0:7.5f} Vd: {1:7.4f}".format(self.getNd(),self.getVd())
            
            


class GradedIndex(RefractiveIndex):
    """
    Class to implement a graded index with a underlying base index and a radially symmeetric variation
    that depend on radial distance from an origin.
    """
    
    def __init__(self,pt,index,coef):
        """
        param pt, two dimensional point giving the location of the origin, 
        param index, the base refrative index
        param coef the for radial polynomial if form 1,r^2,r^4 .... as a list of floats.
        """
        self.point = Vector2d(pt)
        self.index = index
        self.coef = list(coef)


    def copy(self):
        """
        Make fully copy
        """
        return GradedIndex(self.point,self.index.copy(),self.coef)

    def __repr__(self):
        return "wavelength.GradedIndex: pt : {0:s} coef : {1:s}\n base : {2:s}".format(str(self.point),\
                                        str(self.coef),repr(self.index))


    def getValue(self, ray_or_wave):
        """
        Method to get the value, overloads method in super classes
        param ray_or_wave either ray or wavelength, if wavelength gets based value
        if ray get position dependand index.
        """
        base = self.index.getValue(ray_or_wave)
        if hasattr(ray_or_wave,"position"):    # Its a ray
            p = ray_or_wave.position
            dx = p.x - self.point.x
            dy = p.y - self.point.y
            rsqr = dx*dx + dy*dy
            r = 1.0
            weight = self.coef[0]
            for c in self.coef[1:]:
                r *= rsqr
                weight += c*r  
                return base*weight
        else:
            return base                       # its a scalar

        



    def getGradient(self,ray):
        """
        Method to calcualte the Gradient as specifed ray position
        """
        p = ray.position
            
        n = self.index.getValue(ray)
        
        x = p.x - self.point.x
        y = p.y - self.point.y
        rsqr = x*x + y*y

        dx = 2*n*x*self.coef[1]
        dy = 2*n*y*self.coef[1]
        """
        for i in range(2,len(self.coef)):
            r = math.pow(rsqr,i-1)
            dx += 2*x*n*i*coef[i]*r
            dy += 2*y*n*i*coef[i]*r
        """


        return Vector3d(dx,dy,0.0)

class Spectrum(WaveLength):
    """
    Base Sepectrum class, implments a constant spectrum.

    :param bright: the intensity brightness (Default = 1.0)
    :type bright: float

    """
    #  
    #
    def __init__(self,bright = 1.0):
        """
        Set only brighnness
        param bright, float, the brighness,(defaults to 1.0)
        """
        WaveLength.__init__(self)
        self.brightness = float(bright)

    def __str__(self):
        return " b: {0:8.4f}".format(self.brightness)

    #
    def __getNewValue__(self,wave):
        """
        Get the new value, always returns brighntess
        """
        return self.brightness

        

class GaussianSpectrum(Spectrum):
    """
    Class for a Guassian Spectrum with specifed peak, width and brightness.

    :param peak: the peak of the spectrum in microns.
    :type peak: float
    :param width: the width to =/1 e{-1} point in microms
    :type width: float
    :param bright: the peak brightness, (Defaults to 1.0)
    :type bright: float

    """ 
    def __init__(self,peak,width,bright = 1.0):
        """
       
        """
        Spectrum.__init__(self,bright)
        self.peak = float(peak)
        self.width = float(width)
        self.title = "Gaussian Spectrum"

    def __str__(self):
        """
        The str function
        """
        return "p: {0:7.5f} w: {1:7.5f} b: {2:7.5f}".format(self.peak,self.width,self.brightness)

    
    def __getNewValue__(self,wave):
        """ 
        Get the new value at specified wavelength
        param wave the wavelength in microns
        """
        d = wave - self.peak
        return self.brightness*math.exp(-(d*d)/(self.width*self.width))


class PlanckSpectrum(Spectrum):
    """
    Class to give the hot body Planck spectrum  at specified temperature and emissitivity

    :param t: the blackbody temperture in Kelvin (Default = 5,000)
    :type t: float
    :param emissitivity: the emmisttivity factor (defaults to 1.0)
    :type emissitivity: float
    
    """
    def __init__(self, t = 5000.0, emissitivity = 1.0):
        """
       
        """
        Spectrum.__init__(self)
        self.t = float(t)
        self.e = float(emissitivity)
        self.c1 = 5.026e6
        self.c2 = 1.44e4

    def __str__(self):
        """
        The str() function
        """
        return "({0:8.6e}, {1:8.6e})".format(self.t,self.e)


    def setTemperature(self,t):
        """  
        Set the tempertaure

        :param t: the tenperture in degrees Kelvin
        :type t: float

        """
        self.t = t

    def __getNewValue__(self,wave):
        """
        Get the new value 
        param wave the wavelength
        """
        val = self.c1*math.pow(1.0/wave,5)*(1.0/(math.exp(self.c2/(wave*self.t)) - 1.0))
        return val*self.e



#            PhotopicSpectrum 
#
class PhotopicSpectrum(GaussianSpectrum):
    """
    The Phototic (high light level) normal spectral response of the eye.

    :param bright: peak brightness, (Default = 1.0)
    :type bright: float


    """
    
    def __init__(self,bright = 1.0):
        """

        """
        GaussianSpectrum.__init__(self,PhotopicPeak,PhotopicWidth,bright)
        self.title = "Photopic Spectrum"

    
    def __str__(self):
        """
        The str() function
        """
        return "b: {0:7.4f}".format(self.brightness) 

#
class ScotopicSpectrum(GaussianSpectrum):
    """
    The Scotopic (dark adapted) specral response of the eye

    :param bright: peak brightness, (Default = 1.0)
    :type bright: float

    """
    def __init__(self,bright = 1.0):
        """

        """
        GaussianSpectrum.__init__(self,ScotopicPeak,ScotopicWidth,bright)
        self.title = "Scotopic Spectrum"

    #
    #        The __str__ method
    #
    def __str__(self):
        """
        The str() function
        """
        return "b: {0:7.4f}".format(self.brightness) 

#             Tricolour spectrum 
#
class TriColourSpectrum(Spectrum):
    """ 
    Implement a three coloured spectrum (RGB), all with same width

    :param red: the red peak, (default 1.0)
    :type red: float
    :param green: the green peak, (default 1.0)
    :type green: float
    :param blue: the blue peak, (default 1.0)
    :type blus: float
    :param bright: the overall brightness, (default 1.0)
    :type bright: float
    :param width: the width of each peak, (default 0.025)
    :type width: float

    """
    def __init__(self,red = 1.0 ,green = 1.0 ,blue = 1.0 ,bright = 1.0,width = 0.025):
        """ The tricoloured spectrum
        
        """
        Spectrum.__init__(self,bright)
        self.red = float(red)
        self.green = float(green)
        self.blue = float(blue)
        self.width = float(width)
        self.title = "Tricoloured Spectrum"

    def __str__(self):
        """
        str function
        """
        return "b : {0:6.3f} [{1:6.3f}, {2:6.3f}, {3:6.3f}] w: {4:6.3f}".\
            format(self.brightness,self.red,self.green,self.blue,self.width)
        

    def __getNewValue__(self,wave):
        """ Get the new value at specified wavelength
        """
        d = wave - Red
        value = self.brightness*self.red*math.exp(-(d*d)/(self.width*self.width))
        d = wave - Green
        value += self.brightness*self.green*math.exp(-(d*d)/(self.width*self.width))
        d = wave - Blue
        value += self.brightness*self.blue*math.exp(-(d*d)/(self.width*self.width))
        return value


def WavelengthColour(wave):
    """
    Class to form as RGB list of floats to represent a wavelength colour.
    Based on 
    <a href="http://www.cox-internet.com/ast305/color.html">fortran code</a> 
    by Dan Bruton, Stephen F Austin State University.

    :param wave: wavelnegth in microns
    :type wave: float
    :return: hexstring of the colour
 
    """
    
    
    rgb = [0.0,0.0,0.0]      # Default to black

    if wave > 0.37 and wave < 0.75:     # there is colour 
            
        #         Take linear multi-point dog-leg
        if wave < 0.44:
            rgb[0] = (0.44 - wave)/(0.44 - 0.37)
            rgb[2] = 1.0
        elif wave < 0.49 :
            rgb[1] = (wave - 0.44)/(0.49 - 0.44)
            rgb[2] = 1.0
        elif wave < 0.51:
            rgb[1] = 1.0
            rgb[2] = (0.51 - wave)/(0.51 - 0.49)
        elif wave < 0.58:
            rgb[0] = (wave - 0.51)/(0.58 - 0.51)
            rgb[1] = 1.0
        elif wave < 0.645 :
            rgb[0] = 1.0
            rgb[1] = (0.645 - wave)/(0.645 - 0.58)
        else:
            self[0] = 1.0
            
            #     Now correct for eye

        gamma = 0.7                    # gamma of eye
        d = wave - 0.56;               # Distance from peak of vision
        scale = 1.0 - d*d/0.03610944;  # Parabola with zeros 0.37 & 0.75

        rgb[0] = math.pow(scale*rgb[0],gamma)     #Scale by gamma (fit to monitor)
        rgb[1] = math.pow(scale*rgb[1],gamma)
        rgb[2] = math.pow(scale*rgb[2],gamma)

        red = int(round(rgb[0]*255))            # Scale in int 0 -> 255
        green = int(round(rgb[1]*255))
        blue = int(round(rgb[2]*255))
        
        #       Do a format 
        return "#{0:02X}{1:02X}{2:02X}".format(red,green,blue)
    


def RefractiveIndexColour(index = 1.5):
    """
    Class to form as RGB list of floats to represent a refratcive index colour for diagrams

    :param index: value of refrative index (Default = 1.5)
    :type index: float
    :return: hexstring with colour

    """
    if isinstance(index,RefractiveIndex):     # Allow for different parameters
        n = index.getValue()
    else:
        n = float(index)
            
    rgb =  [0.5,0.5,1.0]      # Default to grey/blue

    index_min = 1.4
    index_max = 2.3

    delta = (n - index_min)/(index_max - index_min)
    rgb[1] = 1.0 - delta*rgb[1]
    rgb[1] = min(1.0,max(0.0,rgb[1]))
    
    red = int(round(rgb[0]*255))            # Scale in int 0 -> 255
    green = int(round(rgb[1]*255))
    blue = int(round(rgb[2]*255))
        
        #       Do a format 
    return "#{0:02X}{1:02X}{2:02X}".format(red,green,blue)




        

