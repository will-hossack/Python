"""
Set of classes to deal with optical wavelengths and functions associated with 
wavelengths including refractive index. 
It aslo handles the default wavelength of the package.

This is part of the optics package.

Author: Will Hossack, The Univesrity of Edinburgh.
"""
import math
from os import getenv
from matplotlib.pyplot import plot
from vector import Vector2d,Vector3d
#
"""
Definitions of various wavelengths either used in the packages or tracing.
All values in microns.
"""
#                   Basic colours
Blue = 0.460
Green = 0.55
Red = 0.65
#                   Visual limits
BlueLimit = 0.35
RedLimit = 0.7
#                   Colour matching wavelnegth
BlueColourMatch = 0.425
GreenColourMatch = 0.53
RedColourMatch = 0.65
#                    Visual peaks
ScotopicPeak = 0.502819
ScotopicWidth = 0.0555317
PhotopicPeak = 0.5559087
PhotopicWidth = 0.0599179
#                    Spectral wavelengths
Mercury_i = 0.36501
Mercury_h = 0.4046561
Mercury_g = 0.4358343
Cadmium_F = 0.4799914
Hydrogen_F = 0.4861327
Mercury_e = 0.546073
Helium_d = 0.5875618
Sodium_D2 = 0.5889953
Sodium_D = 0.5892938
Sodium_D1 = 0.5895923
Cadmium_C = 0.6438469
Hydrogen_C = 0.6562725
Helium_r = 0.7065188
Potassium_A = 0.7682
Caesium_s = 0.85211
Mercury_t = 1.01398
#
#                       Laser wavelengths
ArgonBlue = 0.4880
ArgonGreen = 0.5145
KryptonRed = 0.6471
HeNeGreen = 0.5435
HeNeYellow = 0.594
HeNeRed = 0.6328
HeCdBlue = 0.441563
HeCdUV = 0.325
NdYagIR = 1.064
NdYagGreen = 0.532
Ruby = 0.6943
DiodeRed = 0.68
DiodeNearIR = 0.785
DiodeMidIR = 0.86
DiodeLongIR = 1.5
#
#
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

Default = getDefaultWavelength()

def setDefaultWavelength(w):
    """
    Function to set the default wavelength from within the package my resetting the global variable Default.
    """
    global Default
    Default = float(w)


FixedAirIndex = False
FixedAirIndexValue = 1.0

def setFixedAirIndex(type,value = None):
    """
    Set the Refrective index for the package, the default is variale, 
    """
    global FixedAirIndex
    global FixedAirIndexValue
    FixedAirIndex = bool(type)
    if FixedAirIndex and value != None:
        FixedAirIndexValue = float(value)
    

#
class WaveLength(object):
    """
    Define Abstract class to deal with functions of wavelength. This class also support
    plotting of the class via a .draw call.

    There are three local variable to control plotting and information
    self.minWavelength = BlueLimit
    self.maxWavelenth = Redlimit
    self.plotPoints = 200
    self.title = None
    """
    #
    #          
    def __init__(self):
        """
        Default constructor to set defaults, typically called by extending classes only
        No-parameters, it just initialises variables.
        """
        self.currentWavelength = float("Nan")     # Default to illegal
        self.currentValue = float("Nan")
        self.dynamic = False                      # Flag to force dynamics calls
        self.minWavelength = BlueLimit            # Plot paramters
        self.maxWavelength = RedLimit
        self.plotPoints = 200
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
        return "{0:s} ".format(self.__class__) + str(self)

    #         
    def getValue(self,wave = Default ):
        """
        Method to get the current value as the specified wavelength
        param wave the wavelength, this is assumes be a float OR any object
        that has .wavelength as a float variable. Defaults to Default wavelength for the package.

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
            self.currentValue = self.getNewValue(self.currentWavelength)
            
        return self.currentValue


    def getDerivative(self,wave):
        """
        Get the derivative dn/dl numerically using 4 point approximation with delta = wave / 2000 (which will have
        negligible errors for a smooth function)
        """
        if isinstance(wave,float) or isinstance(wave,int):
            w = wave
        elif isinstance(wave.wavelength,float):
            w = wave.wavelength
        else:
            raise TypeError("Wavelength.getDerivative(): call with unknown type {0:s}".format(str(wave)))
        
        delta = w  / 2000.0
        
        return (self.getNewValue(w - 2.0*delta) - 8.0*self.getNewValue(w - delta) + 8.0*self.getNewValue(w + delta) - \
           self.getNewValue(w + 2.0*delta))/(12.0*delta)
    
    #
    #         
    #
    def getNewValue(self,wave):
        """
        Abstract method to get the value at a new wavelength (needs to be defined)
        param wave float, the wavelength 
        return float the value at this wavelength
        """
        raise NotImplementedError("wavelength.Wavelength.getNewValue not implemnted, Class is abstract.")

    #         
    #
    def draw(self,key='r', derivative = False):
        """
        Method to get a matlibplot plot between self.minWavelenth and self.maxWavelength
        param key string matlibplot plotting key, (defaults to 'r')
        return matlibpolot "plot"
        """
        x = []
        y = []
        delta = (self.maxWavelength - self.minWavelength)/(self.plotPoints - 1)
        wave = self.minWavelength
        while wave <= self.maxWavelength:
            if derivative:
                val = self.getDerivative(wave)
            else:
                val = self.getValue(wave)
            x.append(wave)
            y.append(val)
            wave += delta

        if self.title == None:
            return plot(x,y,key)
        else:
            return plot(x,y,key,label=self.title)


#
#
class RefractiveIndex(WaveLength):
    """
    Class RefrativeIndex which extends WaveLength to handle differnt types of 
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
        """
        return self.getValue(Helium_d)
    #
    #          
    def getNe(self):
        """
        Method to get the refrative index at the Mercury_e line
        """
        return self.getValue(Mercury_e)
    #          
    #
    def getVd(self):
        """
        Method to get the Abbe or Vd number, calculated at the Helium_d, Hydroden_F and Hydroden_C lines.
        Will return zero if non-dispersive.
        return the Vd number as float
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
        return the Ve as float
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
        return 6 digit int
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
    """
    def __init__(self,formula,wrange,coef,name = "InfoIndex"):
        """
        Constuctor
        param formula int, formula type 5 
        param wrange, list of two float giving validity range
        param coef, list of floats holding the coefficents.
        param name or key (defaults to InfoIndex)
        
        Format is same as in RefrativeIndex.info database.
        """
        RefractiveIndex.__init__(self)
        self.formula = int(formula)
        self.R = list(wrange)
        self.C = list(coef)
        self.title = str(name)
    #
    #
    def __str__(self):
        """
        Implement st() to give full inpormation
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
    
    def getNewValue(self,wave):
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
    
        else:
            raise NotImplementedError("wavelength.InfoIndex: frommula {0:d} not impmented".\
                                      format(self.formula))



class FixedIndex(RefractiveIndex):
    """
    Class to implement a simple fixed index.
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
        return "wavelength.FixedIndex({0:7.5f})".format(self.value)
    
    def getNewValue(self,w = None):
        """
        Get the NewValue, always just the fixed value
        param w the wavelngth, which is ignored.
        """
        return self.value
    
        
#
class AirIndex(InfoIndex):
    """
    Class for AirIndex, this is either fixed or a special cals of InfoIndex with fixed paramers.
    Controlled by Global variiable
    FixedAirIndex which defaults to False.
    """
   
    def __init__(self):
        """
        No parameter conctructor.
        """
        InfoIndex.__init__(self,6,[0.23,1.69],[0, 0.05792105, 238.0185, 0.00167917, 57.362],"air")

    #     Set __str__ 
    def __str__(self):
        if FixedAirIndex:
            return "wavelength.AirIndex() fixed at {0:7.5f}".format(FixedAirIndexValue)
        else:
            return InfoIndex.__str__(self)

    #
    #            Make copy
    def copy(self):
        return AirIndex()
    #
    #            Method to get the new value at specified wavelength
    #
    def getNewValue(self,wave):
        if FixedAirIndex:
            return FixedAirIndexValue
        else:
            return InfoIndex.getNewValue(self,wave)   # Do the full calculation




class SimpleCauchyIndex(InfoIndex):
    """
    Class to implement a simple a + b/lambda^2 + c / lambda^4 Cauchy index
    """
    def __init__(self, a_or_nd, b_or_vd = None, c = None):
        """
              Implement a simple Cauchy index with variable parameter types
              parameter a,c,b three floats giving a + b/lambda^2 + c / lambda^4
              param nd,vd, two floats being the n_d and V_d values
              param typeint, int of format nnnVVV with nd = 1.nnn and Vd = VV.V
        """
        if c != None:                # all three parameters, 
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
            ip = a_or_nd/1000
            ap = a_or_nd - 1000*ip
            nd = 1.0 + ip/1000.0
            vd = ap/10.0
            SimpleCauchyIndex.__init__(self,nd,vd)
        else:
            raise TypeError("wavelnegth.SimpleCauchyIndex(): called with unknown types")

            
            


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
    Base Sepectrum class, implments a constant spectrum
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
        return "wavelength.Spectrum({0:7.5f})".format(self.brightness)

    #
    def getNewValue(self,wave):
        """
        Get the new value,also returns brighntess
        """
        return self.brightness

        

class GaussianSpectrum(Spectrum):
    """
    Class for a Guassian Spectrum
    """ 
    def __init__(self,peak,width,bright = 1.0):
        """
        param peak, the peak of the spectrum in microns.
        param width, the width to =/1 e{-1} point in microms
        param brigthness the peak brightness, defaults to 1.0
        """
        Spectrum.__init__(self,bright)
        self.peak = float(peak)
        self.width = float(width)
        self.title = "Gaussian Spectrum"

    def __str__(self):
        """
        The str function
        """
        return "({0:8.6e}, {1:8.6e}, {2:8.6e})".format(self.peak,self.width,self.brightness)

    
    def getNewValue(self,wave):
        """ Get the new value at specified wavelength
        param wave the wavelength in microns
        """
        d = wave - self.peak
        return self.brightness*math.exp(-(d*d)/(self.width*self.width))


class PlanckSpectrum(Spectrum):
    """
    Class to give the hot body Planck spectrum  at specified temperature and emissitivity
    """
    def __init__(self, t = 5000.0, emissitivity = 1.0):
        """
        param t the blackbody temperture in Kelvin
        param emissitivity the emmisttivity factor (defaults to 1.0)
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
        param t the tenperture"
        """
        self.t = t

    def getNewValue(self,wave):
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
    The Phototic (high light level) normal spectral response of the eye
    """
    
    def __init__(self,bright = 1.0):
        """
        param brightness, the brightness or peak value
        """
        GaussianSpectrum.__init__(self,PhotopicPeak,PhotopicWidth,bright)
        self.title = "Photopic Spectrum"

    
    def __str__(self):
        """
        The str() function
        """
        return "({0:8.6e})".format(self.brightness) 

#
class ScotopicSpectrum(GaussianSpectrum):
    """
    The Scotopic (dark adapted) specral response of the eye
    """
    def __init__(self,bright = 1.0):
        """
        param brightness, the brightness or peak value
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
        return "({0:8.6e})".format(self.brightness) 

#             Tricolour spectrum 
#
class TriColourSpectrum(Spectrum):

    #
    #               Constructor
    #
    def __init__(self,red = 1.0 ,green = 1.0 ,blue = 1.0 ,bright = 1.0,width = 0.03):
        Spectrum.__init__(self,bright)
        self.red = float(red)
        self.green = float(green)
        self.blue = float(blue)
        self.width = float(width)
        self.title = "Tricoloured Spectrum"

    #
    #          getNextValue
    #
    def getNewValue(self,wave):
        d = wave - Red
        value = self.brightness*self.red*math.exp(-(d*d)/(self.width*self.width))
        d = wave - Green
        value += self.brightness*self.green*math.exp(-(d*d)/(self.width*self.width))
        d = wave - Blue
        value += self.brightness*self.blue*math.exp(-(d*d)/(self.width*self.width))
        return value


class WavelengthColour(list):
    """
    Class to form as RGB list of floats to represent a wavelength colour.
    Based on 
    <a href="http://www.cox-internet.com/ast305/color.html">fortran code</a> 
    by Dan Bruton, Stephen F Austin State University.
    """
    
    #
    def __init__(self,wave):
        """
        parar wave, float is microns, 
        """
        list.__init__(self)
        self += [0.0,0.0,0.0]      # Default to black

        if wave > 0.37 and wave < 0.75:     # there is colour 
            
            #         Take linear multi-point dog-leg
            if wave < 0.44:
                self[0] = (0.44 - wave)/(0.44 - 0.37)
                self[2] = 1.0
            elif wave < 0.49 :
                self[1] = (wave - 0.44)/(0.49 - 0.44)
                self[2] = 1.0
            elif wave < 0.51:
                self[1] = 1.0
                self[2] = (0.51 - wave)/(0.51 - 0.49)
            elif wave < 0.58:
                self[0] = (wave - 0.51)/(0.58 - 0.51)
                self[1] = 1.0
            elif wave < 0.645 :
                self[0] = 1.0
                self[1] = (0.645 - wave)/(0.645 - 0.58)
            else:
                self[0] = 1.0
            
            #     Now correct for eye

            gamma = 0.7                    # gamma of eye
            d = wave - 0.56;               # Distance from peak of vision
            scale = 1.0 - d*d/0.03610944;  # Parabola with zeros 0.37 & 0.75

            self[0] = math.pow(scale*self[0],gamma)     #Scale by gamma (fit to monitor)
            self[1] = math.pow(scale*self[1],gamma)
            self[2] = math.pow(scale*self[2],gamma)

    #
    def hexString(self):
        """
        Method to return as a HTML hex string in #rrggbb where rr / gg / bb are the colours in Hex
        """
        red = int(round(self[0]*255))            # Scale in int 0 -> 255
        green = int(round(self[1]*255))
        blue = int(round(self[2]*255))
        
        #       Do a format 
        return "#{0:02X}{1:02X}{2:02X}".format(red,green,blue)







        

