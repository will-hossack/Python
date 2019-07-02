"""
Set of classes to hold rays for optical ray tracing. This includes Paraxial and Intsity rays.

Author: Will Hossack, The Univesrity of Edinburgh
"""
import math
from vector import Vector3d,Vector2d,Unit3d,Angle
from optics.wavelength import Default,Spectrum,AirIndex,WavelengthColour
from optics.matrix import ParaxialMatrix,ParaxialGroup,ParaxialPlane
from matplotlib.pyplot import plot
#                
#



class Ray(object):
    """
    Base Ray class which just hold wavelength, intensity and refractive index all other (useful) 
    classes extend this Base class. This class has internal variables of

    - self.intensity (float) the intensity
    - self.refractiveindex the current Refactive Index
    - self.monitor the RayMonitor to track ray progress.

    :param wavelength: wavelength in microns (defaults to package Default (0.55um))
    :type wavelength: float
    :param intensity: float or Spectrum (defaults to 1.0)
    :type intensity: float or Spectrum
    :param index: refrctive index (Default = None)
    :type index: RefrativeIndex

    """
    
    def __init__(self,wavelength = Default, intensity = 1.0, index = None):
        """
        Constuctor with two optional arguments
        """
        self.wavelength = float(wavelength)
        if isinstance(intensity,Spectrum):
            self.intensity = intensity.getValue(self.wavelength)
        else:
            self.intensity = float(intensity)
    
        self.refractiveindex = index
        self.monitor = None                       # Hook for ray monitor

    def __str__(self):
        """
        Implement srt() to give basic imformation, typically overloaded by extending class.
        """
        return "l: {0:7.4e} i: {1:7.54} n: {2:s}".format(self.wavelength,self.intensity,repr(self.refractiveindex))
        
    #
    def __repr__(self):
        """
        Implement repr() to give more detailed information by including class name. 
        """
        return "{0:s} ".format(self.__class__.__name__) + str(self)
    #
    # 
    def setInvalid(self):
        """
        Method to set Ray as inValid, need to be defined in extending classes.
        """
        print("Ray.selInvalid needs to be defines")
    #
    #
    #
    def addMonitor(self,mon = None):
        """
        Method to add a RayMonitor to the ray, if called with None (the default) , it will 
        switch off the monitor. 

        :param mon: the ray monitor, for example RayPath() with will alow for plotting.
        :type mon: RayMonitor

        """
        self.monitor = mon
        return self.updateMonitor()

    def updateMonitor(self):
        """
        Update the monitor if it is exits. This is typicvally called auntomatically by the ray when its position is updated.
        """
        if self.monitor != None:
            self.monitor.update(self)
        return self
    
    def isValid(self):
        """
        Method to test if Ray is valid, needs to be defined in extending classes.

        :rerturn: must be True / False

        """
        print("Ray.isValid needs to be defined")
        

    def __bool__(self):
        """
        Implement logical __bool__ test of Ray is valid
        """
        return self.isValid()
    #
    def __iadd__(self,d):
        """
        Implement _iadd_ to propagate a distance d so implements rray += d
        """
        self.propagate(d)
        return self
    
    def __add__(self,d):
        """
        Impement __add__ to  implement newray = self + d
        """
        r = self.copy()
        r.propagate(d)
        return r

    
    def __radd__(self,d):
        """
        Impement __add__ to  implement ner newray = d + self
        """
        r = self.copy()
        r.propagate(d)
        return r
    

    def __imul__(self,s):
        """
        Implement ___imul__ to multiply by a surface / or matrix This the main ray tracing interface by impementing 
        ray *= s where s is an optical surface or for paraxial rays, a paraxial matrix.
        """
        self.propagateThrough(s)
        return self
        

    def __mul__(self,s):
        """
        Implement ___mul__ to multiply by a suraface
        """        
        r = self.copy()
        r.propagateThrough(s)
        return r

    
    def __rmul__(self,s):
        """
        Implement ___mul__ to multiply by a suraface
        """
        r = self.copy()
        r.propagateThrough(s)
        return r


    def draw(self):
        """
        Implement a draw the ray if there is a a suitable RayPath minotor attacked.
        """
        if self.monitor != None and isinstance(self.monitor,RayPath):
            self.monitor.draw()

#    
class ParaxialRay(Ray):
    """
    Class to implement Paraxial Rays with height and angle.

    :param height: height from optical axis (Default = 0.0)
    :type height: float
    :param angle:  angle in radians from optical axis (Defaults = 0.0)
    :type angle: float
    :param plane: location of plane along optical axis (Default = 0.0)
    :type plane: float
    :param wavelength: wavelength in microns defined in (Default = optics.wavelength.Default)
    :type wavelength: float
    :param intensity: intensity of ray (Default = 1.0)
    :type intensity: float

    """
    #    
    def __init__(self,height = 0.0, angle = 0.0, plane = 0.0, wavelength = Default ,intensity = 1.0 ):
        """
        Constuctor with 5 optional arguments
        """
        Ray.__init__(self,wavelength,intensity)   # Set wavelength and intensity
        self.h = float(height)                    # Ray height               
        self.u = float(angle)                     # Paraxial angle
        self.z = float(plane)                     # Position along optical axis

    
    def __str__(self):
        """
        Implment str() print parameters
        """
        return "h: {0:7.4f} u: {1:7.4f} p:  {2:7.4f} {3:s}".format(self.h,self.u,self.z,Ray.__str__(self))
    #            
    def copy(self):
        """
        Method to make a deep copy of a current ParaxialRay

        :return: ParaxialRay being copy of current ray

        """        
        return ParaxialRay(self.u,self.u,self.z,self.wavelength,self.intensity)
    #           
    def setInvalid(self) :
        """
        Method to set an ParaxialRay to inValid by settin the angle to "nan".
        """
        self.u = float("nan")       # set angle to be NaN 
        return self

    #           
    def isValid(self):
        """
        Method to test of a Paraxial Ray is valid, so test if angle is a "nan"
        
        :return: True if valid, else False

        """
        return not math.isnan(self.u)

    #            
    def propagate(self, distance):
        """
        Method to propagate a ray a specified distance with testing for validity. 
        This is a direct calculation ands not use matrix.

        :param distance: the distance the ray is propagated.
        :type distance: float
        :return: True is sucessful, False is Ray is invalid

        The Monitor.update will be called to record update the ray position.

        """
        if self:
            self.z += distance                # update plane
            self.h += self.u*distance         # calculate new height
            self.updateMonitor()
            return True
        else:
            return False

    
    def propagateTo(self, plane):
        """
        Method to propagate the rays to a specified position along the optical axis. 
        Will fail is the plane is infinite or ray is invalid.

        :param plane: the location of plane
        :type plane: float
        :return: True is sucessful, False if Ray is invalid

        """
        if  not math.isinf(plane):
            distance = plane - self.z         # distance to propagate
            return self.propagate(distance)       # do the propagatation
        else:
            return False                          # plane at inf, so fail

    
    def mult(self,m):
        """
        Method to multiply ParaxialRay by ParaxialMatrix and retuen a new array. 

        :param m: the ParxialMatrix
        :type m: optics.matrix.ParaxialMatrix
        :return: ParaxialRay

        Noramlly called via "*" operator
        """
        if self:
            h = self.h*m.A + self.u*m.B
            a = self.h*m.C + self.u*m.D
            p = self.z + m.thickness
            return ParaxialRay(h,a,p,self.wavelength,self.intensity)
        else:
            return self            
    
    
    def multBy(self,m):
        """
        Method to multiply ParaxialRay by ParaxialMatrix in place. This also automatically calls updateMonitor() is active.

        :param m: the ParxialMatrix
        :type m: optics.matrix.ParaxialMatrix
        :return: the updates ParaxialRay

        Normally called via \*= operator.

        """
        if self:
            h = self.h*m.A + self.u*m.B
            a = self.h*m.C + self.u*m.D
            self.h = h
            self.u = a
            self.z += m.thickness
            self.updateMonitor()
            return True
        else:
            return False

    #            
    #            
    def propagateThrough(self, surface):
        """
        Method to propagate Paraxial Ray through a surface or list of surfaces.
        this the main paraxial ray tracing method with same call as for the skew IntensityRays
        
        :param surface: ParaxialMatrix, ParaxialGroup and list of these
        :type surface: optics.matrix.ParaxialMatrix or optics.maytrix.ParaxialGroup
        :return: the updated ray

        Normally called by \*= operator.

        """

        if isinstance(surface,list):
            for s in surface:                # process each surface in the list in turn
                b = self.propagateThrough(s)
                if not b:                    # failed of surface, don't do anymore
                    break
                
            return b

            
        if isinstance(surface,ParaxialGroup):
            self.propagateTo(surface.inputPlane())
            if self.h > surface.inputPlaneHeight:
                self.setInvalid()
                return False
            else:
                return self.multBy(surface)

        if isinstance(surface,ParaxialMatrix):
            return self.multBy(surface)

        if isinstance(surface,float):
            return self.propagateToPlane(surface)
            

        #
        #           get interaction info as a list
        #
        info = surface.getParaxialInteraction(self)

        if math.isnan(info[1]):           # Distance failed
            return False                  # Exit now
           
        self.z += info[1]                 # Update position and height
        self.h = info[2]
        
        if self.monitor != None:          # Update monitor if it exits
            self.monitor.update(self)

        #          Check if curvature is valid
        if math.isnan(info[3]):
            self.setInvalid()                   # Set ray is invalid              
            return False                        # exit false

        #          Deal with different surface types
        #
        if info[0] == 0 :                            # Clear surface (do nothing)
            return True
        elif surface.type == 1:                      # Refratcion
            nl = self.refractiveindex.getValue(self) # Current refractive index
            nr = info[4].getValue(self)              # Refractive index after surface
                                                     # Calcualte new angle using paraxial relation
            self.u = self.h*info[3]*(nl - nr)/nr + nl*self.u/nr
            self.refractiveindex = info[4]     #     Update ray refrective index 
            return True                        #     success
        elif surface.type == 2:                #     Reflection
            self.u = 2.0*self.h*info[3] - self.u     # Calcualte relfecltion angle, will be -ve
            return True
        else:
            print("Intnesity Ray surface of unknown type implemented")

        return False                                   # If here we have failed (somehow), return false

    #
    #           
    def crossesZero(self):
        """
        Method to locate where ray crosses optical axis is global coordinates

        :return: float where ray cosses optical axis

        """
        if self.h == 0.0 :
            return self.plane            # already there
        elif self.u == 0:
            return float('inf')          # Trap infinity
        else:
            return self.z - self.h/self.u

    #         
    def crosses(self,other):
        """
        Method to locate where two ParaxialRays cross in global coordinates

        :param other: the other ParaxialRay
        :type other: ParaxialRay
        :return: float, position where they cross

        """
        dtheta = other.u - self.u
        if dtheta == 0.0:
            return float('inf')        #   Trap infinity where rays are parallel
        else:
            return (self.h - other.h - self.z*self.u + other.z*other.u)/dtheta
    


class SourcePoint(Vector3d):
    """
    Class implement a source point being a Position wih an attched intensity or spectrum.
    """
    #
    def __init__(self,pos,s_or_i = 1.0):
        """
        Create a 3d source point
        param pos, Position, the three-D Posistion
        param spectrum the attched spectrum
        """
        Vector3d.__init__(self,pos)
        self.spectrum = None                    # Add null to allow for testing
        if isinstance(s_or_i,Spectrum):      # Add spectrum if given
            self.spectrum = s_or_i
        else:
            self.intensity = float(s_or_i)      # else record intensity as a float


    def __str__(self):
        """
        Implement str
        """
        if self.spectrum == None:
            return Vector3d.__str__(self) + " i: " + str(self.intensity)
        else:
            return Vector3d.__str__(self) + " s: " + str(self.spectrum)

    def __repr__(self):
        """
        Implment repr()
        """
        return "ray.SoucePoint" + str(self)


    def copy(self):
        """
        Make a copy of the SourcePoint
        """
        if self.spectrum == None:
            return SourcePoint(self,self.intensity)
        else:
            return SourcePoint(self,self.spectrum)


    def clone(self,pos,s_or_i = 1.0):
        """
        Implement the a clone.
        """
        if self.spectrum == None:
            return SourcePoint(self,self.intensity)
        else:
            return SourcePoint(self,self.spectrum)

    def getIntensity(self,wave = Default):
        """
        Get the intensity as specified wavelength
        """
        if self.spectrum == None:
            return self.intensity
        else:
            return self.spectrum.getValue(wave)
#          

       
#          
class IntensityRay(Ray):
    """
    Class to form a Intensity Ray full vector ray tracing. 

    :param pos:  the starting position of the ray, or ParaxialRay
    :type pos: vector.Vector3d or ParaxialRay
    :param dirn: the starting direction of the ray (defaults to (0,0,1))
    :type dirn: vector.Unit3d or vector.Angle
    :param wavelength: the wavelenth (defaults to Default)
    :type wavelength: float
    :param intensity: intensity (defaults = 1.0)
    :type intensity: float or optics.wavelength.Spectrum
    :param index: RefractiveIndex, (defaults to AirIndex())
    :type index: optics.wavelength.RefractiveIndex

    """
    
    def __init__(self, pos = 0.0, dirn = 0.0 , wavelength = Default, intensity = 1.0, index = AirIndex()):
        """
        Consrctructor for to set parameters
        
        """
        if isinstance(pos,ParaxialRay):
            Ray.__init__(self,pos.wavelength,pos.intensity,pos.refractiveindex)
            self.position = Vector3d(0.0,pos.h,pos.z)
            self.director = Unit3d(Angle(pos.u))
        else:
            Ray.__init__(self,wavelength,intensity,index)  # Set wavelnegth intensity and index in super
            if isinstance(pos,float) or isinstance(pos,int):
                self.position = Vector3d(0,0,pos)
            else:
                self.position = Vector3d(pos)                  # Make localcopy of Position and Director (since they get updated)
            if isinstance(dirn,float) or isinstance(dirn,int):
                self.director = Unit3d(Angle(dirn)) 
            else:
                self.director = Unit3d(dirn)
        self.pathlength = None                             # Set opl to none (not calculated)




    def __str__(self):    
        """
        Implement str() to give detailed report an all variables for checking.
        """
        return "p : {0:s} u: {1:s} ".format(repr(self.position),repr(self.director)) + Ray.__str__(self)

    
    def copy(self):
        """
        Return a (deep) copy of the current Intensity

        :return: Deep copy of current Ray.
        """
        r = IntensityRay(self.position,self.director,\
                         self.wavelength,self.intensity,\
                         self.refractiveindex.copy())
        r.pathlength = self.pathlength
        return r

    
    def setInvalid(self):
        """
        Method to set the ray to inValid (sets the Unit3d as invalid)
        
        :rerturn: self

        """
        self.director.setInvalid()       # Set director to be invalid
        return self
    
    def isValid(self):
        """
        Method to check the ray is Valid (checks the Director is valid)

        :return: bool, True / False

        """
        return self.director.isValid()
    
    def getPhaselength(self):
        """
        Method to get the phase length, being 2*pi*pathelength/wavelength
        If it is not being calcualted, then will return None.
        
        :return: float, or None

        """
        if self.pathlengh == None:
            return None
        else:
            return 2000.0*math.pi*self.pathlength/self.wavelength
           
    def propagate(self,distance):
        """
        Method to propagate the ray a specifed distance using its own current direction.
        This also upadtes the pathlength and is the main method to propgate rays.
        param distance float the distance to propagate the ray.

        :param distance: distance to be propagated
        :type distance: float

        Normally called via "+=" operator with a float.

        Method clecks if the ray is valid, and if sp propagates it.
        """

        if self :
            self.position += distance*self.director

            if self.pathlength != None:
                self.pathlength += distance*self.refractiveindex.getValue(self.wavelength)

            self.updateMonitor()            # Uupdated the monitor     
       

        return self

    #
    #
    def propagateThrough(self,surface):
        """
        Method to propagate a ray through surface, or list of ray surfaces.
        This is the main method that does most of the work, it first propagate the
        ray to the surface and then depending on the surface it will block, reflect 
        or refract the ray depending on the type of surface. If the ray is refracted its 
        refractive index is also updated.

        :param surface:  Surface or list of list(Surface), if list each one is dealt with in order.
        :type surface: optics.surface.Surface or list of Surfces.
        :return: bool true is passed through, false if blocked. 

        Normalled called via the "\*=" operator.

        """
        #
        #      Deal with list.
        if isinstance(surface,list):
            for s in surface:                # process each surface in the list in turn
                b = self.propagateThrough(s)
                if not b:                    # failed of surface, don't do anymore
                    break
                
            return b

        #
        #           get interaction info as a list
        #
        info = surface.getSurfaceInteraction(self)
        if math.isnan(info.distance):           # Distance failed
            return False                  # Exit now
            
        self.position = info.position       # Update ray position
        
        if self.pathlength != None:       # Update pathlength if valid
            self.pathlength += info.distance*self.refractiveindex.getValue(self.wavelength)
        if self.monitor != None:          # Update monitor if it exits
            self.monitor.update(self)

        #          Check if surface normal is valid
        if not info.normal:
            self.setInvalid()                   # Set ray as invalid              
            return False                        # exit false

        #          Deal with different surface types
        #
        if info.type == 0 :                            # Clear surface (do nothing)
            return True

        elif info.type == 1:                                      # Refratcion
            nl = self.refractiveindex.getValue(self)              # Current refractive index
            nr = info.refractiveindex.getValue(self)              # Refractive index after surface
            ratio = nr/nl
            b = self.director.refraction(info.normal,ratio)       # Do the refratcion
            if b:
                self.refractiveindex = info.refractiveindex       #  Update ray refrective index 
                return True                        
            else:
                return False

        elif info.type == 2:                                      # Reflection
            return self.director.reflection(info.normal)          # Do reflection

        else:            
            raise TypeError("IntensityRay from unknow surface type {0:d}".format(info.type))

        return False                                   # If here we have failed (somehow), return false

    
    def pointInPlane(self,plane):
        """
        Method to calcualte where this ray striked an optical place
        This does NOT alter othe current ray.
        
        :param plane: the OpticalPlane or z the location on the optical axis
        :return: vector.Vector2d, the point in the plane.

        """
        if isinstance(plane,float):
            pt = Vector3d(0.0,0.0,plane)
        else:
            pt = plane.getPoint()
        if self:
            d = (pt.z - self.position.z)/self.director.z
            return Vector2d(self.position.x + d*self.director.x - pt.x,\
                            self.position.y + d*self.director.y - pt.y)
        else:
            return Vector2d().setInvalid()

        
        
#

class RayMonitor(object):
    """
    Class to monitor the progress of rays during the tracing process. The extending classes are used to record paths
    for printing / drawing etx.
    
    :param wavelength: Wavelength of ray, (Default = optics.wavelength.Default)
    :type wavelength: float

    """
    def __init__(self,wavelength = Default):
        """
        Create a defaut monitor that just hold wavelnegth
        """
        self.wavelength = wavelength

    def __repr__(self):
        """
        More detail with class name
        """
        return "{0:s} ".format(self.__class__.__name__) + str(self)

class PrintPath(RayMonitor):
    """
    Class to print changes of a ray math in real time.
    
    :param wavelength: Wavelength of ray, (Default = optics.wavelength.Default)
    :type wavelength: float
    
    """
    def __init__(self,ray = None, wavelength = Default):
        RayMonitor.__init__(self,wavelength)
        if ray != None:                              # If ray given add the monitor to the ray
            ray.addMonitor(self)
            self.wavelength = ray.wavelength

        
    def update(self,ray):
        """
        udate method,  prints the current ray position. Called automatically when the ray is propagated.

        :param ray: the current ray.

        """
        print(repr(ray))

class RayPath(RayMonitor):
    """
    Class to record a ray path. the path in held in three lists x[], y[] and z[]

    :param wavelength: Wavelength of ray, (Default = optics.wavelength.Default)
    :type wavelength: float

    """
    def __init__(self,ray = None, wavelength = Default):
        RayMonitor.__init__(self,wavelength)
        self.x = []
        self.y = []
        self.z = []
        if ray != None:                              # If ray given add the monitor to the ray
            ray.addMonitor(self)
            self.wavelength = ray.wavelength
        
    def copy(self):
        """
        Method to form a copy of itself
        """
        return RayPath()

    def __str__(self):
        """
        srt, return wavelength and num number of points in the path
        """
        return "l: {0:7.4f} n: {1:d}".format(self.wavelength,len(self.x))
    

    def getInfo(self):
        """
        Get the path formated as a string
        """
        s = repr(self)
        for i in range(len(self.x)):
            s += "\nx: {0:7.4f} y: {1:7.4f} z: {2:7.4f}".format(self.x[i],self.y[i],self.z[i])
        return s

    def update(self,ray):
        """
        udate method,  records the x,y,z position of the  ray is lists Called automatically as the ray is propagated.
        """
        self.wavelength = ray.wavelength
        if isinstance(ray,ParaxialRay):
            self.x.append(0.0)
            self.y.append(ray.h)
            self.z.append(ray.z)
        else:
            self.x.append(ray.position.x)
            self.y.append(ray.position.y)
            self.z.append(ray.position.z)
        

    def draw(self):
        """
        Plot to current axis with colour of ray given by its wavelength defined in wavelength.WavelengthColour.

        """
        col = WavelengthColour(self.wavelength)
        plot(self.z,self.y,col)
    

#
class RayPencil(list):
    """
    Class to hold a list of rays and implement methods to propagate on-mass; this will work for any
    type of Ray; at present IntensityRay and ParaxialRays (of a  mixture of the two).

    :param \*args: list of rays to be added to the Pancil (Default = None)

    """
    
    def __init__(self, *args):
        """
        Make a ray pencil with optional set of rays to be appended.
        """
        list.__init__(self)

        for r in args:
            self.append(r)
    
    def addCollimatedBeam(self,ca,u,key = "vl" ,nrays = 10 ,wave = Default, intensity = 1.0):
        """
        Method to add a collimated beam of IntensityRays

        :param ca: circular aperture to filled
        :type ca: optics.surface.CircularAperture
        :param u: direction of rays (can be Unit3d, Angle or float)
        :type u: vector.Unit3d or vector.Angle
        :param key: method of fill, allowed keys as "vl", "hl" and "array",(default is "vl")
        :type key: str
        :param nrays: number or rays across radius, (default = 10)
        :type nrays: int
        :param wave: the wavelength, (default = Default)
        :type wave: float
        :param intensity: the ray intensity, (default = 1.0)
        :type intensity: float or optics.wavelenth.Spectrum

        """

        if not hasattr(ca, "maxRadius"):
            ca = ca.entranceAperture()
        pt = ca.getPoint()         # Reference point
        radius = ca.maxRadius
        dr = radius/(nrays + 0.1)
        if isinstance(u,float) or isinstance(u,int):
            u = Unit3d(Angle(u))
        else:
            u = Unit3d(u)

        
        jmin = 0                  # Set default to central ray only
        jmax = 1
        imin = 0
        imax = 1
        
        if key == "vl":           # Vertical 
            jmin = -nrays
            jmax = nrays + 1
        elif key == "hl":
            imin = -nrays
            imax = nrays + 1
        elif key == "array":
            jmin = -nrays
            jmax = nrays + 1
            imin = jmin
            imax = jmax
        else:
            print("ray.RayPencil.addCollimatedBeam: illegal key {0:s}".format(str(key)))

        for j in range(jmin,jmax):
            for i in range(imin,imax):
                y = j*dr
                x = i*dr
                if x*x + y*y <= radius*radius:
                    dist = radius + x*u.x + y*u.y
                    p = Vector3d(pt.x + x, pt.y + y, pt.z)     # Point in aperture
                    p -= dist*u                                # Propagate back
                    ray = IntensityRay(p,u,wave,intensity)
                    self.append(ray)
        
        return self


    def addCollimatedParaxialBeam(self,ca,u,nrays = 10 ,wave = Default, intensity = 1.0):
        """
        Method to add a collimated paraxial beam
        param ca  aperture to fill
        param u direction of rays
        param nrays, number or rays aross radius, (default = 10)
        param wave, the wavelength, (default = Default)
        param intensity, the ray intensity, (default = 1.0)
        """

        if hasattr(ca, "maxRadius"):
            if isinstance(ca.maxRadius,float):
                radius = ca.maxRadius
            else:
                radius = ca.maxRadius()
        else:
            radius = 10.0
        pt = ca.getPoint().z
        dr = radius/(nrays + 0.1)

        jmin = -nrays
        jmax = nrays + 1

        for j in range(jmin,jmax):
            y = j*dr
            dist = radius + y*u
            ray = ParaxialRay(y,u,pt,wave,intensity)
            ray.propagate(-dist)
            self.append(ray)
        
        return self


    #
    def addSourceBeam(self, ca, source, key = "vl" ,nrays = 10 ,wave = Default):
        """
        Method to add beam from a source
        param ca circular aperture to fill
        param u direction of rays
        param key method of fill, allowed keys as "vl", "hl" and "array",(default is "vl")
        param nrays, number or rays aross radius, (default = 10)
        param wave, the wavelength, (default = Default)
        """
        if not hasattr(ca, "maxRadius"):
            ca = ca.entranceAperture()
        pt = ca.getPoint()         # Reference point
        radius = ca.maxRadius
        dr = radius/(nrays + 0.1)
        s = Vector3d(source)

        if hasattr(source,"intensity"):
            intensity = source.intensity
        else:
            intensity = 1.0  
        
        jmin = 0                  # Set default to central ray only
        jmax = 1
        imin = 0
        imax = 1
        
        if key == "vl":           # Vertical 
            jmin = -nrays
            jmax = nrays + 1
        elif key == "hl":
            imin = -nrays
            imax = nrays + 1
        elif key == "array":
            jmin = -nrays
            jmax = nrays + 1
            imin = jmin
            imax = jmax
        else:
            print("ray.RayPencil.addSourceBeam: illegal key {0:s}".format(str(key)))

        for j in range(jmin,jmax):
            for i in range(imin,imax):
                y = j*dr
                x = i*dr
                if x*x + y*y <= radius*radius:
                    p = Vector3d(pt.x + x, pt.y + y, pt.z)     # Point in aperture
                    u = Unit3d(p - s)
                    ray = IntensityRay(s,u,wave,intensity)
                    self.append(ray)
        
        return self


    def addSourceParaxialBeam(self,pg, height, sourceplane, nrays = 10 ,wave = Default, intensity = 1.0):
        """
        Add a ray paraxial ray pencil from a rounce to input of a paraxial group
        """
        if isinstance(sourceplane,ParaxialPlane):
            z = sourceplane.inputPlane()
        else:
            z = float(sourceplane)

        dist = pg.inputPlane() - z

        jmin = -nrays
        jmax = nrays + 1
        dy = pg.maxRadius()/(nrays+0.1)

        for j in range(jmin,jmax):
            y = dy * j                # Height in input aperture
            u = (y - height)/dist     # Calcualte angle using paraxial approx
            ray = ParaxialRay(height,u,z,wave,intensity)
            self.append(ray)
            
        return self
    #
    #
    def removeInvalid(self):
        """
        Method to remove invalid rays from the pencil
        """
        for r in self:
            if not r :
                self.remove(r)

    #
    def propagate(self,distance):
        """
        Method to propagate all rays a equals distance
        """
        for r in self:
            r.propagate(distance)
        return self
    #
    #
    def propagateThrough(self,sur):
        """
        Propagate the whole pencil hrough the surface
        """
        for r in self:             # For each ray
            if r:                  # in  each ray is valid
                r.propagateThrough(sur)

        return self
    #
    #
    #
    def __imul__(self,surface):
        """
        Implement ___rmul__ to multiply by a suraface
        """
        return self.propagateThrough(surface)
    #
    def addMonitor(self,monitor):
        """
        Method to add/remove a copy of the monitor to each ray
        """
        for r in self:
            if monitor == None:
                r.addMonitor()
            else:
                r.addMonitor(monitor.copy())
        return self
    #
    #
    def draw(self):
        """
        Draw each ray in turn.
        """
        for r in self:
            r.draw()
        
#               


        
#

class GaussianBeam(Ray):
    """
    Class to work with Gaussian Beams, class uses the underlying Ray class
    """
    def __init__(waist, wave = Default, intensity = 1.0):
        """
        Construct Gaussan beam
        param waist   beam waist in mm
        param wave wavelength in micoms
        param intesity the intensity
        """
        ParaxialRay.__init__(self,0.0,0.0,0.0,wave,intensity)     # Paraxial Ray class
        self.beam = complex(0.0,-math.pi*waist*waist*1000/self.wavelength)

    def getCurvature(self):
        """
        The beam curvature
        """
        return self.beam.real

    def getRadius(self):
        """
        Get the current beam radius
        """
        a = abs(self.beam.imag)
        return math.sqrt(self.wavelength/(math.pi*1000*a))

    def getWaist(self):
        """
        Get the minumum beam waist
        """
        c = self.getCurvature()
        a = abs(self.beam.imag)
        g = a/(c*c + a*a)
        return math.sqrt(self.wavelength*g/(math.pi*1000))

    def getWaistLocation(self):
        """
        Get the waist location in global coordinates
        """
        c = self.getCurvature()
        a = self.beam.imag
        z = -c*(c*c + a*a)
        return self.plane + z
        

    def getDivergence(self):
        """
        Get the divergence angle
        """
        w = self.getWaist()
        return 2.0*self.wavelength/(math.pi*1000*w)

    
        


        




