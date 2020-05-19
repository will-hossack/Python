"""
Set of classes to hold rays for optical ray tracing. This includes Paraxial and Intensity rays.
"""
import math
from vector import Vector3d,Vector2d,Unit3d,Angle
from optics.wavelength import Spectrum,AirIndex,WavelengthColour,getDefaultWavelength
from optics.matrix import ParaxialMatrix,ParaxialGroup,ParaxialPlane
from matplotlib.pyplot import plot



#   Global Current Angle (mainly used by GUI)
CurrentAngle = Unit3d(0.0,0.0,1.0)

def getCurrentAngle():
    """
    Get the current angle used in the package, typically called from GUI interface.
    The Default is (0,0,1) so along the optical axis.
    
    :return: The current angle global value
    """
    return CurrentAngle

def setCurrentAngle(u):
    """
    Method to set a current angle, typivcally used by the GUI intreface
    
    :param u: the current new current angle, any format accepted by Unit3d
    :type u: Unit3d or any format accepted or float

    """
    global CurrentAngle
    if isinstance(u,float):
        u = Angle(u)
    CurrentAngle = Unit3d(u)

def getCurrentSourcePoint(self):
    """
    Method to get the current SourcePoint, typically used by GUI intreface
    """
    return CurrentSource

def setCurrentSourcePoint(self,p):
    """
    Method to set the current SourcePoint, typically usde by GUI interface
    """
    global CurrentSource
    CurrentSource = p.copy()
    


class SourcePoint(Vector3d):
    """
    Class implement a source point being a Position with an attached fixed intensity or a spectrum.

    :param pos: the position of the source point
    :type: vector.Vector3d
    :param s: Spetrum or intensity (Default = 1.0)
    :type s_or_i: :class:`optics.wavelength.Spectrum` or float
    """
    #
    def __init__(self,pos,s = 1.0):
        """
        Create a 3d source point
        """
        if isinstance(pos,float) or isinstance(pos,int):
            Vector3d.__init__(self,[0,0,pos])
        else:
            Vector3d.__init__(self,pos)
        if isinstance(s,Spectrum):         # Add spectrum if given
            self.spectrum = s
        else:
            self.spectrum = Spectrum(s)    # Set constant spectrum


    def __str__(self):
        """
        Implement str
        """
        return Vector3d.__str__(self) + " s: " + str(self.spectrum)

    def __repr__(self):
        """
        Implment repr()
        """
        return "{0:s}".format(__class__.__name__) + str(self)


    def copy(self):
        """
        Make a copy of the SourcePoint

        :return: copy of current `SourcePoint`

        """
        return SourcePoint(self,self.spectrum)


    def getIntensity(self,wavelength = None):
        """
        Get the intensity as specified wavelength
        """
        if wavelength == None:
            wavelength = getDefaultWavelength()
        
        return self.spectrum.getValue(wavelength)


CurrentSource = SourcePoint(0.0)


class Disc(object):
    """
    Simple disc object for making simple beams without needing the complexities of
    optical surfaces.
    
    :param pt: the centre of the Disc is global coordinates (Default = 0,0,0) 
    :type pt: Vector3d or float
    :param radius: The radius (Default = 1.0)
    :type radius: float                                                      )
    """
    def __init__(self,pt = 0.0,radius = 1.0):
        
        if isinstance(pt,float) or isinstance(pt,int):
            self.point = Vector3d(0.0,0.0,pt)
        else:
            self.point = Vector3d(pt)
            
        self.maxRadius = float(radius)
         
    def getPoint(self):
        """
        Method to get the reference point
        :return: the reefrence point as Vector3d
        """
        return self.point
    
    def getRadius(self):
        """
        Method of get the radius
        :return: the radius as a float
        """
        return self.maxRadius
    

class Ray(object):
    """
    Base Ray class which just hold wavelength, intensity and refractive index all other (useful) 
    classes extend this Base class. This class has internal variables of

    - self.intensity the intensity as `float`
    - self.refractiveindex the current :class:`optics.wavelength.RefactiveIndex`
    - self.monitor the `RayMonitor` to record ray progress.

    :param wavelength: wavelength in microns (defaults to package Default = 0.55um)
    :type wavelength: float
    :param intensity: the intensity  (defaults to 1.0)
    :type intensity: `float` or :class:`optics.wavelength.Spectrum`
    :param index: refrative index (Default = None)
    :type index: :class:`optics.wavelength.RefractiveIndex`

    """
    
    def __init__(self,wavelength = None, intensity = 1.0, index = None):
        """
        Constuctor with two optional arguments
        """
        self.wavelength = getDefaultWavelength(wavelength)
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
    
    
    def getIntensity(self):
        """
        Get the intetsity 
        
        :return: the intensity as a float
        """
        return self.intensity
    
    def getWavelength(self):
        """
        Get the wavelength of the ray in microns
        
        :return: the wavelength as a float
        """
        return self.wavelength
    
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
        Update the monitor if it is exits. This is typically called auntomatically by the ray when its position is updated
        and not by the user.

        :return: self 
        """
        if self.monitor != None:
            self.monitor.update(self)
        return self
    
    def isValid(self):
        """
        Method to test if Ray is valid, needs to be defined in extending classes.

        :return: (bool) so True / False

        Also impletents __bool__

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
        Implement _iadd_ to propagate a distance d so implements ray += d
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
    def __init__(self,height = 0.0, angle = 0.0, plane = 0.0, wavelength = None ,intensity = 1.0 ):
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
    
       
#          
class IntensityRay(Ray):
    """
    Class to form a Intensity Ray full vector ray tracing. 

    :param pos:  the starting position of the ray, or ParaxialRay
    :type pos: vector.Vector3d , SourcePoint or ParaxialRay
    :param dirn: the starting direction of the ray (defaults to (0,0,1))
    :type dirn: vector.Unit3d or vector.Angle
    :param wavelength: the wavelenth (defaults to Default)
    :type wavelength: float
    :param intensity: intensity (defaults = 1.0) If SourcePoint given, intensity is calculated from wavelength.
    :type intensity: float or optics.wavelength.Spectrum
    :param index: RefractiveIndex, (defaults to AirIndex())
    :type index: optics.wavelength.RefractiveIndex

    """
    
    def __init__(self, pos = 0.0, dirn = 0.0 , wavelength = None, intensity = 1.0, index = AirIndex()):
        """
        Constructor for to set parameters
        
        """
        if isinstance(pos,ParaxialRay):
            Ray.__init__(self,pos.wavelength,pos.intensity,pos.refractiveindex)
            self.position = Vector3d(0.0,pos.h,pos.z)
            self.director = Unit3d(Angle(pos.u))
        else:
            if isinstance(pos,SourcePoint):
                intensity = pos.spectrum
            Ray.__init__(self,wavelength,intensity,index)     # Set wavelnegth intensity and index in super
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
    
    def getAngle(self):
        """
        Get the director an an Angle
        :return: Angle if valid, else None
        """
        if self.isValid() :
            return Angle(self.director)
        else:
            return None
    
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
        
            
        
    def rotateAboutX(self,angle,origin = None):
        """
        Roate the ray by specified andgle and Origin. If origin is None, then
        (0,0,0) is assumes
        
        :param angle: rotation angle in radians
        :type angle: float
        :param origin: Rotation origin (Default = None)
        :type origin: Vector3d or None
        """
        
        #            Deal with position 
        
        if self:             # Cleck the ray is valid
        
            if origin != None:
                self.position -= origin         # Move to orgin is givem
            self.position.rotateAboutX(angle)   # Rotate position
            if origin != None:
                self.position += origin         # Put origin back in
            self.director.rotateAboutX(angle)   # Rotate director
            self.updateMonitor()                # Updated the monitor     
        return self
        
        
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
            self.position += self.director*distance

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
        self.updateMonitor()                # Update the monitor
        
        if self.pathlength != None:       # Update pathlength if valid
            self.pathlength += info.distance*self.refractiveindex.getValue(self.wavelength)
       

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
        Method to calcualte where the ray will striked a specified optical plane
        This does NOT alter the current ray.
        
        :param plane: the OpticalPlane or z the location on the optical axis.
        :type plane: :class:`optics.surface.OpticalPlane` or float
        :return: :class:`vector.Vector2d`, the point in the plane relative to the plane reference point.
        """
        if isinstance(plane,float) or isinstance(plane,int):
            pt = Vector3d(0.0,0.0,plane)
        else:
            pt = plane.getPoint()
        if self:
            d = plane.getDistance(self.position,self.director)
            #(pt.z - self.position.z)/self.director.z
            return Vector2d(self.position.x + d*self.director.x - pt.x,\
                            self.position.y + d*self.director.y - pt.y)
        else:
            return Vector2d().setInvalid()

        
        
#

class RayMonitor(object):
    """
    Class to monitor the progress of rays during the tracing process. The extending classes are used to record paths
    for printing / drawing etc.
    
    :param wavelength: Wavelength of ray, (Default = None, take current default)
    :type wavelength: float

    """
    def __init__(self,wavelength = None):
        """
        Create a defaut monitor that just hold wavelnegth
        """
        self.wavelength = getDefaultWavelength(wavelength)

    def __repr__(self):
        """
        More detail with class name
        """
        return "{0:s} ".format(self.__class__.__name__) + str(self)

    def update(self,ray):
        """
        Update method, overloaded by inheriting classes
        """
        print("Ray Monitor update method not defined")

    def draw(self):
        """
        Blank method, overloaded in extending classes if used\
        """
        return None

class PrintPath(RayMonitor):
    """
    Class to print changes of a ray math in real time.
    
    :param wavelength: Wavelength of ray, (Default = optics.wavelength.Default)
    :type wavelength: float
    
    """
    def __init__(self,ray = None, wavelength = None):
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

    :param wavelength: Wavelength of ray, (Default = None, give package default)
    :type wavelength: float

    """
    def __init__(self,ray = None, wavelength = None):
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
    
    def add(self,ray):
        """
        Add for compatibility with othe classes, just appends thge ray to the pencil
        :param ray: the Ray to be appended
        :type ray: Ray
        
        """
        self.append(ray)
        
        
    def addRays(self,pt,u,wavelengths,intensity = 1.0):
        """
        Add a set of rays with specified wavelngths
        
        :param pt: the positions of each ray
        :param u: the direction of each ray
        :param wavelengths: the list or np.array of wavelengths
        :param intensity: the intensities, can also be list of same length aw wavelenths
        """
        
        for i,wave in enumerate(wavelengths):
            if isinstance(intensity,float) or isinstance(intensity,int):
                ival = float(intensity)
            else:   
                ival = intensity[i]
            ray = IntensityRay(pt,u,wave,ival)
            self.add(ray)
        
        return self
        

    def addBeam(self, ca, source, key = "vl", nrays = 10, wavelength = None, intensity = 1.0, index = AirIndex(), path = False):
        """
        Method to add a beam if intensity rays being either Collimated or Source Beam. The beam will fill the given circular aperture and will
        either come from a single SourcePoint or at a specified angle.

        :param ca: circular aperture to filled (any object with maxRadius attribute)
        :type ca: optics.surface.CircularAperture
        :param source: source or rays, either a SourcePoint or angle. 
        :type source: SourcePoint or vector.Unit3d or vector.Angle or float
        :param key: method of fill, allowed keys as "vl", "hl" and "array",(default is "vl")
        :type key: str
        :param nrays: number or rays across radius, (default = 10)
        :type nrays: int
        :param wavelenth: the wavelength, (default = Default)
        :type wavelength: float
        :param intensity: the ray intensity, (default = 1.0) only used for Collimated beam; for SourceBeam picked up from SourcePoint
        :type intensity: float
        :param index: the refratcive index, (Default = AirIndex())
        :type index: RefractiveIndex
        :param path: record pathlength, (default = False) is pathlength of each ray recorded
        :type path: bool
        :return: self

        """

        #          Sort out aperture to fill.
        if not hasattr(ca, "maxRadius"):
            ca = ca.entranceAperture()
        pt = ca.getPoint()         # Reference point
        radius = ca.maxRadius
        dr = radius/(nrays + 0.1)
        
        if isinstance(source,SourcePoint):        # Rays from a source
            s = Vector3d(source)
            intensity = source.getIntensity(wavelength)
        else:
            s = Vector3d().setInvalid()             # Set s unvalid (will be used for testing)
            if isinstance(source,float) or isinstance(source,int):
                u = Unit3d(Angle(source))
            else:
                u = Unit3d(source)

        jmin = 0                  # Set default to central ray only
        jmax = 1
        imin = 0
        imax = 1
        #                         Sort out range of ray positions in aperture
        if key == "vl":           # Vertical 
            jmin = -nrays
            jmax = nrays + 1
        elif key == "hl":         # Horizontal
            imin = -nrays
            imax = nrays + 1
        elif key == "array":      # array
            jmin = -nrays
            jmax = nrays + 1
            imin = jmin
            imax = jmax
        else:
            print("ray.RayPencil.addBeam: illegal key {0:s}".format(str(key)))

        # Scan through making the rays
        for j in range(jmin,jmax):
            for i in range(imin,imax):
                y = j*dr                                       # x/y in aperture plane  in local coordinates      
                x = i*dr
                if x*x + y*y <= radius*radius:                 # Ignore if outside radius of aperture
                    p = Vector3d(pt.x + x, pt.y + y, pt.z)     # Point in aperture in global coordinates

                    if s:                                      # From source
                         u = Unit3d(p - s)                              
                         ray = IntensityRay(s,u,wavelength,intensity,index)    # Make source ray
                    else:                                      # Collimated beam
                        dist = radius + x*u.x + y*u.y
                        p -= dist*u                                # Propagate point to make it look nicer
                        ray = IntensityRay(p,u,wavelength,intensity,index)     # Make collimated 
                    if path:
                        ray.pathlength = 0.0
                    self.append(ray)                           # Append to self

        return self
        
    
    def addCollimatedBeam(self,ca,u,key = "vl" ,nrays = 10 ,wave = getDefaultWavelength(), intensity = 1.0, index = AirIndex(), path = False):
        """
        Method to add a collimated beam of IntensityRays that fills a specified input apeture.

        :param ca: circular aperture to filled (any object with maxRadius attribute)
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
        :param path: record pathlength, default = False

        """
        if not hasattr(ca, "maxRadius"):
            ca = ca.entranceAperture()
        pt = ca.getPoint()         # Reference point
        radius = ca.maxRadius
        dr = radius/(nrays + 0.1)

        #            Sort out angle (needed internally)
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
        elif key == "hl":         # Horizontal
            imin = -nrays
            imax = nrays + 1
        elif key == "array":      # array
            jmin = -nrays
            jmax = nrays + 1
            imin = jmin
            imax = jmax
        else:
            print("ray.RayPencil.addCollimatedBeam: illegal key {0:s}".format(str(key)))


        # Scan through making the rays
        for j in range(jmin,jmax):
            for i in range(imin,imax):
                y = j*dr                                       # x/y in aperture plane  in local coordinates      
                x = i*dr
                if x*x + y*y <= radius*radius:                 # Ignore if outside radius of aperture
                    p = Vector3d(pt.x + x, pt.y + y, pt.z)     # Point in aperture in global coordinates
                    dist = radius + x*u.x + y*u.y
                    p -= dist*u                                # Propagate point to make it look nicer
                    ray = IntensityRay(p,u,wave,intensity,index)     # Make the ray
                    if path:
                        ray.pathlength = 0.0
                    self.append(ray)                           # Append to self
        
        return self


    def addCollimatedParaxialBeam(self,ca,u,nrays = 10 ,wave = getDefaultWavelength() , intensity = 1.0):
        """
        Method to add a collimated paraxial beam
        param ca  aperture to fill
        param u direction of rays
        param nrays, number or rays aross radius, (default = 10)
        param wave, the wavelength, (default = Default)
        param intensity, the ray intensity, (default = 1.0)
        param path, record path length
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
    def addSourceBeam(self, ca, source, key = "vl" ,nrays = 10 ,wave = getDefaultWavelength(), index = AirIndex(),path = False):
        """
        Method to add beam from a source point that fills an aperture.
        
        :param ca: circular aperture to fill.
        :type ca: CircularAperture o
        :param source: The source point, if Vector3d, intensity will default to 1.0.
        :type source: SourcePoint or Vector3d
        :param key: method of fill, allowed keys as "vl", "hl" and "array",(default is "vl")
        :type key: str
        :param nrays: number or rays aross radius, (Default = 10)
        :type nrays: int
        :param wave: the wavelength, (default = Default)
        :type wave: float
        :param index: Refrative index (Default = AirIndex())
        :type index: RefratciveIndex

        """
        if not hasattr(ca, "maxRadius"):    
            ca = ca.entranceAperture()
        pt = ca.getPoint()                # Reference point
        radius = ca.maxRadius
        dr = radius/(nrays + 0.1)
        s = Vector3d(source)              # Local copy of source location

        #            Sort out intensity of source
        if isinstance(source,SourcePoint):
            intensity = source.getIntensity(wave)
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
                y = j*dr                                          # x/y in local coordinates
                x = i*dr
                if x*x + y*y <= radius*radius:
                    p = Vector3d(pt.x + x, pt.y + y, pt.z)         # Point in aperture in global coordinates
                    u = Unit3d(p - s)                              # Direction of ray
                    ray = IntensityRay(s,u,wave,intensity,index)    # Make ray
                    self.append(ray)                                # Add to pencil
                    if path:
                        ray.pathlength = 0.0
        
        return self


    def addSourceParaxialBeam(self,pg, height, sourceplane, nrays = 10 ,wave = getDefaultWavelength(), intensity = 1.0):
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

    def rotateAboutX(self,angle,origin = None):
        """
        Rotate all the valid rays in the pencil about the X axis by specfied
        angle with respect to given origin (None = (0,0,0))
        
        :param angle: angle in radians
        :type angle: float
        :param origin: rotation origin (Defaault = None)
        :type origin: Vector3d or None
        """
        for r in self:
            if r:
                r.rotateAboutX(angle,origin)
                
        return self

    def propagate(self,distance):
        """
        Method to propagate all rays a equals distance, normaally called via the += operator
        
        :param distance: the distance
        :type distance:  float

        """
        for r in self:
            r.propagate(distance)
        return self
    

    def propagateThrough(self,sur):
        """
        Propagate the whole pencil through the an optical surface, or OpticalGroup. Normally called via \*= operator.

        :param sur: the Surface of OpticalGroup
        :type sur: Surface or OpticlGroup

        """
        for r in self:             # For each ray
            if r:                  # in  each ray is valid
                r.propagateThrough(sur)

        return self

    
    def __imul__(self,surface):
        """
        Implement ___rmul__ to multiply by a suraface
        """
        return self.propagateThrough(surface)
    
    def __iadd__(self,d):
        return self.propagate(d)

    
    def addMonitor(self,monitor = None):
        """
        Method to add/remove a copy of the monitor to each ray.

        :param monitor: The RayMonitor, if None, then the monitor is removed.

        """
        for r in self:
            if monitor == None:
                r.addMonitor()
            else:
                r.addMonitor(monitor.copy())
        return self


    def draw(self):
        """
        Draw each ray in turn to the current plot axis assuming that a RayPath monitor have been added, else it will do nothing.

        """
        for r in self:
            r.draw()
        
               


        
#

class GaussianBeam(Ray):
    """
    Class to work with Gaussian Beams, class uses the underlying Ray class
    """
    def __init__(waist, wave = getDefaultWavelength(), intensity = 1.0):
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

    
        


        




