"""
Set of classes to hold rays for optical ray tracing and paraxial
matrix methods.

Author: Will Hossack, The Univesrity of Edinburgh
"""
import math
from vector import Vector3d,Vector2d,Unit3d
from wavelength import Default,Spectrum,AirIndex,WavelengthColour
from matplotlib.pyplot import plot
#                
#

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

#                Base Ray class whih just hold wavelength and intensity
#                all other (useful) classes extend this Base class
class Ray(object):
    """
    Base Ray class whih just hold wavelength, intensity and refractive index all other (useful) 
    classes extend this Base class
    """
    #            Constuctor with two optional arguments
    #            wavelength (defaults to Defalt (0.55um))
    #            intensity with float of Spectrum, defaults to 1.0
    #            index the local refractive index, defaults to AirIndex
    def __init__(self,wavelength = Default, intensity = 1.0, index = None):
        """
        Constuctor with two optional arguments
        param wavelength float in microns (defaults to Default (0.55um))
        param intensity with float OR Spectrum, (defaults to 1.0)
        param index the refractive index (defaults to AirIndex())
        """
        self.wavelength = float(wavelength)
        if isinstance(intensity,Spectrum):
            self.intensity = intensity.getValue(self.wavelength)
        else:
            self.intensity = float(intensity)
        if index == None:
            self.refractiveindex = AirIndex()      
        else:
            self.refractiveindex = index
        self.monitor = None                 # Hook for ray monitor

    #           toString method to print out basic information
    def __str__(self):
        """
        Implement srt() to give basic imformation, typically overloaded by extending class.
        """
        return "l: {0:8.5e} i: {1:8.5e}".format(self.wavelength,self.intensity)
        
    #
    def __repr__(self):
        """
        Implement repr() to give more detailed information, typiaclly overloaded by extending class.
        """
        return "ray.Ray({0:8.5e}, {1:8.5e})".format(self.wavelength,self.intensity)

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
        param mon the ray monitor, for example RayPath() with will alow for plotting.
        """
        self.monitor = mon
        if self.monitor != None:
            self.monitor.update(self)
        return self
    #
    #
    def isValid(self):
        """
        Method to test if Ray is valid, needs to be defined in extending classes.
        """
        print("Ray.isValid needs to be defined")
    #
    #
    def __nonzero__(self):
        """
        Implement logical __nonzero__ test of Ray is valid
        """
        return self.isValid()
    #
    def __iadd__(self,d):
        """
        Implement _iadd_ to propagate a distance d
        """
        if isinstance(d,float) or isinstance(d,int):
            self.propagate(d)
            return self
        else:
            raise TypeError("rays.Ray __iadd_ type error")
    #
    #

    def __add__(self,d):
        """
        Impement __add__ to  implement a = self + d
        """
        if isinstance(d,float) or isinstance(d,int):
            r = self.copy()
            r.propagate(d)
            return r
        else:
            raise TypeError("ray.Ray __add__ type error")
    #
    #
    def __radd__(self,d):
        """
        Impement __add__ to  implement a = self + d
        """
        if isinstance(d,float) or isinstance(d,int):
            r = self.copy()
            r.propagate(d)
            return r
        else:
            raise TypeError("ray.Ray __radd__ type error")
        
    #
    #
    def __imul__(self,surface):
        """
        Implement ___rmul__ to multiply by a suraface
        """
        self.propagateThrough(surface)
        return self
        

    def __mul__(self,surface):
        """
        Implement ___mul__ to multiply by a suraface
        """        
        r = self.copy()
        r.propagateThrough(surface)
        return r
        

    def __rmul__(self,surface):
        """
        Implement ___mul__ to multiply by a suraface
        """
        r = self.copy()
        r.propagateThrough(surface)
        return r

    def draw(self):
        """
        Implemeent a draw if there is a RayPath minotor in place
        """
        if self.monitor != None and isinstance(self.monitor,RayPath):
            self.monitor.draw()
    
        
#          
class IntensityRay(Ray):
    """
    Class to form a Intensity Ray full vector ray tracing. 
    """
    #
    #
    def __init__(self, pos, dirn = Unit3d(0,0,1), wavelength = Default, intensity = 1.0, index = None):
        """
        Consctructor for to set parameters
        param pos Vector3d, the starting position of the ray, or ParaxialRay
        param dirn Unit3d or Angle, the starting direction of the ray (defaults to (0,0,1))
        param wavelength float (defaults to Default)
        param intensity float or Spectrum (defaults = 1.0)
        param index RefractiveIndex, (defaults to AirIndex())
        """
        if isinstance(pos,ParaxialRay):
            Ray.__init__(self,pos.wavelength,pos.intensity,pos.refractiveindex)
            self.position = Vector3d(0.0,pos.h,pos.z)
            self.director = Unit3d(Angle(pos.u))
        else:
            Ray.__init__(self,wavelength,intensity,index)  # Set wavelnegth intensity and index in super
            self.position = Vector3d(pos)                  # Make localcopy of Position and Dirctor
            self.director = Unit3d(dirn)
        self.pathlength = None                             # Set opl to zero 
    #
    #    
    def __repr__(self):
        """
        Implement repr() to give detailed report an all variables for checking.
        """
        return "IntensityRay: l: {0:8.5f} i: {1:8.5f}\n{2:s}\n{3:s}\nopl: {4:s} n: {5:s}"\
            .format(self.wavelength,self.intensity,repr(self.position),\
                   repr(self.director),str(self.pathlength),repr(self.refractiveindex))
    #
    #
    def copy(self):
        """
        Return a (deep) copy of the current IntesnityRay.
        """
        r = IntensityRay(self.position,self.director,\
                         self.wavelength,self.intensity,\
                         self.refractiveindex.copy())
        r.pathlength = self.pathlength
        return r

    #
    #
    def setInvalid(self):
        """
        Method to set the ray to inValid (sets the Unit3d as invalid)
        """
        self.director.setInvalid()       # Set director to be invalid
        return self
    #
    #
    def isValid(self):
        """
        Method to check the ray is Valid (checks the Director is valid)
        """
        return self.director.isValid()
    #
    #
    def getPhaselength(self):
        """
        Method to get the phase length, being 2*pi*pathelength/wavelength
        """
        return 2000.0*math.pi*self.pathlength/self.wavelength
    #
    #       
    def propagate(self,distance):
        """
        Method to propagate the ray a specifed distance using its own current direction.
        This also upadtes the pathlength and is the main method to propgate rays.
        param distance float the distance to propagate the ray.

        Method clecks if the ray is valid, and if sp propagates it.
        """

        if self :
            self.position += distance*self.director

            if self.pathlength != None:
                self.pathlength += distance*self.refractiveindex.getValue(self.wavelength)

            if self.monitor != None:              # Update monitor of it exits
                self.monitor.update(self)
       

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
        param surface, Surface or list of list(Surface), if list each one is dealt with in order.
        return boolean, true is passed through, false if blocked. If blocked the ray will
        be invalid and it position will be set to the intersetcion point here it was blocked.
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
            return self.director.reflection(info.normal)    # Do reflection

        else:            
            raise TypeError("IntensityRay from unknow surface type {0:d}".format(info.type))

        return False                                   # If here we have failed (somehow), return false

    #
    #
    def pointInPlane(self,plane):
        """
        Method to calcualte where this ray striked an optical place
        This does NOT alter othe current ray.
        param plane, the OpticalPlane or z the location on the optical axis
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
class RayPath(object):
    """
    Class to record a ray path. the path in held in three lists x[], y[] and z[]
    """
    def __init__(self,ray = None):
        self.x = []
        self.y = []
        self.z = []
        self.wavelength = Default                    # The default
        if ray != None:
            ray.addMonitor(self)
            self.wavelength = ray.wavelength
        
    def copy(self):
        """
        Method to form a copy of itself
        """
        return RayPath()

    def __str__(self):
        return "Ray Monitor"

    def update(self,ray):
        """
        udate method,  records the x,y,z position of the  ray is lists
        """
        self.wavelength = ray.wavelength
        if isinstance(ray,ParaxialRay):
            self.y.append(ray.h)
            self.z.append(ray.z)
        else:
            self.x.append(ray.position.x)
            self.y.append(ray.position.y)
            self.z.append(ray.position.z)
        

    def draw(self):
        """
        Do a plot to pyplot with colour of ray given by its wavelength defiined in wavelength.WavelengthColour
        """
        col = WavelengthColour(self.wavelength).hexString()
        return plot(self.z,self.y,col)
    

#
class RayPencil(list):
    """
    Class to hold a list of rays and implement methods to propagate on-mass; this will work for any
    type of Ray; at present IntensityRay and ParaxialRays (of a  mixture of the two)
    """
    
    def __init__(self, *args):
        """
        Make a ray pencil with optional set of rays to be appended.
        """
        list.__init__(self)

        for r in args:
            self.append(r)
    #
    def addCollimatedBeam(self,ca,u,key = "vl" ,nrays = 10 ,wave = Default, intensity = 1.0):
        """
        Method to add a collimated beam
        param ca circular aperture to fill
        param u direction of rays
        param key method of fill, allowed keys as "vl", "hl" and "array",(default is "vl")
        param nrays, number or rays aross radius, (default = 10)
        param wave, the wavelength, (default = Default)
        param intensity, the ray intensity, (default = 1.0)
        """

        if not hasattr(ca, "maxRadius"):
            ca = ca.entranceAperture()
        pt = ca.getPoint()         # Reference point
        radius = ca.maxRadius
        dr = radius/(nrays + 0.1)
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
        param ca circular aperture to fill
        param u direction of rays
        param nrays, number or rays aross radius, (default = 10)
        param wave, the wavelength, (default = Default)
        param intensity, the ray intensity, (default = 1.0)
        """

        if not hasattr(ca, "maxRadius"):
            ca = ca.entranceAperture()
        pt = ca.getPoint().z
        radius = ca.maxRadius
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
        Draw each ray
        """
        for r in self:
            r.draw()

#               
class ParaxialRay(Ray):
    """
    Class for Paraxial Rays
    """
    #    
    def __init__(self,height = 0.0, angle = 0.0, plane = 0.0, \
                 wavelength = Default ,intensity = 1.0, index = None):
        """
        Constuctor with 5 optional arguments
        param height (defaults to 0.0) height from optical axis
        param angle (defaults to 0.0) angle in radians from optical axis
        param plane (defaults to 0.0) location of plane along optical axis
        param wavelength (defaults to Default) wavelength in microns
        paramintensity (defaults to 1.0) intensity of ray
        """
        Ray.__init__(self,wavelength,intensity, index)   # Set wavelength and intensity
        self.h = height               
        self.u = angle
        self.z = plane
    #
    # 
    def __repr__(self):
        """
        Implement repr() details with full name
        """
        return "ray.ParaxialRay{0:s}".format(str(self))
    #
    def __str__(self):
        """
        Implment str() print parameters
        """
        return "({0:7.5f}, {1:7.5f}, {2:7.5f}, {3:7.5f}, {4:7.5f})".\
            format(self.h,self.u,self.z,self.wavelength,self.intensity)
    #            
    def copy(self):
        """
        Method to make a deep copy of a current ParaxialRay
        """        
        return ParaxialRay(self.u,self.u,self.z,\
                           self.wavelength,self.intensity)
    #           
    def setInvalid(self) :
        """
        Method to set an ParaxialRay to inValid
        """
        self.u = float("nan")       # set angle to be NaN 
        return self

    #           
    def isValid(self):
        """
        Method to test of a Paraxial Ray is valid
        """
        return not math.isnan(self.u)

    #            
    def propagate(self, distance):
        """
        Method to propagate a ray a specified distance
        param distance, the distance the ray is propagated
        returns True is sucessful, False is Ray is invalid
        """
        if self.isValid() :
            self.z += distance                # update plane
            self.h += self.u*distance         # calculate new height
            if self.monitor != None:          # Update monit of it exits
                self.monitor.update(self)
            return True
        else:
            return False

    
    def propagateToPlane(self, plane):
        """
        Method to propagate the rays to a specified plane
        param plane, the location of plane
        return True is sucessful, False if Ray is invalid
        """
        if  not math.isinf(plane):
            distance = plane - self.plane         # distance to propagate
            return self.propagate(distance)       # do the propagatation
        else:
            return False                          # plane at inf, so fail

    
    def mult(self,m):
        """
        Method to multiply ParaxialRay by ParaxialMatrix and
        return new ParaxialRay
        """
        if self.isValid():
            h = self.h*m.A + self.u*m.B
            a = self.h*m.C + self.u*m.D
            p = self.z + m.thickness
            return ParaxialRay(h,a,p,self.wavelength,self.intensity)
        else:
            return self            
    
    
    def multBy(self,m):
        """
        Method to multiply ParaxialRay by ParaxialMatrix in place.
        """
        if self.isValid():
            h = self.h*m.A + self.u*m.B
            a = self.u*m.C + self.u*m.D
            self.h = h
            self.u = a
            self.z += m.thickness
            if self.monitor != None:              # Update monit of it exits
                self.monitor.update(self)
            return True
        else:
            return False

    #            
    #            
    def propagateThrough(self, surface):
        """
        Method to propagate Paraxial Ray through a surface or list of surfaces.
        this the main paraxial ray tracing method with same call as for the skew IntensityRays
        param surface, or list of surfaces
        """

        if isinstance(surface,list):
            for s in surface:                # process each surface in the list in turn
                b = self.propagateThrough(s)
                if not b:                    # failed of surface, don't do anymore
                    break
                
            return b

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
        other is other ParaxialRay
        """
        dtheta = other.u - self.u
        if dtheta == 0.0:
            return float('inf')        #   Trap infinity where rays are parallel
        else:
            return (self.h - other.h - self.z*self.u + other.z*other.u)/dtheta
#
#            
#
class ParaxialMatrix(object):
    """
     ParaxialMartrix for propagation between input and output planes
    """
    #        
    #
    def __init__(self, a_or_m = 1.0, b = 0.0, c = 0.0, d = 1.0, t = 0.0):
        """
        Constructor with five parameters, the four matrix elements and thickness
        a (default to 1.0)
        b (defaults to 0.0)
        c (defaults to 0.0)
        d (defaults to 1.0)
        t (defaults to 0.0)
        """
        if isinstance(a_or_m,ParaxialMatrix):
            m = a_or_m
            self.A = m.A
            self.B = m.B
            self.C = m.C
            self.D = m.D
            self.thickness = m.thickness
        else:
            a = a_or_m
            self.A = float(a)
            self.B = float(b)
            self.C = float(c)
            self.D = float(d)
            self.thickness = t


    def __repr__(self):
        """
        Convert to a string
        """
        return "ParaxialMatrix : [{0:7.5f},{1:7.5f}]\n                 [{2:7.5f},{3:7.5f}] t: {4:7.5f}"\
            .format(self.A,self.B,self.C,self.D,self.thickness)

    #      
    def copy(self):
        """
        Method to make a copy of the ParaxialMatrix
        """
        return ParaxialMatrix(self)

    #          
    def trace(self):
        """
        Return the trace of the matrix
        """
        return self.A + self.D

    #          
    def determinant(self):
        """
        Return the determinant of the matrix
        """
        return self.A*self.D - self.B*self.C



    def inverse(self):
            """
            Form the Inverse of the matrix, it also set the thickness to be -ve
            """
            det = self.determinant()
            a = self.D/det
            self.B /= -det
            self.C /= -det
            self.D = self.A/det
            self.A = a
            self.thickness = -self.thickness

            return self
    #          
    def scale(self,a):
        """
        Method to scale the matrix by a factor a
        element B and thickness as scaled by a, C is dvided by a, A/D not changed
        """
        self.B *= a
        self.C /= a
        self.thickness *= a
        return self

    #          
    def backPower(self):
        """
        Method to get the back power
        """
        return float(-self.C)

    #          
    def backFocalLength(self):
        """
        Method to get the back focal length
        """
        return 1.0/self.backPower()

    #         
    def backFocalPlane(self):
        """
         Method to get the back focal plane (relative to output plane)
        """
        return -self.A/self.C

    #          
    def backPrincipalPlane(self):
        """
        Method to get the back principal plane (relative to the output plane)
        """
        return (1.0 - self.A)/self.C

   #         
    def frontPower(self):
        """
         Method to get the front power
        """
        return float(self.C/self.determinant())

    #          
    def frontFocalLength(self):
        """
        Method to get the front focal length
        """
        return 1.0/self.frontPower()

    #          
    def frontFocalPlane(self):
        """
        Method to get the front focal plane (relative to input plane)
        """
        return self.D/self.C

    #          
    def frontPrincipalPlane(self):
        """
        Method to get the front principal plane (refaltive to the input plane)
        """
        return (self.D - self.determinant())/self.C

    #          
    def setFrontFocalLength(self, f):
        """
        Method to set the front focal length by scaling 
        """
        self *= f/self.frontFocalLength()

    #          
    def setBackFocalLength(self, f):
        """
        Method to set the back focal length by scaling 
        """
        self *= f/self.backFocalLength() 


    #         
    def __mul__(self,m):
        """
        Method to pre-multiply the current matrix by a another Paraxialmatrix, or if a float
        it is scaled 
        return new ParaxialMatrix
        """
        if isinstance(m,ParaxialMatrix):
            a = m.A*self.A + m.B*self.C
            b = m.A*self.B + m.B*self.D
            c = m.C*self.A + m.D*self.C
            d = m.C*self.B + m.D*self.D
            t = m.thickness + self.thickness
        elif isinstance(m,float) or isinstance(m,int):    # Scale
            a = self.A
            b = self.B * m
            c = self.C / float(m)
            d = self.D
            t = self.thickness * m
        else:
            raise TypeError("ParaxialMatrix *= call with unknown type " + str(m))
        
        return ParaxialMatrix(a,b,c,d,t)

    def mult(self,m):
        """
        Method to pre-multiply the current matrix by a another Paraxialmatrix
        return new ParaxialMatrix
        """
        a = m.A*self.A + m.B*self.C
        b = m.A*self.B + m.B*self.D
        c = m.C*self.A + m.D*self.C
        d = m.C*self.B + m.D*self.D
        t = m.thickness + self.thickness
        return ParaxialMatrix(a,b,c,d,t)

    #          
    def __imul__(self,m):
        """
        Method to pre-multiply the currenr matrix by a another Paraxialmatrix in place
        """
        if isinstance(m,ParaxialMatrix):
            a = m.A*self.A + m.B*self.C
            b = m.A*self.B + m.B*self.D
            c = m.C*self.A + m.D*self.C
            d = m.C*self.B + m.D*self.D
            self.A = a
            self.B = b
            self.C = c
            self.D = d
            self.thickness += m.thickness
        elif isinstance(m,float) or isinstance(m,int):    # Scale
            self.scale(m)
        else:
            raise TypeError("ParaxialMatrix *= call with unknown type " + str(m))
        return self

    def multBy(self,m):
        """
        Method to pre-multiply the current matrix by a another Paraxialmatrix in place
        """
        a = m.A*self.A + m.B*self.C
        b = m.A*self.B + m.B*self.D
        c = m.C*self.A + m.D*self.C
        d = m.C*self.B + m.D*self.D
        self.A = a
        self.B = b
        self.C = c
        self.D = d
        self.thickness += m.thickness
        return self

    #
    #
    def __add__(self,d):
        """
        Short hand to implement a propagation matrix 
        param d float the distance
        """
        m = PropagationMatrix(d)
        return self*m

    def __iadd__(self, d):
        """
        Short hand to implement a propagation martix in place
        param d the distance
        """
        if isinstance(d,float) or isinstance(d,int):
            m = PropagationMatrix(d)
            self *= m
            return self
        else:
            raise TypeError("ParaxialMatrix += call with illegal type : " + str(d))

#
#            
class ThinLensMatrix(ParaxialMatrix):
    """
    ParaxialMatrix for a thin lens 
    """
    #       
    def __init__(self,f_or_cl,n=None,cr=None):
        """
        Constuctor with one or three paramters
        param f_or_cl float flocal length or left curvature
        param n refractive index (defaults None)
        param cr reight curvaure
        """
        if n == None:                   # Only one parameter
            ParaxialMatrix.__init__(self,1.0 , 0.0 , -1/f_or_cl , 1.0, 0.0)
        else:
            a = DielectricMatrix(1.0,n,f_or_cl) 
            b = DielectricMatrix(n,1.0,cr) 
            ParaxialMatrix.__init__(self,a*b)

#
#        
class ThickLensMatrix(ParaxialMatrix):
    """
     ParaxialMatrix for a thick lens.
    """
    #
    def __init__(self, cl , n , t , cr):
        """
        Constructor with 4 parameters, all reqired
        param cl float left curvature
        param n  float refractive index
        param t  float thickness of lens
        param cr float right curvature
        """
        a = DielectricMatrix(1.0,n,cl)     # Front surface
        b = PropagationMatrix(t)           # Thickness
        c = DielectricMatrix(n,1.0,cr)     # Back surface
        s = a * b * c                      # Total matrix
        ParaxialMatrix.__init__(self,s)    # Set self

#
#           
class PropagationMatrix(ParaxialMatrix):
    """
     ParaxialMatrix for propagation
    """
    #      
    def __init__(self,d):
        """
        param d float the propagation distance
        """
        ParaxialMatrix.__init__(self,1.0 , d , 0.0 , 1.0, d)

#
#             
class DielectricMatrix(ParaxialMatrix):
    """
    DielectricMatrix for flat or curved interface
    """
    #       
    def __init__(self, nl, nr, c = 0.0):
        """
        Constructor with three parameters
        param nl float refractive index on left of interface
        param nr float refractive index on right of interface
        param c float  curvature of interface. (default = 0.0)
        """
        ParaxialMatrix.__init__(self, 1.0 , 0.0, c*(nl - nr)/nr , nl/nr, 0.0)


#            
class ParaxialGroup(ParaxialMatrix):
    """
    ParaxialGroup class to represent paraxial group with input plane, ParaxialMatrix and input and outplane heights.
    """
    #
    #
    def __init__(self, p = 0.0, matrix = ParaxialMatrix(), inh = float("inf"), outh = None):
        """
        Create the ParaxialGroup, 
        param p float the input plane of the matrix (default = 0.0)
        param matrix the paraxial matrix, (defaults to unit matrix)
        param inh float height of input plane, (defaults to float("inf"))
        param inh float height of outut plane, (deafults to None, which is copy of input)
        """
        ParaxialMatrix.__init__(self,matrix)
        self.inputPlane = float(p)             # Input plane
        self.inputPlaneHeight = inh
        if outh == None:
            self.outputPlaneHeight = self.inputPlaneHeight
        else:
            self.outputPlaneHeight = outh


    #    
    def __repr__(self):
        """
        Implment repr()
        """
        return "ray.ParaxialGroup : i : {0:7.5f} hi : {1:7.5f} ho : {2:7.5f}\n{3:s}".\
            format(self.inputPlane, self.inputPlaneHeight, self.outputPlaneHeight, \
                   ParaxialMatrix.__repr__(self))

    #          Method to make a deep copy of the current Paraxial Group
    def copy(self):
        return ParaxialGroup(self.inputPlane,self,self.inputPlaneHeight,self.outputPlaneHeight)

    #          
    def getInputPlane(self):
        """
        Method to get input plane
        """
        return self.inputPlane

    #          
    def getOutputPlane(self):
        """
        Method to get output plane
        """
        return self.inputPlane + self.thickness    # always calculate

    #          
    def scale(self,a):
        """
        Method to scale the matrix and plane heights 
        """
        ParaxialMatrix.scale(self,a)
        self.inputPlaneHeight *= a
        self.outputPlaneHeight *= a
    #          
    #          
    def backFocalPlane(self):
        """
        Method to get the back focal plane in global coodinates
        """
        return self.getOutputPlane() + ParaxialMatrix.backFocalPlane(self)

    def backNodalPoint(self):
        """
        Get the back Nodal point, (normally same as back principal plane)
        """
        return self.backFocalPlane() + self.frontFocalLength()

    #       
    def backPrincipalPlane(self):
        """
        Method to get the back principal plane 
        """
        return self.getOutputPlane() + ParaxialMatrix.backPrincipalPlane(self)

    #          
    def frontFocalPlane(self):
        """
        Method to get the front focal plane in global coodinates
        """
        return self.getInputPlane() + ParaxialMatrix.frontFocalPlane(self)

    #          
    def frontPrincipalPlane(self):
        """
        Method to get the front principal plane in global coordinates
        """
        return self.getInputPlane() + ParaxialMatrix.frontPrincipalPlane(self)

    def frontNodalPoint(self):
        """
        Get Front Modal point, (normall the same as front Prinicpal Plane)
        """
        return self.frontFocalPlane() + self.backFocalLength()
    #          
    def imagePlane(self, op):
        """
        Method to get the image plane for specified object plane using geometric lens fomula
        and properties of the current group.
        param op float location on object plane on optical axis
        return float position on image plane.
        """
        u = self.frontPrincipalPlane() - float(op)      # distance from front principal plane
	v = u/(self.backPower()*u - 1.0)           # distance from back principal plane
	return  self.backPrincipalPlane() + v          # where the image is
    
    #
    def planePair(self,mag):
        """
        Calcualte the object image plane pair for specifed magnification in global coordinates.
        param mag float the magnification (note most imaging system mag is -ve)
        return list of [op,ip] being the location of the obeect and image places respectively
        """
        f = self.backFocalLength()
        u = f*(1.0 - 1.0/mag)
        v = f*(1.0 - mag)
        op = self.frontPrincipalPlane() - u
        ip = self.backPrincipalPlane() + v
        return op,ip


    def imagePoint(self,op):
        """
        Method to three-dimensional image of a point in obejct space in global coordinates
        geometric optics.
        param op Position, object point, can also ve Unit3d or Angle where it will assume an onfinite object
        return Position the image point.
        """
        if isinstance(op,Unit3d):                         # Infinte object
            p = Vector3d(0,0,self.backNodalPoint())
            return Vector3d(p + op*(self.backFocalLength()/op.z))

        elif isinstance(op,Angle):                        # Also infinite object
            return self.pointImage(Unit3d(op))

        elif isinstance(op,Vector3d):                       # Finite object
            u = self.frontPrincipalPlane() - op.z
            v = u/(self.backPower()*u - 1)
            ip = self.backPrincipalPlane() + v               # Image location
            mag = -v/u
            return Vector3d(mag*op.x , mag*op.y , ip)
       
        else:
            raise TypeError("ray.ParaxialGroup.pointImage: called with unknown type {0:s}".format(str(op))) 

    #
    def cardinalPoints(self):
        """ 
        Method to get the 6 cardinal points as a list, order is:
        0:     Front Focal Point
        1:     Back Focal Point
        2:     Front principal plane
        3:     Back principal plane
        4:     Front nodal point (front principal plane for system in air)
        5:     Back nodal point (back principal plane for system in air)
        """
        return [self.frontFocalPlane(),self.backFocalPlane(),\
                self.frontPrincipalPlane(),self.backPrincipalPlane(),\
                self.frontNodalPoint(),self.backNodalPoint()]


   
    def draw(self):
        """
        Draw 4 cardinal planes to plot
        """
        y = [-self.inputPlaneHeight,self.inputPlaneHeight]        # The plane heights
        ff = self.frontFocalPlane()       
        zff = [ff,ff]
        fp = self.frontPrincipalPlane()
        zfp = [fp,fp]
        bp = self.backPrincipalPlane()
        zbp = [bp,bp]
        bf = self.backFocalPlane()
        zbf =  [ bf,bf]
        plot(zff,y,"b",zfp,y,"r",zbp,y,"g",zbf,y,"b")



#            Class to make a ParaxialAperture
class ParaxialAperture(ParaxialGroup):

    #         Constuctor with only position and height
    def __init__(self,p,h):
        m = ParaxialMatrix()          # Default identity matrix
        ParaxialGroup.__init__(self,m,p,h)  # set matrix, position and input height
        

    

#            ParaxialSystem with extends list
class ParaxialSystem(list):
    
    #
    #        Constructor with optional single ParaxialGroup
    def __init__(self,pg=None):
        list.__init__(self)
        if isinstance(pg,list) :        # Supplied a list
            for g in pg:
                self.append(g)        # Appled each

        if pg != None:                # Append single suppiled Group
            self.append(pg)

    #     Method to get input plane (input plane of first group)
    def getInputPlane(self):
        return self[0].getInputPlane()
        
    #    Method to get the output plane (output plane of last element)
    def getOutputPlane(self):
        return self[-1].getOutputPlane()

    #    Method to get the overall ParaxialGroup of the system.
    def getParaxialGroup(self):
        gr = self[0].copy()         # copy of first element
        for g in self[1:]:          # 
            d = g.getInputPlane() - gr.getOutputPlane() # Distance to input
            gr.matrix.multBy(PropagationMatrix(d))      # Propagate to input
            gr.matrix.multBy(g.matrix)                  # do mult of matrix
            gr.outputPlaneHeight = g.outputPlaneHeight  # set ouput height

        return gr                   # Return the group

        
#

class GaussianBeam(ParaxialRay):
    """
    Class to work with Gaussian Beams
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

    
        


        




