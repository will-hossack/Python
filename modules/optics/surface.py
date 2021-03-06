"""
Set of classes to implement various types of optical surface.

"""
from vector import Vector2d,Vector3d,Unit3d
from optics.ray import SourcePoint
from optics.matrix import ParaxialGroup
import math
from matplotlib.pyplot import plot

"""
Define the three types of surface.
"""

Clear = 0        #: Defeine a clear surface.
Refracting = 1   #: Define a refracting surface.
Reflecting = 2   #: Define a refelcting surface,

"""
Define surface plotting parameter
"""
SurfacePlotPoints = 10
#
"""
Define global definition of a blocked ray.
"""
Blocked = Unit3d()


class SurfaceInteraction(object):
    """
    Class to hold the interaction of a vector skew ray with a general surface.
    It contains all the information needed to calualte the interaction and update the ray.

    :param type: the type of surface, Clear = 0, Refracting = 1, Reflecting = 2.
    :type type: int
    :param point: the surface reference point in global coordinate.
    :type point: vector.Vector3d
    :param distance: distance from current ray position to the surface (in the direction of the ray)
    :type distance: float
    :param position: interaction point on surface. (where ray strikes surface)
    :type position: vector.Vector3d
    :param normal: surface normal at surface interaction point 
    (this should be set to invalid if ray is blocked by the surface).
    :type normal: vector.Unit3d
    :param index: refractive index on the image side of surface, (this should be None unless the surface type 
    is refracting.)
    :type index: optics.wavelength.RefractiveIndex

    """
    
    def __init__(self,type,point,distance,position,normal,index):
        """
       
        """
        self.type = type           
        self.point = point
        self.distance = distance
        self.position = position
        self.normal = normal
        self.refractiveindex = index

    #
    def __repr__(self):
        """
        Detailed reprecentation for debug purposes.
        """
        s = "surface.SurfaceInteraction: Type : {0:d}\nSurfacePoint: {1:s}\nDistance : {2:8.4e}\nPosition : {3:s}\n"\
            .format(self.type,str(self.point),self.distance,str(self.position))
        s += "Normal : {0:s}\nRefractiveIndex : {1:s}".format(str(self.normal),repr(self.refractiveindex))
        return s
        
 

class Surface(object):
    """
    Base class for an arbitrary surface which needs to be extended to be useful. 
    This class defines.surface reference point, surface type
    refractive index on image size (may be None) and membership of optical group.

    :param pt: the surface reference point, defaults to (0,0,0)
    :type pt: vector.Vector3d or float
    :param type: the surface type (default to Clear 0)
    :type type: int
    :param index: the refractive index on the image side of the surface (default to None)
    :type index: optics.wavelength.RefractiveIndex

    """
    #
    #
    def __init__(self,pt = 0.0 ,type = Clear, index = None):
        """
        Constructor to form a basic surface
        
        """
        
        self.setPoint(pt)
        self.type = type
        #
        self.group = None
        self.refractiveindex = index
        


    def __str__(self):
        """     Basic info as str
        """
        return "spt: {0:s} type: {1:d} n: {2:s}".format(str(self.point),self.type,str(self.refractiveindex))
        
    
    def __repr__(self):
        """
        Implement the repr() method
        """
        return "{0:s} ".format(self.__class__.__name__) + str(self)
        

    def setPoint(self,pt = 0.0):
        """
        Method to set the surface reference point in a consistent way for all surfaces. The specified point can be copy from
        Surface, Vector3d, list, truple, Vector2d or scalar. 
       
        :param pt: surface point can be Surface, Vector3d, 3 component list/truple, Vector2d with z = 0, float giving (0,0,z), (Default = 0.0) 
        :type z_or_v: vector.Vector3d, list[], Vector2d or float
        :return: self
                           
        """
        if isinstance(pt,Surface):
            self.point = pt.point.copy()
        elif isinstance(pt,Vector3d) or isinstance(pt,list) or isinstance(pt,tuple):
            self.point = Vector3d(pt)
        elif isinstance(pt,Vector2d):
            self.point = Vector3d(pt.x,pt.y,0.0)
        elif isinstance(pt,float) or isinstance(pt,int):
            self.point = Vector3d(0.0,0.0,pt)
        else:
            raise TypeError("surface.Surface.setPoint: called with unknown type")
        
        return self


    def movePoint(self,delta):
        """
        Move the reference point by about delta

        :param delta: distance moved
        :type delta: Vector3d or float
        """
        if isinstance(delta,float) or insinstance(delta,int):
            self.point += Vector3d(0.0,0.0,float(delta))
        else:
            self.point += delta
        self.paraxial = None


    def scale(self,a):
        """
        Scale surface, this will scale the surface point only, this is 
        typically extended for more complex sufaces
        
        :param a: scale value
        :type a: float
        :return: self
        """
        self.point *= a
        return self

    def getPoint(self):
        """
        Method to get the surface point in global coordinates taking account
        that the Surface my belong to an optics.lens.OpticalGroup
        
        :return: vector.Vector3d the reference point in global coordinates.

        Note: this should always be used rather than direct reference to surface.point
        """
        if self.group == None:
            return self.point
        elif self.group.group == None:
            return self.group.point + self.point
        else:
            return self.group.group.point + self.group.point + self.point

    #
    #
    def makeStandAlone(self):
        """
        Method to make the Surface standalone, it is removed from any OpticalGroup
        and the reference point is re-set in gobal coordinates.
        
        :return: self
        """
        if self.group != None:
            self.point = self.getPoint()
            self.group.remove(self)
            self.group = None
        return self

    #
    #
    def getNormal(self,pt):
        """
        Abstract class to get normal at specified point which is assume to
        be on the surface.
        
        :param pt: Point where surface normal is calcualted, assumed to be on the surface)
        :type pt: vector.Vector3d
        :return: vector.Unit3d  set to Invalid.
        
        """
        return Unit3d()
    
    def getDistance(self,r,u):
        """
        Astract class to get distance from specified poin
        """
        return float("nan")

    
    def getSurfaceInteraction(self,ray):
        """
        Method to get get  surface interaction information for a ray.
        
        :param ray: the Ray
        :type ray: optics.ray.Ray
        :return: SurfaceInteration with reference point and other parameters invalid
        
        This is abstract surface, so only type, surface point and refarctive index will be valid
        """
        pt = self.getPoint()

        return SurfaceInteraction(self.type,pt,float("nan"),Vector3d().setInvalid(),Blocked,self.refractiveindex)



    def draw(self,option = None):
        """
        Abstract method to draw surface, needs to be defined.
        """
        print("surface.Surface.draw error Need to implement draw")

#
#
class FlatSurface(Surface):
    """
    Class to implement an Flat Surface with specifed surface normal.

    :param pt: reference point, (defaults to (0,0,0))
    :type pt: vector.Vextor3d or float
    :param normal: the surface normal (Default = (0,0,1)
    :type normal: vector.Unit3d
    :param type: Surface type  (defaults to Clear = 0)
    :type type: int
    :param index: the refractive index (defaults to None)
    :type index: optics.wavelength.RefractiveIndex

    """
    #
    #
    def __init__(self,pt = 0.0, normal = Unit3d(0,0,1), type = Clear,index = None): 
        """
        Constructor
       
        """
        Surface.__init__(self,pt,type,index)
        self.normal = Unit3d(normal)
        self.curvature = 0.0


    def __str__(self):
        """
        Implement str()
        """
        return "spt: {0:s} type: {1:d} u: {2:s} n: {3:s}".\
            format(str(self.point),self.type,str(self.normal),str(self.refractiveindex))

        
        
    def getNormal(self,pt):
        """
        Get surface normal at any point.
        """
        return self.normal

    #
    #    
    def getDistance(self,r,u):
        """
        Get the distance from specifed Postition to the surface

        :param r: the position
        :type r: vector.Vector3d
        :param u: the direction
        :type u: vector.Unit3d
        :return: float, the distance

        Not normally called by used, usually called via surfaceInteraction()
        
        """
        p = self.getPoint()      # Reference point
        dv = p - r
        dist = self.normal.dot(dv)/self.normal.dot(u)
        return dist
        

    def getSurfaceInteraction(self,ray):
        """
        Method to get back the surface interaction information for a ray.

        :param ray: the ray
        :type ray: optics.ray.IntensityRay
        :return: SurfaceInteraction

        """
        distance = self.getDistance(ray.position,ray.director)
        pos = ray.position.propagate(distance,ray.director)
        pt = self.getPoint()
        
        return SurfaceInteraction(self.type,pt,distance,pos,self.normal,self.refractiveindex)
       
#
#
class OpticalPlane(FlatSurface):
    """
    Class to implement a flat optical place normal to the optical, this is simpler and
    version of FlatSurface with a fixed normal along the optical (z)-axis.

    :param pt: the surface point, (Defaults to (0.0,0.0,0.0)
    :type pt: vector.Vector3d or float
    :param type: surface type, defaults to Clear (0)
    :type type: int
    :param index: refrative index image side of surface, (Default = None)
    :type index: `optics.wavelength.RefractiveIndex`
    
    """
    
    def __init__(self,pt = 0.0 ,type = Clear, index = None):
        """
        Conststrucor for Optical plane
       
        """
        FlatSurface.__init__(self,pt,Unit3d(0,0,1),type,index)

    def __str__(self):
        """
        Implement str()
        """
        return "spt: {0:s} type: {1:d} n: {2:s}".format(str(self.point),self.type,str(self.refractiveindex))
        
    
    def surfaceVector(self,pos):
        """
        Method to get the Vector2d in the plane for a specified point, assumed to
        be in the plane.

        :param pos: three dimensional point in global coordinates.
        :type pos: :class:`vector.Vector3d`
        :return: :class:`vector.Vector2d` point on plane

        """
        p = pos - self.getPoint()
        return Vector2d(p.x,p.y)


    def getSourcePoint(self,pt_or_x,y = None,intensity = 1.0):
        """
        Get the SourcePoint for a specified point in the plane.
        
        :param pt_or_x: Vector2d point in the plane or x component
        :type pt_or_x: Vector2d or float
        :param y: y component (Default = Null)
        :type y: float
        :param intensity:

        """

        if isinstance(pt_or_x,Vector2d):
            x = pt_or_x.x
            y = pt_or_x.y
        else:
            x = float(pt_or_x)
            y = float(y)
        
        refpt = self.getPoint()
        x += refpt.x
        y += refpt.y
        return SourcePoint([x,y,refpt.z],intensity)
    
    
    
    def getDistance(self,r,u):
        """
        Get the distance from specifed Postition to the surface

        :param r: the position
        :type r: vector.Vector3d
        :param u: the direction
        :type u: vector.Unit3d
        :return: float, the distance

        Not normally called by used, usually called via surfaceInteraction()
        
        """
        pt = self.getPoint()
        d = (pt.z - r.z)/u.z
        return d
    
    def getSurfaceInteraction(self,ray):
        """
        Method to get the SurfaceInteraction information with a Ray
        
        :param ray: the input ray
        :type ray: :py:class:`optitcs.ray.IntensityRay`
        :return: `SurfaceInteraction`

        """
        pt = self.getPoint()
        #
        #       Distance calcaultion simplified.
        distance = (pt.z - ray.position.z)/ray.director.z
        pos = ray.position.propagate(distance,ray.director)
 
        return SurfaceInteraction(self.type,pt,distance,pos,self.normal,self.refractiveindex)
        

    def getParaxialInteraction(self,ray):
        p = self.getPoint()
        distance = p.z - ray.z
        height = ray.h + distance*ray.u

        return [self.type,distance,height,0.0,self.refractiveindex]

    def paraxialGroup(self, option = None):
        """
        Get the ParaxialGroup of the surface for compatibility, it will be reference point
        and unit matrix.
        """
        pt = self.getPoint()
        return ParaxialGroup(pt.z)
        


    def incrementSurface(self,delta):
        """
        Get a new optical surface where the reference point is incremented a distance along
        the optical (z-axis). 

        :param delta: The distance incremented along the optical axis
        :type delta: float
        :return: new OpticalPlane

        """
        return OpticalPlane([self.point.x,self.point.y,self.point.z + delta],self.type,self.refractiveindex)
        

    def draw(self,height=10.0,option = None):
        """
        Define Draw with extra height parameter.
        """
        p = self.getPoint()
        yvals =  [p.y + height, p.y - height]
        zvals =  [p.z,p.z]
        return plot(zvals,yvals,"k",lw=2.0)

    
class CircularAperture(OpticalPlane):
    """
    Class to  implement a circular aperture of fixed radius, rays outside the aperure will be blocked.

    :param pt: the plane reference point, (Defaults =(0.0,0.0,0.0))
    :type pt: Vector3d or float
    :param radius: radius of aperture (Default = 1.0)
    :type radius: float

    """
    def __init__(self,pt = 0.0 , radius = 1.0):
        """
        Constuctor for a circular aperture
       
        """
        OpticalPlane.__init__(self,pt)
        self.outerRadius = radius
        self.maxRadius = radius                          # add max radius to standardise external refs

    def __str__(self):
        """
        Implements str()
        """
        return "spt: {0:s} rad: {1:7.3f}".format(str(self.point),self.outerRadius)
    
    def getRadius(self):
        """
        Get the radius, for thiis it is maxRadius
        """
        return self.maxRadius

    def scale(self,a):
        """
        Scale aperture, will scale position and radius. if a < 0, then abs(a) is used.
        
        :param a: Scale factor
        :type a: float

        """
        a = abs(a)                       # Allow for 
        OpticalPlane.scale(self,a)
        self.outerRadius *= a
        self.maxRadius *=a
        return self

    #
    def getNormal(self,pos):
        """
        Method to get normal at specifed point, this will block rays outside the aperture
        """
        r = self.surfaceVector(pos)
        if r.absSquare() <= self.outerRadius*self.outerRadius:
            return self.normal
        else:
            return Blocked



    def getSurfaceInteraction(self,ray):
        """
        Method to get back the surface interaction information with a Ray
        
        :param ray: the input ray
        :type ray: :py:class:`optitcs.ray.IntensityRay`
        :return: `SurfaceInteraction`
        
        """
        pt = self.getPoint()
        distance = (pt.z - ray.position.z)/ray.director.z
        pos = ray.position.propagate(distance,ray.director)

        #            dx and dy are position relative to surface point
        dx = pos.x - pt.x
        dy = pos.y - pt.y
        #             Test if within aperture
        if dx*dx + dy*dy <= self.outerRadius*self.outerRadius:
            u = self.normal
        else:
            u = Blocked

        return SurfaceInteraction(self.type,pt,distance,pos,u,self.refractiveindex)


    def getParaxialInteraction(self,ray):
        p = self.getPoint()
        distance = p.z - ray.z
        height = ray.h + distance*ray.u

        dy = height - p.y
        if abs(dy) <= self.outerRadius:
            c = 0.0
        else :
            c = float("nan")

        return [self.type,distance,height,c,self.refractiveindex]
            
        


    def draw(self,option = None):
        """
        Draw the aperture wih two small chevrons to the current plot axis
        """
        p = self.getPoint()
        bar = self.outerRadius/5
        ytop = [p.y + self.outerRadius + bar/2 ,p.y + self.outerRadius + bar/2,\
                p.y + self.outerRadius,p.y + self.outerRadius + bar/2]
        zvals = [p.z - bar/2, p.z + bar/2, p.z , p.z - bar/2]
        ylower = [p.y - self.outerRadius - bar/2 ,p.y - self.outerRadius - bar/2,\
                  p.y - self.outerRadius,p.y - self.outerRadius - bar/2]

        return plot(zvals,ytop,"k",zvals,ylower,"k",lw=2.0)



class AnnularAperture(CircularAperture):
    """
    Class to  implement an annular aperture of fixed inner and outer radi rays outside the aperure will be blocked.

    :param pt: the plane reference point, (Defaults =(0.0,0.0,0.0))
    :type pt: Vector3d or float
    :param innerradius: inner radius of aperture (Default = 0.0)
    :type innerradius: float
    :param outerradius: outerradius of aperture (Default = 1.0)
    :type outerradius: float

    """
    
    def __init__(self,pt = 0.0, innerradius = 0.0, outerradius = 1.0):
        """
        Constuctor for a annular circular aperture
        
        """
        CircularAperture.__init__(self,pt,outerradius)
        self.innerRadius = innerradius

    def __str__(self):
        """
        Implements str()
        """
        return "spt: {0:s} inner: {1:7.3f} outer: {2:7.3f}".\
            format(str(self.point),self.innerRadius,self.outerRadius)
    

    def scale(self,a):
        """
        Scale aperture, will scale position and radius. if a < 0, then abs(a) is used.

        :param a the scale parameter
        :type a: float

        """
        CircularAperture.scale(self,a)
        self.innerRadius *= abs(a)
        return self

    #
    def getNormal(self,pos):
        """
        Method to get normal at specifed point, this will block rays outside the aperture
        """
        r = self.surfaceVector(pos)
        rsqr = r.absSquare()
        if rsqr <= self.outerRadius*self.outerRadius and rsqr >= self.innerRadius*self.innerRadius:
            return self.normal
        else:
            return Blocked



    def getSurfaceInteraction(self,ray):
        """
        Method to get the surface interaction information with a Ray
        
        :param ray: the input ray
        :type ray: :py:class:`optitcs.ray.IntensityRay`
        :return: `SurfaceInteraction`

        """
        pt = self.getPoint()
        distance = (pt.z - ray.position.z)/ray.director.z
        pos = ray.position.propagate(distance,ray.director)

        dx = pos.x - pt.x
        dy = pos.y - pt.y
        rsqr = dx*dx + dy*dy
        #
        #            Test that ray in is transparent region
        if rsqr <= self.outerRadius*self.outerRadius and rsqr >= self.innerRadius*self.innerRadius :
            u = self.normal
        else:
            u = Blocked
            
        return SurfaceInteraction(self.type,pt,distance,pos,u,self.refractiveindex)


    def getParaxialInteraction(self,ray):
        p = self.getPoint()
        distance = p.z - ray.z
        height = ray.h + distance*ray.u

        dy = height - p.y
        if abs(dy) <= self.outerRadius and abs(dy) >= self.innerRadius:
            c = 0.0
        else:
            c = float("nan")

        return [self.type,distance,height,c,self.refractiveindex]

    #
    def draw(self,option = None):
        """
        Draw the aperture wih two small chevrons and blocked centre.
        """
        p = self.getPoint()
        bar = self.outerRadius/5
        ytop = [p.y + self.outerRadius + bar/2 ,p.y + self.outerRadius + bar/2,\
                p.y + self.outerRadius,p.y + self.outerRadius + bar/2]
        zvals = [p.z - bar/2, p.z + bar/2, p.z , p.z - bar/2]
        ylower = [p.y - self.outerRadius - bar/2 ,p.y - self.outerRadius - bar/2,\
                  p.y - self.outerRadius,p.y - self.outerRadius - bar/2]
        yblock = [p.y + self.innerRadius, p.y - self.innerRadius]
        zvalb = [p.z,p.z]
        return plot(zvals,ytop,"k",zvals,ylower,"k",zvalb,yblock,"k",lw=2.0)

    #


class IrisAperture(CircularAperture):
    """
    Class to  implement a circular aperture of variable radius, rays outside the aperure
    will be blocked.

    :param pt: the plane reference point, (Defaults =(0.0,0.0,0.0))
    :type pt: Vector3d or float
    :param radius: radius of aperture (Default = 1.0)
    :type radius: float
    :param ratio: the faction aperture is opened, (Default = 1.0)
    
    """    
    
    def __init__(self,pt = 0.0, radius = 1.0, ratio = 1.0):
        """
        Constuctor for a iris aperture, being a variable radius circular aperture
       
        """
        CircularAperture.__init__(self,pt,radius)
        self.ratio = ratio                   # this typically controlled externally

    def __str__(self):
        """
        Implements str()
        """
        return "spt: {0:s} rad: {1:7.3f} ratio: {2:6.2f}".\
            format(str(self.point),self.outerRadius,self.ratio)

    def setRatio(self,ratio = 1.0):
        """
        Set the ratio of the aperture

        :param ratio: the new aperture ration, (Default = 1.0)
        :type ratio: float

        """
        self.ratio = float(ratio)
        return self
    
    def getRadius(self):
        """
        Get the radius inclduing the ratio.
        """
        return self.outerRadius*self.ratio
        
    #
    def getNormal(self,pos):
        """
        Method to get normal at specifed point,
        """
        r = self.surfaceVector(pos)
        radius = self.outerRadius*self.ratio          # Take into account the ratio
        if r.absSquare() <= radius*radius:
            return self.normal
        else:
            return Blocked

    def getSurfaceInteraction(self,ray):
        """
        Method to get back the surface interaction information with a Ray
        
        :param ray: the input ray
        :type ray: :py:class:`optitcs.ray.IntensityRay`
        :return: `SurfaceInteraction`
    
        """
        pt = self.getPoint()
        distance = (pt.z - ray.position.z)/ray.director.z
        pos = ray.position.propagate(distance,ray.director)

        dx = pos.x - pt.x
        dy = pos.y - pt.y
        
        radius = self.getRadius()        # Take into account the ratio
        if dx*dx + dy*dy <= radius*radius:
            u = self.normal
        else:
            u = Blocked

        return SurfaceInteraction(self.type,pt,distance,pos,u,self.refractiveindex)
        


    def getParaxialInteraction(self,ray):
        p = self.getPoint()
        distance = p.z - ray.z
        height = ray.h + distance*ray.u

        dy = height - p.y
        if abs(dy) <= self.outerRadius*self.ratio:
            c = 0.0
        else:
            c = float("nan")

        return [self.type,distance,height,c,self.refractiveindex]

    def draw(self,option = None):
        """
        Draw the aperture, same as aperure but with extra bars to mark current radius
        """
        p = self.getPoint()
        bar = self.outerRadius/5
        ytop = [p.y + self.outerRadius + bar/2 ,p.y + self.outerRadius + bar/2,\
                p.y + self.outerRadius,p.y + self.outerRadius + bar/2]
        zvals = [p.z - bar/2, p.z + bar/2, p.z , p.z - bar/2]
        ylower = [p.y - self.outerRadius - bar/2 ,p.y - self.outerRadius - bar/2,\
                  p.y - self.outerRadius,p.y - self.outerRadius - bar/2]

        radius = self.outerRadius*self.ratio
        ytopb = [p.y + radius , p.y + self.outerRadius]
        zvalb = [p.z,p.z]
        ylowerb = [p.y - radius, p.y - self.outerRadius]
        
        return plot(zvals,ytop,"k",zvals,ylower,"k",zvalb,ytopb,"k",zvalb,ylowerb,"k",lw=2.0)




class KnifeAperture(CircularAperture):
    """
    Class to give a circular aperture with a moveable "knife edge" or "wire" used in 
    optical testing.

    :param pt: the surface point 
    :type pt: Vector3d or float
    :param radius: radius of the aperture (default to 1.0)
    :type radius: float
    :param knife: distance of knife edge from axis (Default = 0.0)
    :type knife: float
    :param theta: angle of knife edge wrt to y-axis in radians (Default = 0.0)
    :type theta: float
    :param shift: shift along the z-axis from point position
    :type shift: float

    """
    def __init__(self,pt = 0.0, radius= 1.0 ,knife = 0.0, theta = 0.0, shift = 0.0):
        """
        
        """
        CircularAperture.__init__(self,pt,radius)
        self.setKnife(knife,theta,shift)
        self.setWire(False)

    def setKnife(self,knife = 0.0, theta = 0.0, shift = 0.0):
        """
        Set the knife parameters
        
        :param knife: distance of knife edge from axis (Default = 0.0)
        :type knife: float
        :param theta: angle of knife edge wrt to y-axis in radians (Default = 0.0)
        :type theta: float
        :param shift: distance along the z axis from the aperture position (Default = 0.0)
        :type shift: float

        """
        self.knife = float(knife)
        self.theta = float(theta)
        self.shift = Vector3d(0.0,0.0,shift)        # Hold shift as vector
        return self

    def setWire(self,wire = True, thickness = 0.01):
        """
        Switch to wire rather than knife
        
        :param wire: switch to wire mode (Default = True)
        :type wire: boolean
        :param thickness: thickness of wire (Default = 0.01 mm)
        :type thickness: float
        
        """
        self.wire = wire
        self.thickness = thickness
        return self

    def __str__(self):
        """
        Implement str()
        """
        return "spt: {0:s} rad: {1:7.3f} knife: {2:7.3f} angle: {3:7.3f} shift: {4:s}".\
            format(str(self.point),self.outerRadius,self.knife,self.theta,str(self.shift))
        

    def getSurfaceInteraction(self,ray):
        """
        Method to get back the surface interaction information with a Ray
        
        :param ray: the input ray
        :type ray: :py:class:`optitcs.ray.IntensityRay`
        :return: `SurfaceInteraction`

        """
        pt = self.getPoint() + self.shift   # Include local shift
        distance = (pt.z - ray.position.z)/ray.director.z
        pos = ray.position.propagate(distance,ray.director)

        #            dx and dy are position relative to surface point
        dx = pos.x - pt.x
        dy = pos.y - pt.y
        #             Test if within aperture
        if dx*dx + dy*dy <= self.outerRadius*self.outerRadius:

            #         Test wrt to knife edge
            r =  dy*math.cos(self.theta) - dx*math.sin(self.theta)
            if self.wire:
                if r > self.knife - self.thickness/2 and r < self.knife + self.thickness/2 :
                    u = Blocked
                else:
                    u = self.normal
            else:
                if r < self.knife:
                    u = self.normal
                else:
                    u = Blocked
        else:
            u = Blocked

        return SurfaceInteraction(self.type,pt,distance,pos,u,self.refractiveindex)


class ImagePlane(OpticalPlane):
    """
    Class to form an image plane (either input or output) of specific location and size.

    :param pt: the surface point,(Default to (0,0,0))
    :type pt: Vector3s or flrat
    :param xsize: the horizontal size, Defaults to 100.0mm
    :type xsize: float
    :param ysize: the vertical size, defaults to xsize (square)
    :type ysize: float

    Note the size does NOT affect the surface interaction, it only alters the .draw() method.
    """

    def __init__(self,pt = 0.0 ,xsize = 100.00, ysize = None):
        """
        Constructor for ImagePlane
       
        """
        OpticalPlane.__init__(self,pt)
        self.setSize(xsize,ysize)


    def __str__(self):
        """
        Implement str()
        """
        return OpticalPlane.__str__(self) + " xs: {0:6.3f} ys : {1:6.3f}".format(self.xsize,self.ysize)


    def scale(self,a):
        """
        Scale plane, scales both position and size, if a < 0,abs(a) is used.
        
        :param a: the scale parameter.
        :type a: float

        """
        a = abs(a)
        OpticalPlane.scale(self,a)
        self.xsize *= a
        self.ysize *= a
        return self

    def setSize(self,x,y = None):
        """
        Set the size of plane but without scaling the underlying OpticalPlane.

        :param x: xsize OR list of [xsize,ysize]
        :type x: float or list[float,float]
        :param y: ysize (if None defaults to xsize, so square)
        :type y: float or None

        """
        if isinstance(x,list) or isinstance(x,tuple) :
            self.xsize = float(x[0])
            self.ysize = float(x[1])
        else:
            if y == None:
                y = x
            self.xsize = float(x)
            self.ysize = float(y)
        return self

    def getSize(self):
        """
        Return the image size as truple
        """
        return self.xsize,self.ysize

    def draw(self,option = None):
        """
        Define Draw plane of correct height parameter.
        """
        p = self.getPoint()
        yvals =  [p.y + self.ysize/2, p.y - self.ysize/2]
        zvals =  [p.z,p.z]
        return plot(zvals,yvals,"k",lw=2.0)
        


class QuadricSurface(OpticalPlane):
    """
    Method to implement a quadric surface determined by curvature and epsilon, the qudric factor

    :param pos: the plane reference point
    :type pos: Vector3d or float
    :param curve: the curcature.
    :type curve: float
    :param epsilon: the quadric parameter
    :type epsilon: float
    :param radius: the maxradius
    :type radius: float
    :param index:  the Refratcive index. If present, surface is refracting, if None it is reflecting.
    :type index: RefractiveIndex or None

    """
                                                                    
    #
    def __init__(self, pos, curve, epsilon, radius , index = None):
        """
        Constructor
       
        """
        if index == None:
            OpticalPlane.__init__(self,pos,Reflecting)       # Reflecting
        else:
            OpticalPlane.__init__(self,pos,Refracting,index) # Refracting
        self.curvature = float(curve)
        self.epsilon = float(epsilon)
        self.maxRadius = float(radius)


    def __str__(self):
        """
        Implement str()
        """
        return "spt: {0:s} type: {1:d} c: {2:7.5f} e: {3:7.5f} rad: {4:6.4f} n: {5:s}".\
            format(str(self.point),self.type,self.curvature,self.epsilon,self.maxRadius,str(self.refractiveindex))
        
        
    def scale(self,a):
         """
         Scale surface, this will scale the position of the OpticalPlane, the
         maxradius and the curvature. If a < 0, then abs(a) will be applied to the plane
         and radius, but the sign of the curfature will be reversed.
         
         :param a:  the scale factor
         :type a: float

         """
         OpticalPlane.scale(self,abs(a))
         self.maxRadius *= abs(a)
         self.curvature /= a
         return self
     
    def getSourcePoint(self,pt_or_x,y = None,intensity = 1.0):
         """
         Get the SourcePoint for a specified point in the plane allowing for curvature.
        
         :param pt_or_x: Vector2d point in the plane or x component
         :type pt_or_x: Vector2d or float
         :param y: y component (Default = Null)
         :type y: float
         :param intensity:

         """

         if isinstance(pt_or_x,Vector2d):
            x = pt_or_x.x
            y = pt_or_x.y
         else:
            x = float(pt_or_x)
            y = float(y)
            
         c = self.curvature
         e = self.epsilon
         rsqr = x*x + y*y
        
         a = 1.0 - c*c*e*rsqr
         if a < 0.0 :
            raise ValueError("ray.QuadricSurface.pointin Plane: impossible surface: c: {0:8.5e} e: {1:8.5} r: {2:8.5e}"\
                             .format(c,e,rsqr))
         else:
            z =  c*rsqr/(1.0 + math.sqrt(a))
        
         refpt = self.getPoint()    # Ref point of plane
         x += refpt.x
         y += refpt.y
         z += refpt.z
        
         return SourcePoint([x,y,z],intensity)
    

    def getDistance(self,r,u):
        """
        Method to get the distance from r in direction u
        to the surface (not use in tratracing, use getSurfaceInteraction())

        :param r: the point in space.
        :type r: Vector3d
        :param u: the director at that point
        :type u: Unit3d
        :return: the distance to the surface as float

        """
        p = self.getPoint()     # Reference pt
        #
        d = (p.z - r.z)/u.z     # Distance to plane (direct calculation)
        x = r.x + d*u.x - p.x   # x/y Pt in plane
        y = r.y + d*u.y - p.y
        
        c = self.curvature
        eps = self.epsilon

        #
        #         form distance from plane to surface that will NOT fail for curve = 0    
        f = c*(x*x + y*y)
        g = u.z - c*(x*u.x + y*u.y)
        e = c*(1.0 + (eps - 1.0)*u.z*u.z)
        a = g*g - e*f
        #
        #          trap case of ray missing surface totally....... 
        if a < 0:
            return float("nan")           # Missed surface
        else:
            d += f/(g + math.copysign(math.sqrt(a),g))    # add extra distance
            return d

    #
    #
    def getNormal(self,r):
        """
        Method to get the surface normal at point r on the surface
        param r the point on the surface
        return Director, the surface normal
        """
        p = self.getPoint()

        x = r.x - p.x       # x/y in plane
        y = r.y - p.y

        c = self.curvature
        eps = self.epsilon

        if x*x + y*y > self.maxRadius*self.maxRadius:
            return Blocked        # Blocked by radius less
        else:
            x *= -c
            y *= -c
            z = 1.0 - c*eps*(r.z - p.z)
            return Unit3d(x,y,z)   # Will be auto normalsied
    #
    #
    #
    def edgePlane(self):
        """
        Get the plane at the edge of the surface,it will fail with ValueError
        if this is an impossible surface that has no edge....
        """
        c = self.curvature
        e = self.epsilon
        r = self.maxRadius
        
        a = 1.0 - c*c*e*r*r
        if a < 0.0 :
            raise ValueError("ray.QuadricSurface.edgePlane: impossible surface: c: {0:8.5e} e: {1:8.5} r: {2:8.5e}"\
                             .format(c,e,r))
        else:
            return c*r*r/(1.0 + math.sqrt(a))


    #
    #
    def eccentricity(self):
        """
        Method to get the eccentricity of the surface, not this uniques measures as it depend on the
        sign of epsilon, only of historical use.
        """
        if self.epsilon > 0:
            return math.sqrt(1.0 + self.epsilon)
        else:
            return math.sqrt(1.0 - self.epsilon)

    def getSurfaceInteraction(self,ray):
        """
        Method to get back the surface interaction information with a Ray
        
        :param ray: the input ray
        :type ray: :py:class:`optitcs.ray.IntensityRay`
        :return: `SurfaceInteraction`

        """
        p = self.getPoint()

        d = (p.z - ray.position.z)/ray.director.z     # Distance to plane (direct calculation)
        x = ray.position.x + d*ray.director.x - p.x   # x/y Pt in plane
        y = ray.position.y + d*ray.director.y - p.y

        
        c = self.curvature                            # local variable for c to simply typing.
        eps = self.epsilon

        #
        #         form distance from plane to surface that will NOT fail for curve = 0    
        f = c*(x*x + y*y)
        g = ray.director.z - c*(x*ray.director.x + y*ray.director.y)
        e = c*(1.0 + (eps - 1.0)*ray.director.z*ray.director.z)
        a = g*g - e*f
        #
        #          trap case of ray missing surface totally....... 
        if a > 0:
            d += f/(g + math.copysign(math.sqrt(a),g))     # add extra distance, retaining sign
        else:                                              # Flag distance as NaN
            return SurfaceInteraction(self.type,p,float("nan"),p,Blocked,self.refractiveindex)


        #          Get position of ray on surface
        pos = ray.position.propagate(d,ray.director)
        
        #          Now get surface normal at this point
        x = pos.x - p.x       # x/y in plane
        y = pos.y - p.y

        if x*x + y*y > self.maxRadius*self.maxRadius:
            u =  Blocked        # Blocked by max radius 
        else:
            x *= -c
            y *= -c
            z = 1.0 - c*eps*(pos.z - p.z)
            u =  Unit3d(x,y,z)   # Will be auto normalsied

        #
        #     Return list of information
        return SurfaceInteraction(self.type,p,d,pos,u,self.refractiveindex)

    #
    #
    def getParaxialInteraction(self,ray):
        """
        Get parxial interaction with surface abd return as list:
        type:   the surface type
        distance : distance to surface
        height :   ray height in surface
        curvature: curfature of the surface
        refrativeindex, refearive index on the right.
        """
        p = self.getPoint()
        distance = p.z - ray.z
        height = ray.h + distance*ray.u

        dy = ray.h - p.y
        if abs(dy) <= self.maxRadius:
            c = self.curvature
        else:
            c = float("nan")

        return [self.type,distance,height,c,self.refractiveindex]


    def draw(self,option = None):
        """
        Method to draw surface in matplotlinb plot which is returned
        """
        p = self.getPoint()       # Reference point
        yvals = []                # list of y value
        zvals = []                # list of z vals
        
        c = self.curvature
        e = self.epsilon

        for i in range(-SurfacePlotPoints,SurfacePlotPoints + 1):
            r = i*self.maxRadius/SurfacePlotPoints
            a = 1.0 - c*c*e*r*r
            if a > 0.0 :               # Ignore impossible surface
                z = c*r*r/(1.0 + math.sqrt(a))
                yvals.append(p.y + r)
                zvals.append(p.z + z)

        return plot(zvals,yvals,"k",lw=2.0)           # Plot in black
            
                
class SphericalSurface(QuadricSurface):
    """
    Method to implement a spherical  surface determined by curvature and,this is a QuadricSurface
    with epsilon = 1.0

    :param pos: the plane reference point
    :type pos: Vector3d or float
    :param curve: the curcature.
    :type curve: float
    :param radius: the maxradius
    :type radius: float
    :param index:  the Refratcive index. If present, surface is refracting, if None it is reflecting.
    :type index: RefractiveIndex or None

    """
    
    #
    def __init__(self,pos,curve,radius,index = None):
        """
        Constructor
        """
        QuadricSurface.__init__(self,pos,curve,1.0,radius,index)

        
    def __str__(self):
        """
        Implement str()
        """
        return "spt: {0:s} type: {1:d} c: {2:7.5f} rad: {3:6.4f} n: {4:s}".\
            format(str(self.point),self.type,self.curvature,self.maxRadius,str(self.refractiveindex))


class SphericalImagePlane(SphericalSurface):
    """
    Special spherical image plane to deal w
    """
    def __init__(self,pos,curve,maxRadius):
        """
        Constrcuor
        """
        SphericalSurface.__init__(self,pos,curve,maxRadius)
        self.type = Clear

    def __repr__(self):
        """
        Implement repr() for the surface
        """ 
        return "surface.SphericalImagePlane({0:s}, {1:8.6f}, {2:7.5f})".\
            format(str(self.point),self.curvature,self.maxRadius)
        

class ParabolicSurface(QuadricSurface):
    """
    Method to implement a parabolic surface determined by curvature and. This is a QudaricSurface
    with epsilon = 0.0
    
    :param pos: the plane reference point
    :type pos: Vector3d or float
    :param curve: the curcature.
    :type curve: float
    :param radius: the maxradius
    :type radius: float
    :param index:  the Refratcive index. If present, surface is refracting, if None it is reflecting.
    :type index: RefractiveIndex or None


    """
    
    
    def __init__(self,pos,curve,maxRadius,index = None):
        """
        Constructor
        """
        QuadricSurface.__init__(self,pos,curve,0.0,maxRadius,index)


    def __str__(self):
        """
        Implement str()
        """
        return "spt: {0:s} type: {1:d} c: {2:7.5f} rad: {3:6.4f} n: {4:s}".\
            format(str(self.point),self.type,self.curvature,self.maxRadius,str(self.refractiveindex))




