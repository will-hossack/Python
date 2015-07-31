"""
Set of classes to implement various types of optical surface.

Author: Will Hossack, The University of Edinburgh 
"""
from vector import Vector2d,Vector3d,Unit3d
import math
from matplotlib.pyplot import plot

"""
Define the three types of surface.
"""
Clear = 0
Refracting = 1
Reflecting = 2
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
    Class to return the surface interaction for a skew ray
    """
    
    def __init__(self,type,point,distance,position,normal,index):
        """
        param int type the surface
        param point Vector3d the surface reference point if global coordinate
        param distance float distance to surface
        param position interaction point on surface
        param normal surface normal at surface interaction point
        param index refreatcive index on image side of surface, may be null
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
            .format(self.type,str(self.point),self.distance,str(self.poistion))
        s += "Normal : {0:s}\nRefractiveIndex : {1:s}".format(str(self.normal),repr(self.refractiveindex))
        return s
        
#
#             Basic Surface 
#
class Surface(object):
    """
    Base class for an arbitrary surface, need to be extended to be useful.
    """
    #
    #
    def __init__(self,pt = 0.0 ,type = Clear, refindex = None):
        """
        Constructor to form a basic surface
        param pt Vector3d, the surface reference point, defaults to (0,0,0)
        param type int the surface type (default to Clear 0)
        param refindex the refractive index on the image side of the surface (default to None)
        """
        
        self.setPoint(pt)
        self.type = type
        #
        self.group = None
        self.refractiveindex = refindex
        
    #
    #
    def __repr__(self):
        """
        Implement the repr() method
        """
        return "Surface({0:s} , {1:d})".format(str(self.point),self.type)
        

    def setPoint(self,z_or_v = 0.0):
        """
        Method to set the point in a consistent way for all surfaces
        param z_or_v can be
        Vector3d or Position
        list or tuple with 3 compoents
        Vector2d will give (x,y,0.0)
        float or int will set to (0,0,z)
        default will set to (0,0,0)
                           
        """
        if isinstance(z_or_v,Vector3d) or isinstance(z_or_v,list) or isinstance(z_or_v,tuple):
            self.point = Vector3d(z_or_v)
        elif isinstance(z_or_v,Vector2d):
            self.point = Vector3d(z_or_v.x,z_or_v.y,0.0)
        elif isinstance(z_or_v,float) or isinstance(z_or_v,int):
            self.point = Vector3d(0,0,z_or_v)
        else:
            raise TypeError("surface.Surface.setPoint: called with unknown type")
        
        return self


    def scale(self,a):
        """
        Scale surface
        """
        self.point *= a
        return self
    #          
    #
    def getPoint(self):
        """
        Method to get the surface point in global coordinates taking account
        that the Surface my belong to an OpticalGroup
        return Vector3d the reference point
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
        """
        print("Surface.getNormal needs to be be defined")
        return Unit3d()

    #
    #
    def draw(self):
        print("surface.Surface.draw error Need to implement draw")

#
#
class FlatSurface(Surface):
    """
    Class to implement an Flat Surface with specifed surface normal
    """
    #
    #
    def __init__(self,pt = None, normal = Unit3d(0,0,1), type = Clear,refindex = None): 
        """
        Constructor
        param pt Vector3d, reference point, (defaults to (0,0,0))
        param normal Unit3d, the surface normal
        param type int (defaults to Clear)
        param refindex the refractive index (defaults to None)
        """
        Surface.__init__(self,pt,type,refindex)
        self.normal = Unit3d(normal)
        self.curvature = 0.0
    #
    #
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
        param r the position
        param u the director
        """
        p = self.getPoint()      # Reference point
        dv = p - r
        dist = self.normal.dot(dv)/self.normal.dot(u)
        return dist
        

    def getSurfaceInteraction(self,ray):
        """
        Method to get back the surface interaction information for a ray.
        type:     surface type
        point:    the surface point in global coordinates
        distance: distance from current ray position to surface
        pos :     Position, intration point with surface
        norm:     surface normal at that point
        refrative : refrative index (if refracting surface)
        """
        pt = self.getPoint()
        dv = pt - ray.position
        distance = self.normal.dot(dv)/self.normal.dot(ray.director)
        pos = ray.position.propagate(distance,ray.director)
        
        return SurfaceInteraction(self.type,pt,distance,pos,self.normal,self.refractiveindex)
       
#
#
class OpticalPlane(FlatSurface):
    """
    Class to implement a flat optical place normal to the optical, this is simpler and
    version of FlatSurface with a fixed normal axis.
    """
    #
    #
    def __init__(self,pt = 0.0 ,type = Clear, refindex = None):
        """
        Conststrucor for Optical plane
        param pt the surface point, (defaults to (0.0,0.0,0.0)
        param type surface type, defaults to Clear (0)
        param refindex refrative index of right of surface, defaults to None
        """
        FlatSurface.__init__(self,pt,Unit3d(0,0,1),type,refindex)
    
    #
    #
    def __repr__(self):
        """
        Implement the repr() method
        """
        return "surface.OpticalPlane({0:s} , {1:d})".format(str(self.point),self.type)
    
    #
    def surfaceVector(self,pos):
        """
        Method to get the Vector2d from the reference point to the specifed position, assumed to
        be in the plane.
        """
        p = pos - self.getPoint()
        return Vector2d(p.x,p.y)



    def getSurfaceInteraction(self,ray):
        """
        Method to get back the surface interaction information for a ray.
        type:     surface type
        point:    surface reference point in global coorinates
        distance: distance from current ray position to surface
        pos :     Position, intration point with surface
        norm:     surface normal at that point
        refrative : refrative index (if refracting surface)
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
        
        

    def draw(self,height=10.0):
        """
        Define Draw with exita height parameter.
        """
        p = self.getPoint()
        yvals =  [p.y + height, p.y - height]
        zvals =  [p.z,p.z]
        return plot(zvals,yvals,"k",lw=2.0)
    
class CircularAperture(OpticalPlane):
    """
    Class to  implement a circular aperture of fixed radius, rays outside the aperure
    will be blocked.
    """
    #
    #
    def __init__(self,pt = 0.0 , radius = 1.0):
        """
        Constuctor for a circular aperture
        param pt Position, the plane reference point, dfaults to (0.0,0.0,0.0)
        param radius the radius of the aperture (deafults = 1.0)
        """
        OpticalPlane.__init__(self,pt)
        self.outerRadius = radius
        self.maxRadius = radius                          # add max radius to standardise external refs

    def __repr__(self):
        """
        Implement the repr()
        """
        return "opticalsurface.CircularAperture({0:s} , {1:8.5f})".\
            format(str(self.point),self.outerRadius)

    def scale(self,a):
        """
        Scale aperture, will scale position and radius. if a < 0, then abs(a) is used.
        param a the scale parameter.
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
        Method to get back the surface interaction information for a ray
        type:     surface type
        point:    surface reference point in global coordinates
        distance: distance from current ray position to surface
        pos :     Position, intration point with surface
        norm:     surface normal at that point
        refrative : refrative index (if refracting surface)
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
            
        

    #
    def draw(self):
        """
        Draw the aperture wih two small chevrons
        """
        p = self.getPoint()
        bar = self.outerRadius/5
        ytop = [p.y + self.outerRadius + bar/2 ,p.y + self.outerRadius + bar/2,\
                p.y + self.outerRadius,p.y + self.outerRadius + bar/2]
        zvals = [p.z - bar/2, p.z + bar/2, p.z , p.z - bar/2]
        ylower = [p.y - self.outerRadius - bar/2 ,p.y - self.outerRadius - bar/2,\
                  p.y - self.outerRadius,p.y - self.outerRadius - bar/2]

        return plot(zvals,ytop,"k",zvals,ylower,"k",lw=2.0)

    #

class AnnularAperture(CircularAperture):
    """
    Class to  implement an annular aperture of fixed inner and outer radius.
    """
    #
    #
    def __init__(self,pt = 0.0, innerradius = 0.0, outerradius = 1.0):
        """
        Constuctor for a circular aperture
        param pt Position, the plane reference point
        param radius the radius of the aperture (deafults = 1.0)
        """
        CircularAperture.__init__(self,pt,outerradius)
        self.innerRadius = innerradius

    def __repr__(self):
        """
        Implement the repr()
        """
        return "opticalsurface.AnnularAperture({0:s} , {1:8.5f} , {2:8.5f})".\
            format(str(self.point),self.innerRadius,self.outerRadius)

    def scale(self,a):
        """
        Scale aperture, will scale position and radius. if a < 0, then abs(a) is used.
        param a the scale parameter.
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
        Method to get back the surface interaction information for a ray
        Returns the list
        type:     surface type
        distance: distance from current ray position to surface
        pos :     Position, intration point with surface
        norm:     surface normal at that point
        refrative : refrative index (if refracting surface)
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
        if abs(dy) <= self.outerRadius and abs(y) >= self.innerRadius:
            c = 0.0
        else:
            c = float("nan")

        return [self.type,distance,height,c,self.refractiveindex]

    #
    def draw(self):
        """
        Draw the aperture wih two small chevrons
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
    Class to implement a variable radius circular aperture.
    """    
    #
    #
    def __init__(self,pt = 0.0, radius = 1.0, ratio = 1.0):
        """
        Constuctor for a iris aperture, being a variable radius circular aperture
        param pt Position, the plane reference point
        param radius float the radius of the aperture (deafults = 1.0)
        parar ratio float fration that aperture is open
        """
        CircularAperture.__init__(self,pt,radius)
        self.ratio = ratio                   # this typically controlled externally

    def __repr__(self):
        """
        Implement the repr()
        """
        return "surface.IrisAperture({0:s} , {1:8.5f}, {2:8.5f})".\
            format(str(self.point),self.outerRadius,self.ratio)
        
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
        Method to get back the surafce information information for a ray
        Returns the list
        type:     surface type
        distance: distance from current ray position to surface
        pos :     Position, intration point with surface
        norm:     surface normal at that point
        refrative : refrative index (if refracting surface)
        """
        pt = self.getPoint()
        distance = (pt.z - ray.position.z)/ray.director.z
        pos = ray.position.propagate(distance,ray.director)

        dx = pos.x - pt.x
        dy = pos.y - pt.y
        
        radius = self.outerRadius*self.ratio          # Take into account the ratio
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

    def draw(self):
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


class ImagePlane(OpticalPlane):
    """
    Class to form an image plane (either input or output) of specific location and size.
    """

    def __init__(self,pt = 0.0 ,xsize = 36.00, ysize = 24.0):
        """
        Constructor for ImagePlane
        param pt the surface point, defaults to (0,0,0)
        param xsize the horizontal size, defaults to 36.00 (horizontal size of 35 mm film)
        param ysize the vertical size, defaults to 24.00 (horizontal size of 35 mm film)

        Note the size does NOT affect the surface interaction, it only alters the .draw() method.
        """
        OpticalPlane.__init__(self,pt)
        self.xsize = xsize
        self.ysize = ysize


    def __str__(self):
        """
        Implement str()
        """
        return "({0:s},{1:8.5f},{2:8.5f})".format(str(self.point),self.xsize,self.ysize)

    def __repr__(self):
        """
        Implement repr()
        """
        return "surface.ImagePlane" + str(self)



    def scale(self,a):
        """
        Scale plane, scales both position and size, if a < 0,abs(a) is used.
        param a the scale parameter.
        """
        a = abs(a)
        OpticalPlane.scale(self,a)
        self.xsize *= a
        self.ysize *= a
        return self

    def setSize(self,x,y):
        """
        Set the size of plane but without scaling the underlying OpticalPlane.
        """
        self.xsize = float(x)
        self.ysize = float(y)
        return self


    def draw(self):
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
    """
                                                                    
    #
    def __init__(self,pos,curve,epsilon,maxRadius,index = None):
        """
        Constructor
        param pos the plane reference point
        param curve the curcature
        param epsilon the quadric parameter
        param maxRadius the maxradius
        index the Refratcive index. If present, surface is refracting, if None 
        it is reflecting.
        """
        if index == None:
            OpticalPlane.__init__(self,pos,Reflecting)       # Reflecting
        else:
            OpticalPlane.__init__(self,pos,Refracting,index) # Refracting
        self.curvature = curve
        self.epsilon = epsilon
        self.maxRadius = maxRadius

        
    def scale(self,a):
         """
         Scale surface, this will scale the position of the OpticalPlane, the
         maxradius and the curvature. If a < 0, then abs(a) will be applied to the plane
         and radius, but the sign of the curfature will be reversed.
         param a the scale factor
         """
         OpticalPlane.scale(self,abs(a))
         self.maxRadius *= abs(a)
         self.curvature /= a
         return self
    

    def getDistance(self,r,u):
        """
        Method to get the distance from r in direction u
        to the surface
        param r the point in space
        param u the director at that point
        return the distance to the surface
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
        Method to get back the surface information information for a ray
        Returns the list
        type:     surface type
        distance: distance from current ray position to surface
        pos :     Position, intration point with surface
        norm:     surface normal at that point
        refrative : refrative index (if refracting surface)
        """
        p = self.getPoint()

        d = (p.z - ray.position.z)/ray.director.z     # Distance to plane (direct calculation)
        x = ray.position.x + d*ray.director.x - p.x   # x/y Pt in plane
        y = ray.position.y + d*ray.director.y - p.y

        
        c = self.curvature                            # local variable to 
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


    def draw(self):
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
    """
    
    #
    def __init__(self,pos,curve,maxRadius,index = None):
        """
        Constructor
        param pos the plane reference point, may be Position of float
        param curve the curvature
        param maxRadius the maximum radius
        index the Refractive index
        """
        QuadricSurface.__init__(self,pos,curve,1.0,maxRadius,index)

    #
    #
    def __repr__(self):
        """
        Implement repr() for the surface
        """ 
        return "surface.SphericalSurface({0:s}, {1:8.6f}, {2:7.5f}, {3:s})".\
            format(str(self.point),self.curvature,self.maxRadius,repr(self.refractiveindex))



class SphericalImagePlane(SphericalSurface):
    """
    Special spherical image plane to deal with the humen eye model
    """
    def __init__(self,pos,curve,maxRadius):
        """
        Paramteters of surface
        param pos position, either float or Position
        param curve the curvature
        param maxRadius maxradius of the surface
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
    """
    
    #
    def __init__(self,pos,curve,maxRadius,index = None):
        """
        Constructor
        param pos the plane reference point
        param curve the curvature
        param maxRadius the maxradius
        index the Refratcive index
        """
        QuadricSurface.__init__(self,pos,curve,0.0,maxRadius,index)          

    #      
    #
    #
    def __repr__(self):
        """
        Implement repr() for the surface
        """ 
        return "surface.ParabolicSurface({0:s}, {1:8.6f}, {2:7.5f}, {3:s})".\
            format(str(self.point),self.curvature,self.maxRadius,repr(self.refractiveindex))



