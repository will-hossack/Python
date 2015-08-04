"""
Set of classes for analysis of optical systems
"""
import ray
from surface import OpticalPlane,ImagePlane,SurfaceInteraction,SphericalSurface,KnifeEdgeAperture
from wavelength import Default,WavelengthColour
from vector import Vector2d,Vector3d,Unit3d,Angle
import matplotlib.pyplot as plt
import numpy as np
import math
import array

class TargetPlane(ImagePlane):
    """
    For a target plane, being at ImagePlane with target points or various types
    """
    #
    def __init__(self,pt = 0.0 ,xsize = 36.00, ysize = 24.0):
        """
        Constuctor with
        param pt Poistion or float, the plane position
        param xsize float xsize of plane (defaults to 36mm)
        param ysize float ysize of plane (defaults to 24 mm)
        """
        ImagePlane.__init__(self,pt,xsize,ysize)
        self.targets = []                 # List of targets
        self.wavelength = Default

    def __repr__(self):
        """
        Implement repr()
        """
        return "analysis.TargetPlane" + str(self)

    #
    #
    def add(self,t,key = "global"):
        """
        Add a target
        """
        if isinstance(t,ray.RayPencil):    # Deal with ray pencil
            for r in t:
                if r:
                    self.add(r)
                    
        elif isinstance(t,ray.IntensityRay): # Deal with ray
            pt = self.getPoint()
            self.wavelength = t.wavelength
            v = t.pointInPlane(pt.z)
            self.add(Vector3d(v.x,v.y,pt.z))
        else:
            if key.startswith("g") :       # Use global coordinates
                self.targets.append(t)     # Just add
            elif key.startswith("l"):                         # apply local
                pt = self.getPoint()
                p = Vector3d(t.x - pt.x, t.y - pt.y, pt.z)
                self.targets.append(p)
            else:
                raise TypeError("analysis.TargetPlane.add illegal key")

    #
    def addGrid(self,xn,yn = 0,radius = float("inf")):
        """
        Fill the TargetPlace with a set regular set of targets
        param xn number of targets across horizontal
        param yn numbers of tarets across if 0 or negative, n will
        be set to that the targets are on a square grid
        param radius, float, have targets only in a masked of 
        specified radius.

        Note: number of targets will be rounded to ensure an
        odd number across array so there will always be one
        at (0,0)
        """
        pt = self.getPoint()

        dx = 0.5*self.xsize/(xn/2 + 0.1)
        if yn > 0 :
            dy = 0.5*self.ysize/(yn/2 + 1 + 0.1)
        else:
            dy = dx
            yn = int(round(self.ysize/dy))

        for j in range(-yn/2,yn/2+1):
            for i in range(-xn/2,xn/2+1):
                y = dy*j
                x = dy*i
                if x*x + y*y < radius*radius and \
                   abs(x) < self.xsize/2 and abs(y) < self.ysize/2:
                     self.add(Vector3d(x + pt.x ,y + pt.y ,pt.z),"global")
        return self

    def rayPencil(self,pt_or_u,wave = Default, intensity = 1.0):
        """
        Get an intesnity RayPenci, one ray from each target
        param pt_or_u Vector3d or Unit3d, of Position each 
        ray will pass through this point, if Director, this
        is direction of each ray.
        param wave wavelength, (defaults of Default)
        param intensity intensity, (defaults to 1.0)
        """
        pt = self.getPoint()
        pencil = ray.RayPencil()
        for t in self.targets:
            #                Start position of ray
            pos = Vector3d(pt.x + t.x, pt.y + t.y, pt.z)
            if isinstance(pt_or_u,Unit3d):
                u = Unit3d(pt_or_u)
            else:
                u = Unit3d(pt_or_u - pos)
            r = ray.IntensityRay(pos,u,wave,intensity)
            pencil.append(r)
        return pencil

        
    def draw(self):
        """
        Draw the plane in MatLibPlot
        """
        pt = self.getPoint()
        xframe = [pt.x - self.xsize/2, pt.x + self.xsize/2,\
                  pt.x + self.xsize/2, pt.x - self.xsize/2,\
                  pt.x - self.xsize/2]
        yframe = [pt.y + self.ysize/2, pt.y + self.ysize/2,\
                  pt.y - self.ysize/2, pt.y - self.ysize/2,
                  pt.y + self.ysize/2]

        col = WavelengthColour(self.wavelength).hexString()

        xpt = []
        ypt = []
        
        for t in self.targets:
            xpt.append(t.x)
            ypt.append(t.y)

        plt.plot(xframe,yframe,"k")
        plt.plot(xpt,ypt,linestyle='none',color=col,marker='x')


class OpticalImage(ImagePlane):
    """
    Class to hold an image in a plane with a sampling grid. The actual
    image is held in a nmpy array
    """
    
    def __init__(self,pt = None,xsize = 200.0, ysize = 200.0, xpixel_or_im = 256, ypixel = 256):
        """
        Form the OpticalImage with either blank array of nmpy image array
        param pt the plane point (default = None,(0,0,0))
        param xsize xsize of plane in mm (default = 200.0mm)
        param ysize ysize of plane in mm (default = 200.0mm)
        param xpixel_or_im x-pixel size of image (default = 256) OR nmpy array of floats
        param ypixel y-pixel size of of image (default = 256)
        """
        
        ImagePlane.__init__(self,pt,xsize,ysize)
        if isinstance(xpixel_or_im,int):
            self.image = np.zeros((xpixel_or_im,ypixel),dtype = float)
        else:
            self.image = xpixel_or_im
        self.xpixel,self.ypixel = self.image.shape       # set xpixel and ypixel from image data


    def __str__(self):
        """
        Implement str()
        """
        return "({0:s},{1:8.5f},{2:8.5f},{3:d},{4:d})".format(str(self.point),self.xsize,self.ysize,self.xpixel,self.ypixel)

    def __repr__(self):
        """
        Implement repr()
        """
        return "analysis.OpticalImage" + str(self)



    def getSource(self,i,j):
        """
        Get pixel as i,j as a SourcePoint
        param i the x pixel location
        param y the y pixel location
        return SourcePoint giving x,y,z and intensity of pixel in global coordinates.
        """
        pt = self.getPoint()              # Reference point
        x = self.xsize*(float(i)/self.xpixel - 0.5) + pt.x
        y = self.ysize*(float(j)/self.ypixel - 0.5) + pt.y
        pos = Vector3d(x,y,pt.z)
        return ray.SourcePoint(pos,self.image[i,j])


    def getSurfaceInteraction(self,r):
        """
        Method to get back the surface interaction information for a ray
        Returns the list
        type:     surface type
        distance: distance from current ray position to surface
        pos :     Position, intration point with surface
        norm:     surface normal at that point
        refrative : refrative index (if refracting surface)

        This also add the ray intensity to the cloeset pixel
        """
        
        #       get interaction with super class
        info = ImagePlane.getSurfaceInteraction(self,r)
        
        #            Add ray to pixel
        i = int(round(self.xpixel*(info.position.x + self.xsize/2 - info.point.x)/self.xsize))
        j = int(round(self.ypixel*(info.position.y + self.ysize/2 - info.point.y)/self.ysize))

        #          Check if pixel is in image (note due to distrortions it may not be)
        if i >= 0 and i < self.xpixel and j >=0 and j < self.ypixel:
            self.image[i,j] += r.intensity               # Add it to the image

        #          Retun info to calling object
        return info 


    def getRayPencil(self,ca,i,j,nrays = 5,wave = Default):
        """
        Method to get a RayPencil from the i,j image pixel. 
        param ca circular apereture (or lens) to fill
        param i,j the pixel coordinates
        nrays number of ray across radius (default = 5)
        param wave wavelength of rays (default = Default)
        
        Note will return None of the pixel intensity is 0.0
        """
        source = self.getSource(i,j)
        if source.intensity == 0.0:
            return None
        else:
            return ray.RayPencil().addSourceBeam(ca,source,"array",nrays,wave)

    def getImage(self,lens, mag, nrays = 5, wave = Default, design = Default):
        """
        Method to get the image of OpticalPlane, which is assumes to holds an image,
        through lens as spefied mag with specifed wavelengths.
        param lens the lens system
        param mag the imaging magnifications
        param nrays number of rays on radius
        param wave wavelength of imaging (to do the actual tracing)
        param design wavelength of the design (used to set the object/image location)
        return OpticalImage with same pixel resolution as the object

        Note: The current OpticalImage will be moved to the object plane, this method will
        be slow for large images or complex lenses.
        """
        op,ip = lens.planePair(mag,design)    # get object / image locations
        self.setPoint(op)                     # Move object plane to corerct location

        ximage = abs(mag)*self.xsize          # Size of image plane
        yimage = abs(mag)*self.ysize

        image = OpticalImage(ip,ximage,yimage,self.xpixel,self.ypixel)  # Make image (set to zero)

        #
        #            Go through each pixel in turn and progate it.
        #
        for j in range(0,self.ypixel):
            for i in range(0,self.xpixel):
                pencil = self.getRayPencil(lens, i, j, nrays, wave)
                if pencil != None:                   # Will be None if pixel intensity is zero, so don't bother
                    pencil *= lens
                    pencil *= image

        return image                                 # Return the image

    def addTestGrid(self, xgap = 10, ygap = None ):
        """
        Method to add a test grid being a grid of one pixel wide
        in a grid pattern.
        param xgap gap in x directions between pixels (defaults to 10)
        param ygap gap in y directions between pixels, (defaults to xgap)
        """
        if ygap == None:
            ygap = xgap
        
        xw = xgap*(self.xpixel//xgap)
        yw = ygap*(self.ypixel//ygap)

        xs = (self.xpixel - xw)//2
        ys = (self.ypixel - yw)//2

        for j in range(0,self.ypixel):
            for i in range(0,self.xpixel):
                if j%ygap == ys or i%xgap == xs:
                    self.image[i,j] = 1.0

        return self

    def draw(self):
        """
        Display the image via imshow with gray comlour map and correct etent
        """
        return plt.imshow(self.image,cmap=plt.cm.gray,extent=(-self.xsize/2,self.xsize/2,-self.ysize/2,self.ysize/2))
        



class CurvedOpticalImage(OpticalImage,SphericalSurface):
    """
    Class to hold an curved image in a plane with a sampling grid. The actual
    image is held in a nmpy array
    """
    
    def __init__(self,pt = None,curve = 0.0, xsize = 200.0, ysize = 200.0, xpixel_or_im = 256, ypixel = 256):
        """
        Form the OpticalImage with either blank array of nmpy image array
        param pt the plane point (default = None,(0,0,0))
        param curve the curvature of the plane (taken as spherical)
        param xsize the x size (default = 200)
        param ysize the y size (default  = 200)
        param xpixel_or_im x-pixel size of image (default = 256) OR nmpy array of floats
        """
        

        OpticalImage.__init__(self,pt,xsize,ysize,xpixel_or_im,ypixel)
        SphericalSurface.__init__(self,pt,curve,math.sqrt(0.25*(xsize*xsize + ysize*ysize)))


    def __str__(self):
        """
        Implement str()
        """
        return "({0:s},{1:8.5f},{2:8.5f},{3:8.5f},{4:d},{5:d})".\
            format(str(self.point),self.curvature,self.xsize,self.ysize,self.xpixel,self.ypixel)

    def __repr__(self):
        """
        Implement repr()
        """
        return "analysis.CurvedOpticalImage" + str(self)



    def getSource(self,i,j):
        """
        Get pixel as i,j as a SourcePoint
        param i the x pixel location
        param y the y pixel location
        return SourcePoint giving x,y,z and intensity of pixel in global coordinates.
        """
        pt = self.getPoint()              # Reference point
        x = self.xsize*(float(i)/self.xpixel - 0.5)
        y = self.ysize*(float(j)/self.ypixel - 0.5)
        rsqr = x*x + y*y
        a = 1.0 - self.curvature*self.curvature*rsqr

        
        if a < 0.0 :
            raise ValueError("CurvedOpticalImage.getSource: impossible surface: c: {0:8.5e} e: {1:8.5} r: {2:8.5e}"\
                             .format(c,e,r))
        else:
            z =  c*r*r/(1.0 + math.sqrt(a))

        pos = Vector3d(x + pt.x ,y + pt.y ,z + pt.z)
        return ray.SourcePoint(pos,self.image[i,j])


    def getSurfaceInteraction(self,r):
        """
        Method to get back the surface interaction information for a ray
        Returns the list
        type:     surface type
        distance: distance from current ray position to surface
        pos :     Position, intration point with surface
        norm:     surface normal at that point
        refrative : refrative index (if refracting surface)

        This also add the ray intensity to the cloeset pixel
        """

        info = SphericalSurface.getSurfaceInteraction(self,r)

        


        #     Return list of information
        return SurfaceInteraction(self.type,d,pos,u,self.refractiveindex)


class Moments(object):
    """
    Class to implment two dimensional moments
    """

    def __init__(self,order = 2):
        """
        Constructor with order only, defaults to 2
        """
        self.order = order
        self.xdim = order + 1
        #self.moment = array.array('d',[0.0]*(self.__index(0,order) + 1))
        self.moment = array.array('d',[0.0]*(self.xdim*self.xdim))
        self.points = 0
       
    
    #
    def get(self,m,n):
        """
        Method to get a momment value
        """
        #    return self.moment[self.__index(m,n)]
        return self.moment[m*self.xdim + n]

    #
    def addPoint(self,p,value):
        """
        Method to add a point
        """
        self.points += 1
        for k in range(0,self.order + 1):
            for n in range(0,k + 1):
                m = k - n
                # i = self.__index(m,n)
                v = value*math.pow(p.x,m)*math.pow(p.y,n)
                self.moment[m*self.xdim + n] += v
        return self

    #
    def addRay(self,ray,plane):
        """
        Method add ray or raypencil in speficied place
        """
        if isinstance(ray,list): # If list add all the valid rays
            for r in ray:
                if r:
                    self.addRay(r,plane)
        else:
            if ray:
                v = ray.pointInPlane(plane)
                self.addPoint(v,ray.intensity)
        return self


    #
    def centroid(self):
        m = self.get(0,0)
        return Vector2d(self.get(1,0)/m,self.get(0,1)/m)

    def radius(self):
        c = self.centroid()
        m = self.get(0,0)
        r = (self.get(2,0) + self.get(0,2))/m - c.absSquare()
        return r


    def __index(self,m,n):
        """
        Internal method to look up moments in list
        """
        ord = m + n
        return (ord*(ord + 1)/2 + n)


class MomentsFixed(object):
    """
    Class to implment two dimensional moments
    """

    def __init__(self):
        """
        Constructor with order only, defaults to 2
        """
        self.m00 = 0.0
        self.m10 = 0.0
        self.m01 = 0.0
        self.m20 = 0.0
        self.m11 = 0.0
        self.m02 = 0.0
        self.points = 0
       
    
    #
    def addPoint(self,p,value):
        """
        Method to add a point
        """
        self.points += 1
        self.m00 += value
        self.m10 += value*p.x
        self.m01 += value*p.y
        self.m20 += value*p.x*p.x
        self.m11 += value*p.x*p.y
        self.m02 += value*p.y*p.y
        return self

    #
    def addRay(self,ray,plane):
        """
        Method add ray or raypencil in speficied place
        """
        if isinstance(ray,list): # If list add all the valid rays
            for r in ray:
                if r:
                    self.addRay(r,plane)
        else:
            if ray:
                v = ray.pointInPlane(plane)
                self.addPoint(v,ray.intensity)
        return self


    #
    def centroid(self):
        
        return Vector2d(self.m10/self.m00,self.m01/self.m00)

    def radius(self):
        c = self.centroid()
        r = (self.m20 + self.m02)/self.m00 - c.absSquare()
        return r


class Psf(Vector3d):
    """
    Class to represent a geometric Psf
    """

    #
    #
    def __init__(self,pos = Vector3d(), intensity = 1.0, a = 1.0, b = None, alpha = 0.0):
        """
        Form Psf with parameters
        param pos the position in 3d
        param a radius or major axis
        param b (default None), 
        param alpha angle of ellipse
        """ 
        Vector3d.__init__(self,pos)
        self.intensity = 1.0
        self.major = a
        if b == None:
            self.minor = self.major
        else:
            self.minor = b
        self.alpha = alpha
        
    #
    #
    def __repr__(self):
        """
        Implement repr()
        """
        return "analysis.Psf({0:s},{1:7.5f},{2:8.5e},{3:8.4e},{4:8.4e})".\
            format(str(self),self.intensity,self.major,\
                   self.minor,self.alpha)


    def eccentricity(self):
        """
        Eccenricity of the ellipse
        """
        if a != 0.0 :
            return math.sqrt(1.0 - (self.minor*self.minor)/(self.major*self.major))

    def area(self):
        """
        Area of PSF
        """
        return math.pi*self.major*self.minor
        

    def setWithRays(self,pencil,plane):
        """
        Set PSF from RayPencil in specified plane
        """
        #          Form the moments
        mom = MomentsFixed().addRay(pencil,plane)
        c = mom.centroid()          # Get the centre
        self.set(c.x,c.y,plane)     # Cet position in 3d

        #          Get the moment is local variables
        m00 = mom.m00
        u20 = mom.m20/m00 - c.x*c.x
        u02 = mom.m02/m00 - c.y*c.y
        u11 = mom.m11/m00 - c.x*c.y

        #          Hand code eigen values of covariance marix
        p = u20 + u02
        q = math.sqrt(4.0*u11*u11 + (u20 - u02)**2)

        #          Set the values
        self.intensity = m00
        self.major = math.sqrt(p + q)       # Major axis
        self.minor = math.sqrt(p - q)       # Minor axis
        self.alpha = 0.5*math.atan2(2.0*u11 , (u20 - u02))

        return self
    #
    #
    def optimalArea(self,pencil,plane):
        """
        Method to find the optimal area PSF from a raypencil starting
        as the guess locaion plane """
        
        #    get avarege plane of pencil
        n = 0
        zave = 0.0
        for r in pencil:
            if r:
                zave += r.position.z
                n += 1

        zave /= float(n)

        #            set self with initial condition
        psf = self.setWithRays(pencil,plane)
        area = psf.area()
        zp = plane
        wave = pencil[0].wavelength/1000   # Wavelengh in mm
        delta = -(plane - zave)/100     # Guess at shift
        deltaReduced = True

        #            Loop looking for best
        while True:
            self = Psf().setWithRays(pencil, zp + delta)
            na = self.area()
            if na < area :             # Accept
                zp += delta
                area = na
                deltaReduced = False
                if abs(delta) < wave:
                    break
            else:                      # Dont accept
                if deltaReduced:       # Delta not reduced las time
                    delta = -delta
                    deltaReduced = False
                else:
                    delta *= 0.25      # Reduce delta
                    deltaReduced = True
        #
        #        Found best PSF, in self, just return 
                    
        return self
        
    #
    #
    def draw(self):
        
        n = 20
        dtheta = 2*math.pi/n
        xval = []
        yval = []

        for i in range(0,n + 1):
            theta = i*dtheta
            v = Vector2d(self.major*math.cos(theta),self.minor*math.sin(theta))
            v.rotate(-self.alpha)
            v += Vector2d(self.x,self.y)
            xval.append(v.x)
            yval.append(v.y)

        return plt.plot(xval,yval,"b")


#
#
def aberrationPlot(lens,angle,wave = Default, design =  Default, nrays = 50):
    """
    Method to form the three standard transverse aberrations plots to
    a plt plot.
    param lens the OpticalGroup holding the lens
    paran angle field angle to for the plot
    param wave wavelength of the plot (default to Default)
    param design design wavelength *default to Default)
    param nrays number of rags to trace (default = 50)
    return the three plots as a [] with suitable labels
    """
    
    if isinstance(angle,float):
        u = Unit3d(Angle(angle))                     # direction of beam
    else:
        u = Unit3d(angle)

    ref = lens.imagePoint(u,design)               # Get image point at design wavelength
    ip = OpticalPlane(ref.z)                      # Make back focal plane to proagate to 

    entrance = lens.entranceAperture()
    dr = entrance.maxRadius/(nrays + 0.1)         # Ray separation
    
    rvals = []                # Radius values
    mvals = []                # Meridional
    svalsx = []               # Sagittal x
    svalsy = []               # Sagittal y

    #              Start of loop to make rays
    for i in range(-nrays,nrays + 1):
        r = i*dr                           # Radial poition
        rvals.append(r/entrance.maxRadius) # Record normalsied position
        #
        #         Make the m and s rays at test wavelength
        mray = ray.IntensityRay(entrance.point + Vector3d(0.0, r, 0.0), u, wave)
        sray = ray.IntensityRay(entrance.point + Vector3d(r, 0.0, 0.0), u, wave)
        #
        #         Add to pencil and propagate both back to clear lens
        pencil = ray.RayPencil(mray,sray).propagate(-entrance.maxRadius)
        #         propagate through lens to image surafce
        pencil *= lens
        pencil *= ip

        #            If rays valid (so not blocked), abberations to list
        if mray:
            mpos = mray.position - ref
            mvals.append(mpos.y)
        else:
            mvals.append(float("nan"))
        if sray:
            spos = sray.position - ref
            svalsx.append(spos.x)
            svalsy.append(spos.y)
        else:
            svalsx.append(float("nan"))
            svalsy.append(float("nan"))

    #     Return a list of three plots with suitable labels

    return [plt.plot(rvals,mvals, label="Meridional"),\
            plt.plot(rvals,svalsx,label="Sagittal x"),\
            plt.plot(rvals,svalsy,label="Sagittal y")]


def knifeEdgeTest(lens,angle = 0.0, knife = 0.0, wave = Default, design = Default, optimal = True, nrays = 50):
    """
    Function to give an image of the Focault Knife edge test of a lens as specifed angle.
    param lens, the lens under test
    param angle the angle of incident rays for the test
    param knife height of the knife from the PSF centre
    param wave the wavelength of the test
    param design the design wavelength (for the paraxial calculations)
    parar optimal, is test at optional PSF (else as paraxial)
    nrays number of rays in test
    return a OpticalImage of the test located at 2*flocal length from back nodal point
    """


    if isinstance(angle,float):
        u = Unit3d(Angle(angle))                     # direction of beam
    else:
        u = Unit3d(angle)

    cp = lens.cardinalPoints(design)             # Get the cardinal points
    fl = lens.focalLength(design)                # The focal length
    xsize = 3.0*lens.entranceAperture().maxRadius # Size of output feild
    
    output = OpticalImage(cp[3].propagate(2*fl,u),xsize,xsize) # Output image at 2fl from back nodal
    
    pencil = ray.RayPencil().addCollimatedBeam(lens,u,"array",nrays,wave)   # Make pencil
    pencil *= lens
    if optimal :
        psf = Psf().optimalArea(pencil,cp[1].z)          # Make optimal PSF
    else:
        psf = lens.imagePoint(u,design)                  # Use deign wavelength paraxial approx
    knifeEdge = KnifeEdgeAperture(psf,lens.exitAperture().maxRadius,knife,0.0)  # Make knife edge
    pencil *= knifeEdge                                  # propagate through knife edge
    pencil *= output                                     # Then to output (will give shadow image)

    return output
    
    
#
class WavePoint(Vector2d):
    """
    Call to hold a WavePoint being the optical 
    """
    def __init__(self,r, plane, refpt = None):
        """
        param ray the input ray
        param plane 
        param reference point
        """
        self.wavelength = r.wavelength
        self.set(r.pointInPlane(plane))                      # set self to point in plane
        distance = plane.getDistance(r.position,r.director)  # get distance from ray to plane

        if refpt != None:                # Deal with reference point
            ref = Vector3d(refpt.x, refpt.y, refpt.z - plane.getPoint().z)   # Ref point relative to plane
            
            xc = self.x - ref.x         # Position relative to ref
            yc = self.y - ref.y
            c = 1.0/ref.z               # Curvature
            u = r.director            # direction of ray
            #
            #     Distance to reference sphere (same code as in QuadricSurface)
            f = c*(xc*xc + yc*yc)
            g = u.z - c*(xc*u.x + yc*u.y)
            a = g*g - c*f
            if a < 0:
                raise ValueError("analysis.WavePoint: ray misses reference sphere")
            else:
                distance += f/(g + math.sqrt(a))
            
        # set total pathlength, pathleng of ray + pathlength to plane
        self.pathlength = r.pathlength + distance*r.refractiveindex.getValue(r.wavelength)
    #

    def __repr__(self):
        return "WavePoint: {0:s} wl : {1:7.4f} : pl : {2:8.5e}".format(str(self),self.wavelength,self.pathlength)

    #
    def getPhaselength(self):
        """
        Method to get the phase length, being 2*pi*pathelength/wavelength
        """
        return 2000.0*math.pi*self.pathlength/self.wavelength



#
class WavePointSet(list):
    """
    Class to hold a set of WavePoints
    """
    
    def __init__(self,pencil,plane,refpt = None):
        """
        Create a set of wavepoint from a pencil.
        """
        list.__init__(self)
        self.plane = plane               # Record plane
        self.maxRadius = 0.0             # MaxRadius
        if hasattr(plane, "maxRadius"):
            self.maxRadius = plane.maxRadius # Set to plane maxradius if defined
        for r in pencil:
            if r:
                self.add(WavePoint(r,plane,refpt))
                        
    def add(self,wp):
        """
        Metod to add a WavePoint,it append to WavePointSet and also set the currect
        maxradius of neded
        """
        rad = abs(wp)
        self.maxRadius = max(rad,self.maxRadius)
        self.append(wp)


    def zeroMean(self):
        """
        Method to zerom mean the data by substracting off the average pathlength
        from eack point.
        """
        pathsum = 0.0
        for w in self:
            pathsum += w.pathlength

        pathave = pathsum/len(self)
        for w in self:
            w.pathlength -= pathave

        return self

        
    def leastSqrError(self,zern):
        """
        Method to caculate the least sqr error between the datapoint and a two dimensional
        surface estmate
        """
        zern.radius = self.maxRadius
        
        sm = 0.0
        sumSqr = 0.0
        
        for w in self:
            z = zern.getValue(self)
            diff = w.pathlength - z
            sm += diff
            sumSqr += diff*diff

        n = len(slef)

        return sumSqr/(n*n) - pow(sm/n,2)


        


                    
        
        



    
