"""
Set of classes for analysis of optical systems
"""
import optics.ray as ray
from optics.psf import Psf
from optics.surface import OpticalPlane,ImagePlane,SurfaceInteraction,SphericalSurface,KnifeEdgeAperture
from optics.wavelength import Default,WavelengthColour
from vector import Vector2d,Vector3d,Unit3d,Angle
import matplotlib.pyplot as plt
import numpy as np
import math
import array

class TargetPlane(ImagePlane):
    """
    For a target plane, being at ImagePlane with target points or various types
    """
    
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

        col = WavelengthColour(self.wavelength)

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
    image is held in a nmpy array. If the first parameter is an ImagePlace then the reference point and x/y size will be automaticall y
    taken from the from this and the xsize / ysize parameters ignored.

    :param pt: reference point, or ImagePlane (Default = (0,0,0)
    :type pt: ImagePlane or Vector3d or float
    :param xpixel: xpixel size of image (default = 256) OR nmpy array of floats
    :type xpixel_or_im: numpy.array or int
    :param ypixel: ypixel size of image (Default - 256)
    :type ypixel: int
    :param xsize: x size of plane (Default = 200)
    :type xsize: float
    :param ysize: y size of plane (Default = 200)
    :type ysize: float
    

    """
    
    def __init__(self,pt =  Vector3d() ,xpixel = 256, ypixel = None, xsize = 200, ysize = None):
        """
        Form the OpticalImage with either blank array of nmpy image array
        """

        if isinstance(pt,ImagePlane):                          # Deal with ImagePlane
            ImagePlane.__init__(self,pt.point,pt.xsize,pt.ysize)
        else:
            if ysize == None:
                ysize = xsize
            ImagePlane.__init__(self,pt,xsize,ysize)             # Set underying ImagePlane
            
        if isinstance(xpixel,int):
            if ypixel == None:
                ypixel = xpixel
            self.image = np.zeros((xpixel,ypixel),dtype = float)  # Make array of zeros.
        else:
            self.image = xpixel                             # assume numpy array given
        self.xpixel,self.ypixel = self.image.shape          # set xpixel and ypixel from image data


    def __str__(self):
        """
        Implement str()
        """
        return "pt: {0:s} xpixel: {1:d} ypixel: {2:d} xsize: {3:6.4f} ysize: {4:6.4f}".\
            format(str(self.point),self.xpixel,self.ypixel,self.xsize,self.ysize)


    def getPixelSourcePoint(self,i,j):
        """
        Get pixel as i,j as a SourcePoint.

        :param i: the x pixel location
        :type i: int
        :param j: the y pixel location
        :type j: int
        :return: SourcePoint giving x,y,z and intensity of pixel in global coordinates.

        """
        x = self.xsize*(i/self.xpixel - 0.5)
        y = self.ysize*(j/self.ypixel - 0.5)
        
        return self.getSourcePoint(x,y,self.image[i,j])


    def getSurfaceInteraction(self,r):
        """
        Method to get back the surface interaction information for a ray and also add the ray to the image
        This also add the ray intensity to the cloeset pixel.

        :return: SurfaceInteraction.

        """
        
        #       get interaction with super class
        info = ImagePlane.getSurfaceInteraction(self,r)
        
        #            Add ray to pixel
        if not math.isnan(info.position.x) or not math.isnan(info.position.y) :
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

        :param ca: circular apereture (or lens) to fill
        :param i: the x the pixel coordinates.
        :type i: int
        :param j: the y pixel coordinate
        :type j: int
        :param nrays: number of ray across radius (default = 5)
        :type nray: int
        :param wave: wavelength of rays (default = Default)
        :type wave: float
        
        Note will return None of the pixel intensity is 0.0
        """
        source = self.getPixelSourcePoint(i,j)
        if source.intensity == 0.0:
            return None
        else:
            return ray.RayPencil().addSourceBeam(ca,source,"array",nrays,wave)

    def getImage(self, lens, ip, nrays = 5, wave = Default):
        """
        Method to get the image of OpticalPlane where the image localion is specifed
        by the supplied ImagePlane.
        
        :param lens: the lens system
        :param ip: ImagePlane
        :param nrays: number of rays on radius
        :param wave: wavelength of imaging (to do the actual tracing)
        :return: OpticalImage with same pixel resolution as the object

        """

        image = OpticalImage(ip,self.xpixel,self.ypixel)      # Form image

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

    def getSystemImage(self,lens,mag,nrays = 5, wave = Default, design = None):
        """
        Method to get the image of the object plane from and imaging system with specified lens and magnification.
        The location of the object and image planes are given by paraxial optics.

        :param lens: the imaging lens
        :type lens: OpticalGroup or Lens
        :param mag: The magnification between object and image (normally negative for imaginig system)
        :type mag: float
        :param nrays: number or rays across radius in simulation. (Default = 5)
        :type nrays: int
        :param wave: wavelength of rays in simulation (Default = optics.wavelength.Default)
        :type wave: float
        :param design: wavelength used for the paraxial location of the planes (Default = None) (same as wave)

        """
        if design == None:
            design = wave

        #     Get location of object and image planes and design wavelength
        obj,ima = lens.planePair(mag,self.xsize,self.ysize,design)
        self.setPoint(obj.point)        # Set self to correct location

        im = self.getImage(lens,ima,nrays,wave)     # get the image

        return im

        

    def addTestGrid(self, xgap = 10, ygap = None, intensity = 1.0 ):
        """
        Method to add a test grid being a grid of one pixel wide in a grid pattern.

        :param xgap: gap in x directions between pixels (defaults to 10)
        :type xgap: int
        :param ygap: gap in y directions between pixels, (defaults to xgap)
        :type ygap: int or None
        :param intensity: the intensity
        :type intensity: float

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
                    self.image[i,j] = intensity

        return self

    def draw(self):
        """
        Display the image via imshow with gray comlour map and correct etent
        """
        return plt.imshow(self.image,cmap=plt.cm.gray,\
                          extent=(-self.xsize/2+self.point.x,self.xsize/2+self.point.x,-self.ysize/2+self.point.y,self.ysize/2+self.point.y))
        



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

class AberrationPlot(object):
    """
    Class to form the standard transverse ray aberrations plots

    :param lens: the lens to be tested
    :param wave: the wavelength to apply the tests (Default = optics.wavelength.Default)
    :parm design: the wavelength of the design used in paraxial to give location of plane (Default = None, same as wave)
    :param nrays: number of rays across entrance aperture (Default = 50)

    The actual plot is implemented by the .draw() method
    """

    def __init__(self,lens,wave = Default, design = None,nrays = 50):
        """
        The conststructor.
        """
        self.lens = lens
        self.wavelength = float(wave)
        if design == None:
            self.design = self.wavelength
        else:
            self.design = float(design)

        self.nrays = nrays

    def draw(self,angle = 0.0, colour=["r","g","b"]):
        """
        Form and draw the plots at specified angle
        
        :param angle: the ray angle
        :type angle: float or Angle or Unit3d
        :param colour: line colours is three elemnts list, Default = ["r","g","b"])

        """
        if isinstance(angle,float):
            u = Unit3d(Angle(angle))                     # direction of beam
        else:
            u = Unit3d(angle)

        ref = self.lens.imagePoint(u,self.design)               # Get image point at design wavelength
        ip = self.lens.backFocalPlane(self.design)              # Make back focal plane to proagate to 

        ca = self.lens.entranceAperture()
        dr = ca.maxRadius/(self.nrays + 0.1)                   # Ray separation

        
        rvals = []                # Radius values
        mvals = []                # Meridional
        svalsx = []               # Sagittal x
        svalsy = []               # Sagittal y

        #              Start of loop to make rays
        for i in range(-self.nrays,self.nrays + 1):
            r = i*dr                           # Radial poition
            rvals.append(r/ca.maxRadius)       # Record normalsied position
            #
            #         Make the m and s rays at test wavelength
            mray = ray.IntensityRay(ca.point + Vector3d(0.0, r, 0.0), u, self.wavelength)
            sray = ray.IntensityRay(ca.point + Vector3d(r, 0.0, 0.0), u, self.wavelength)
        
            #       Add to pencil and propagate both back to clear lens
            pencil = ray.RayPencil(mray,sray).propagate(-ca.maxRadius)
            #         propagate through lens to image surafce
            pencil *= self.lens
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
            

        # plots with suitable labels to the current axis.

        plt.plot(rvals,mvals, color = colour[0], label="Meridional")
        plt.plot(rvals,svalsx,color = colour[1], label="Sagittal x")
        plt.plot(rvals,svalsy,color = colour[2], label="Sagittal y")
        


class KnifeEdgeTest(object):
    """
    Class to implement a knife edeg test with methods to deconfigure the knife
    
    :param lens: the lens under test
    :param angle: the angle of the analysis
    :param wave: the test wavelength
    :param design: the design wavelength.

    """
    def __init__(self,lens,angle, wave = Default,design = None):
        """
        The constrcutor
        """
        self.lens = lens
        if isinstance(angle,float):
            self.u = Unit3d(Angle(angle))
        else:
            self.u = Unit3d(angle)
        self.wavelength = float(wave)
        if design == None:
            self.design = self.wavelength
        else:
            self.design = float(design)

        # Set up a basic knife edge aperture at 0,0,0 with default distance and angle.
        self.knife = KnifeEdgeAperture(0.0,self.lens.exitAperture().maxRadius)

    def setKnife(self,knife = 0.0, angle = 0.0):
        """
        Set or reset knife distance and angle.
        
        :param knife: distance from optical axis (Default = 0.0)
        :param angle: angle of knife (Default = 0.0

        """
        self.knife.setKnife(knife,angle)

    def getImage(self,optimal,xpixel = 256, ypixel = None,nrays = 50):
        """
        Get the knife edge image
        """
        cp = self.lens.cardinalPoints(self.design)    # Get the cardinal points
        fl = self.lens.backFocalLength(self.design)
        xsize = 3.0*self.lens.entranceAperture().maxRadius    # Size of output feild
    
       
    
        pencil = ray.RayPencil().addCollimatedBeam(self.lens,self.u,"array",nrays,self.wavelength)   # Make pencil
        pencil *= self.lens                                                                 # Propagate through lens.
        if optimal :
            psf = Psf().optimalArea(pencil,self.lens.backFocalPlane(self.design))            # Make optimal PSF
        else:
            psf = self.lens.imagePoint(self.u,self.design)       # Use deign wavelength paraxial approx
            
        self.knife.setPoint(psf)                              # Set the position of the knife
        output = OpticalImage(psf.propagate(fl,self.u),xpixel,ypixel,xsize,xsize) 
        pencil *= self.knife                                  # propagate through knife edge
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


        


                    
        
        



    
