"""
Set of classes for analysis of optical systems
"""
import optics.ray as ray
from optics.psf import Psf
from optics.surface import OpticalPlane,ImagePlane,SurfaceInteraction,SphericalSurface,KnifeAperture,CircularAperture
from optics.wavelength import Default,WavelengthColour,AirIndex
from vector import Vector2d,Vector3d,Unit3d,Angle
import matplotlib.pyplot as plt
import numpy as np
import math

class TargetPlane(ImagePlane):
    """
    For a target plane, being at ImagePlane with target points or various types. 
    Targets are held as Vector2d in the local plane coordinates.

    :param pt: the reference point for the target plane (Default = 0.0). 
    :type pt: float of Vector3d
    :param xsize: x size of plane (Default = 100mm)
    :type xsize: float
    :param ysize: ysize or target plane (Default = xsize )
    :type ysize: float
    :param wave: wavelength of targets (Default = optics.wavelength.Default)
    :type wave: float

    The inital TragetPlane is empty. use .add() or .addGrid() to add targets.
    """
    
    def __init__(self,pt = 0.0 ,xsize = 100.00, ysize = None,wave = Default):
        """
        Constuctor with
        """
        if isinstance(pt,ImagePlane):
            ImagePlane.__init__(self,pt,pt.xsize,pt.ysize)
        elif isinstance(pt,CircularAperture):
            ImagePlane.__init__(self,pt,2*pt.getRadius())
        else:
            ImagePlane.__init__(self,pt,xsize,ysize)   # Initialse underlying ImagePlane
        self.wavelength = wave
        self.targets = []                          # List of targets to be added



    def __str__(self):
        """
        Update str
        """
        return ImagePlane.__str__(self) + " targets : {0:d}".format(len(self.targets))
        
    def add(self,target,y = None):
        """
        Add a target, or list of targets. 

        :param target: target or list of to be added. 
        :param y: y component if target is x,y pair
        """

        
        if isinstance(target,ray.RayPencil):    # Deal with ray pencil
            for r in target:
                if r:
                    self.add(r)

        elif isinstance(target,list):
            for t in target:
                self.add(t)
        elif isinstance(target,ray.IntensityRay): # Deal with ray
            self.wavelength = target.wavelength
            v = target.pointInPlane(self)
            self.targets.append(v)
        elif isinstance(target,Vector2d):
            self.targets.append(Vector2d(target))
        elif isinstance(target,float):          # Assume x,y given
            self.tragets.append(Vector2d(target,y))
                     
        else:
            raise TypeError("analysis.TargetPlane.illegal type")

        return self

    #
    def addGrid(self,xn,yn = 0,radius = float("inf")):
        """
        Fill the TargetPlace with a set regular set of targets

        :param xn: number of targets across horizontal
        :param yn: numbers of tarets across if 0 or negative, n will  be set to that the targets are on a square grid
        :param radius: float, have targets only in a masked of  specified radius.

        Note: number of targets will be rounded to ensure an
        odd number across array so there will always be one
        at (0,0)
        """
        
        dx = self.xsize/(xn - 1 + 0.1)
        #        dx = 0.5*self.xsize/(xn/2 + 0.1)
        if yn > 0 :
            dy = self.ysize/(yn - 1 + 0.1)
            # dy = 0.5*self.ysize/(yn/2 + 0.1)
        else:
            dy = dx
            yn = int(round(self.ysize/dy))

        for j in range(-yn//2,yn//2+1):
            for i in range(-xn//2,xn//2+1):
                y = dy*j
                x = dx*i
                if x*x + y*y < radius*radius and \
                   abs(x) < self.xsize/2 and abs(y) < self.ysize/2:
                     self.add(Vector2d(x,y))
        return self

    def rayPencil(self,pt_or_u,wave = Default, intensity = 1.0):
        """
        Get an intensity RayPenci, one ray from each target
        
        :param pt_or_u:  Vector3d or Unit3d, of Position each ray will pass through this point or direction of rays, this        
        :param wave: wavelength, (defaults of Default)
        :param intensity: intensity, (defaults to 1.0)
        """
        pencil = ray.RayPencil()
        for t in self.targets:
            #                Start position of ray
            pos = self.getSourcePoint(t)
            if isinstance(pt_or_u,Unit3d):
                u = Unit3d(pt_or_u)
            else:
                u = Unit3d(pt_or_u - pos)
            r = ray.IntensityRay(pos,u,wave,intensity)
            pencil.append(r)
        return pencil

    
    def getPencils(self,ca,key = "array", nrays = 10, index = AirIndex()):
        """
        Method to get RayPencils from each target in trem in an itterator

        :param ca: Circular aperture to be filled
        :param key: the pencil key, (Default = "array")
        :param nrays: numner of arrays across radius (Default = 10)
        :param index: Starting index, (Default = AirIndex())

        The wavelength is set by self.wavelength
        """
        for t in self.targets:
            if t:                              # Check target os valid
                pt = self.getSourcePoint(t)    # Tarhget in glabal
                pencil = ray.RayPencil().addSourceBeam(ca,pt,key,nrays,self.wavelength,index)
                yield pencil
        

        
    def draw(self):
        """
        Draw the plane to the current axis.
        """
        pt = self.getPoint()

        #      Make a frame round the plane and plot it in black
        xframe = [pt.x - self.xsize/2, pt.x + self.xsize/2,\
                  pt.x + self.xsize/2, pt.x - self.xsize/2,\
                  pt.x - self.xsize/2]
        yframe = [pt.y + self.ysize/2, pt.y + self.ysize/2,\
                  pt.y - self.ysize/2, pt.y - self.ysize/2,
                  pt.y + self.ysize/2]
        plt.plot(xframe,yframe,"k")
        #       Now plot the targets at "x"
        xpt = []
        ypt = []
        
        for t in self.targets:
            pt = self.getSourcePoint(t)
            xpt.append(pt.x)
            ypt.append(pt.y)

        col = WavelengthColour(self.wavelength)
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


class KnifeTest(object):
    """
    Class to implement a knife edeg test with methods to deconfigure the knife
    
    :param lens: the lens under test
    :param angle: the angle of the analysis
    :param wave: the test wavelength
    :param design: the design wavelength.

    """
    def __init__(self,lens,angle,wave=Default,design=None):
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

        # Set up a basic knife edge aperture at 0,0,0 with default knife distance, angle and shift
        self.knife = KnifeAperture(0.0,self.lens.exitAperture().maxRadius)

    def setKnife(self,knife = 0.0, angle = 0.0, shift = 0.0):
        """
        Set or reset knife distance, angle and shift
        
        :param knife: distance from optical axis (Default = 0.0)
        :type knife: float
        :param angle: angle of knife (Default = 0.0)
        :type angle: float
        :param shift: shift along the zaxis relative to reference point
        :type shift: float

        """
        self.knife.setKnife(knife,angle,shift)
        return self


    def setWire(self,wire,thickness = 0.01):
        """
        Set into wire mode
        """
        self.knife.setWire(wire,thickness)
        return self

    def getImage(self,refopt = 1,xpixel = 256, ypixel = None,nrays = 50):
        """
        Get the knife edge image
        """
        cp = self.lens.cardinalPoints(self.design)    # Get the cardinal points
        fl = self.lens.backFocalLength(self.design)
        xsize = 3.0*self.lens.entranceAperture().maxRadius    # Size of output feild
    
       
    
        pencil = ray.RayPencil().addCollimatedBeam(self.lens,self.u,"array",nrays,self.wavelength)   # Make pencil
        pencil *= self.lens    # Propagate through lens.

        if refopt == 0:
            psf = self.lens.imagePoint(self.u,self.design)       # Use deign wavelength paraxial approx
        elif refopt == 1:
            psf = Psf().setWithRays(pencil,self.lens.backFocalPlane(self.design))
        elif refopt == 2:
            psf = Psf().optimalArea(pencil,self.lens.backFocalPlane(self.design))            # Make optimal PSF
        else:
            print("Illegal refopt")
            
            
        self.knife.setPoint(psf)                              # Set the position of the knife
        output = OpticalImage(psf.propagate(fl,self.u),xpixel,ypixel,xsize,xsize) 
        pencil *= self.knife                                  # propagate through knife edge
        pencil *= output                                     # Then to output (will give shadow image)

        return output
    



    
