"""
Set of classes to handle wavefronts, analysis and interferoneters.
"""
from vector import Vector2d,Vector3d,Unit3d,Angle
import math
import cmath
from optics.surface import OpticalPlane
from optics.psf import Psf
from optics.wavelength import Default
from optics.zernike import opticalZernike,opticalZernikeName
import optics.ray as ray
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import tio



class WavePoint(Vector2d):
    """
    Class to hold a WavePoint being the optical path length or phase at a point in a plane.
    """

    def __init__(self,pt = Vector2d() ,pathlength = 0.0 , wavelength = Default):
        """
        Basic constructor for a wavepoint

        :param pt: point in the plane
        :type pt: Vector2d
        :param pathlength : the pathlength (default = 0.0)
        :type pathlength : float 
        :param wavelength : the wavelength (default = wavelength.Default)
        :type wavelength :float
        """
        self.set(pt)
        self.pathlength = pathlength
        self.wavelength = wavelength

    def __str__(self):
        return "{0:s} wl : {1:7.4f} : pl : {2:8.5e}".format(Vector2d.__str__(self),self.wavelength,self.pathlength)



    def getPathLength(self):
        """
        Method to get the pathlength in mm
        """
        return  self.pathlength
    
    def getPhaseLength(self):
        """
        Method to get the phase length, being 2*pi*pathelength/wavelength. 
        """
        return 2000.0*math.pi*self.pathlength/self.wavelength

    def setWithWaveFront(self,wf,pt = None):
        """
        Set a single WavePoint with a wavefront. 
        If pt is None, then the current value position of the wavepoint is used.

        :param wf: The WaveFront
        :type wf: WaveFront
        :param pt: the point, of None then the current value of self is iused

        """
        if pt != None:
            self.set(pt)
        p = wf.getValue(self)                                  # Get the phase
        self.pathlength = p*self.wavelength/(2000.0*math.pi)   # convert to distance in mm  
        return self

    
    def setWithRay(self,ray, plane, refpt = None):
        """
        Set the value in a specifed plane using a ray and optional reference point. This is
        the main method used to cancluate wavefront abberations in a plane wrt to a specificied
        refererence point.
        
        :param ray: The intensity ray
        :param plane: The flat plan in which to form the wavepoint 
        :param refpt: The reference point. This assumed to be the image point.

        """

        if isinstance(plane,float):                               # If plane as float, make a plane
            plane = OpticalPlane(plane)
        
        self.wavelength = ray.wavelength
        self.set(ray.pointInPlane(plane))                          # set self to point in plane.
        distance = plane.getDistance(ray.position,ray.director)    # get distance from ray to plane.

        if refpt != None:                                          # Deal with reference point
            ref = refpt - plane.getPoint()                         # Ref point relative to centre of plane.
            
            xc = self.x - ref.x         # Position relative to ref
            yc = self.y - ref.y
            c = 1.0/ref.z               # Curvature
            u = ray.director             # direction of ray
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
        self.pathlength = ray.pathlength + distance*ray.refractiveindex.getValue(ray.wavelength)

        return self



class WavePointSet(list):
    """
    Class to hold a list of WavePoints, consists of a list of WavePoints with manipulation methods.
    """

    def __init__(self,radius = 0.0, *args):
        """
        Set max radius and append any WavPoints given.

        :param radius: Radius distribution in plane, (default = 0.0) 
        This is auto-expanded as WavePoints are added if needed.
        :type radius: float
        :param \*args: WavePoint to be added on creation.

        """
        list.__init__(self)
        self.maxRadius = radius
        self.wavelength = Default
        for wp in args:
            self.add(wp)

    
    def setWithRays(self,pencil,plane,refpt = None):
        """
        Create a set of wavepoint from a pencil, in a specified plane with an optional reference point
        
        :param pencil: The raypencil containing a list of arrays
        :param plane: the plane of the wavepoints.
        :param refpt: The reference point, may be None.

        """
       
        if isinstance(plane,float):      # If plane as float, make a plane
            plane = OpticalPlane(plane)
        self.plane = plane               # Record plane
        if hasattr(plane, "maxRadius"):
            self.maxRadius = plane.maxRadius # Set to plane maxradius if defined
        for r in pencil:
            if r:                        # Only take valid rays.
                self.add(WavePoint().setWithRay(r,plane,refpt))

        return self
                        
    def add(self,wp):
        """
        Method to add a WavePoint,it append to WavePointSet and also extend the currect
        maxradius if needed.

        :param wp: the WavePoint to be added
        :type wp: WavePoint

        """
        rad = abs(wp)
        self.maxRadius = max(rad,self.maxRadius)
        self.append(wp)
        self.wavelength = wp.wavelength


    def zeroMean(self):
        """
        Method to zero-mean the data by substracting off the average pathlength from each point.

        :return: self
        """
        pathsum = 0.0
        for w in self:
            pathsum += w.pathlength

        pathave = pathsum/len(self)
        for w in self:
            w.pathlength -= pathave

        return self



    def setWithWaveFront(self,wf):
        """
        Set the current set of wavepoints with a WaveFront
        """
        zern.radius = self.maxRadius
        for w in self:
            w.setWithWaveFront(zern)

        return self

    def getPhaseValues(self):
        """
        Get the Phase length Values as a numpy array. This is an internal method used 
        in fitting and not nornmally called bu the user.

        :return: np.array of floats 
        """
        y = []
        for w in self:
            y.append(w.getPhaseLength())
        return np.array(y)

    def fitZernikeFunctionNine(self,x,a,b,c,d,e,f,g,h,i):
        """
        Fitting function for a 9 parameter (4th order) Zernike expansion
        in a form that can be used by the SciPy curve\_fit 

        :return: np.array of values at the WavePoints.
        """
        ze = ZernikeWaveFront(self.maxRadius,self.wavelength,a,b,c,d,e,f,g,h,i)
        y = []
        for w in self:
            y.append(ze.getValue(w))
        return np.array(y)
    
    def fitZernikeFunctionSixteen(self,x,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p):
        """
        Fitting function for a 16 parameter (6th order) Zernike expansion
        in a form that can be used by the SciPy curve\_fit 

        :return: np.array of values at the WavePoints.
        """
        ze = ZernikeWaveFront(self.maxRadius,self.wavelength,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)
        y = []
        for w in self:
            y.append(ze.getValue(w))
        return np.array(y)

    def fitZernikeFunctionTwentyFive(self,xv,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,\
                                     q,r,s,t,u,v,w,x,y):
        """
        Fitting function for a 25 parameter (8th order) Zernike expansion
        in a form that can be used by the SciPy curve\_fit 

        :return: np.array of values at the WavePoints.
        """
        
        ze = ZernikeWaveFront(self.maxRadius,self.wavelength,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,\
                              q,r,s,t,u,v,w,x,y)
        y = []
        for w in self:
            y.append(ze.getValue(w))
        return np.array(y)


    def fitZernike(self,order = 4):
        """
        Fit a Zernike Expansion to the WaveFront to specified order.
        This used the SciPi curve\_fit methiod to do the actual fitting work. 
        Note it will Zero Mean first.
        
        :param order: Order of expansion, 4,6 and 8 implemented (Default = 4)
        :type order: int
        :return: the fitted ZernikeWaveFront of unit radius.

        Also self.zerr is avaialble as a np.array of the errors terms from the
        fit.
        """
        self.zeroMean()
        #                    Get the phase values as a np array
        y = self.getPhaseValues()

        #      Do a fix with using curve_fit
        if order == 4:
            popt,pcov = curve_fit(self.fitZernikeFunctionNine,self,y)
        elif order == 6:
            popt,pcov = curve_fit(self.fitZernikeFunctionSixteen,self,y)
        elif order == 8:
            popt,pcov = curve_fit(self.fitZernikeFunctionTwentyFive,self,y)
        else:
            print("Order error in fitZernike")
            return None
            
        # record the errors to they can be read if required 
        self.zerr = np.sqrt(np.diag(pcov))

        #      Return the result as a Zernike Wavefront with unit radius.
        ze = ZernikeWaveFront(1.0,self.wavelength,*popt)
        return ze



    def fitSeidelFunction(self,xv,a,b,c,d,e,f):
        """
        Fit function for Seidel aberrations with 6 parameters
        """
        se = SeidelWaveFront(self.maxRadius,self.wavelength,self.fieldangle,a,b,c,d,e,f)
        y = []
        for w in self:
            y.append(se.getValue(w))
        return np.array(y)

    
    def fitSeidel(self,angle = 0.0):
        """
        Method to fit the seidel abberations 
        """
        self.fieldangle = float(angle)
        
        self.zeroMean()           # Zero mean the data
        y = self.getPhaseValues() # Get phase values as np array
        #                         # Do the fit
        popt,pcov = curve_fit(self.fitSeidelFunction,self,y)
        self.serr = np.sqrt(np.diag(pcov))    

        se = SeidelWaveFront(1.0,self.wavelength,self.fieldangle,*popt)
        return se
        
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

        n = len(self)

        return sumSqr/(n*n) - pow(sm/n,2)


class WaveFrontAnalysis(object):
    """
    Class to implement a set of wavefront analysis of a specified lens.

    :param lens: The lens to be tested
    :type lens: optics.lens.Lens
    :param design: the design wavelenth (Default = 0.55)
    :type design: float

    """

    def __init__(self,lens, design = Default):
        """
        Constructor
        """
        self.lens = lens
        self.design = design
        self.refpt = Vector3d()
        self.ip = self.lens.backFocalPlane(self.design)              # Make back focal plane to proagate to 

    def getWavePointSet(self,u,wave = Default,nrays = 10,refopt = 1):
        """
        Get the wavepointset for collimated beam in the exitpupil of the lens

        :param u: the angle of input beam
        :type u: Unit3d or float
        :param wave: analysis wavelength
        :type wave: float
        :param nrays: number of raays across input aperture, (Default = 10)
        :type nrays: int
        :param refopt: reference point option, 0 = paraxial, 1 = centre of PSF in image plane, 2 = optimal area PSF
        :type refopt: int
        """
        #       Sort out angle
        if isinstance(u,float) or isinstance(u,int):
            u = Unit3d(Angle(u))
        else:
            u = Unit3d(u)
        
        #      Make pencil with array and record path
        pen = ray.RayPencil().addCollimatedBeam(self.lens,u,"array",nrays=nrays,wave=wave,path=True)
        ep = self.lens.exitPupil(self.design)         # Exit pupil of lens

        pen *= self.lens                              # Penil through lens

        if refopt  == 0:
            self.refpt = self.lens.imagePoint(u,self.design)    # Paraxial point location
        elif refopt == 1:
            self.refpt = Psf().setWithRays(pen,self.ip)              # Centre of PSF in image plane
        elif refopt == 2:
            self.refpt = Psf().optimalArea(pen,self.ip)              # Optimal area PSF, not in image plane
        else:
            print("Illegal ref type")

        #     Form the wavefront
        wf = WavePointSet().setWithRays(pen,ep,self.refpt)
        return wf
    


    def fitZernike(self,u,wave = Default, order = 4, refopt = 1):
        """
        Fit zernike to wavefront in exit pupil with sensible defaults

        :param u: angle of collimated beam
        :type u: Unit3d or float
        :param wave: analysis wavelength
        :type wave: float
        :param order: order of Zernike, must by 4,6,8 only
        :type order: int
        :param refopt: option for reference point
        :type refopt: int

        """
        nrays = 3*order

        wf = self.getWavePointSet(u,wave,nrays,refopt)
        ze = wf.fitZernike(order)

        return ze

    def fitSeidel(self,angle,wave=Default,refopt = 1):
        """
        Fit the Seidel aberrations 
        """
        nrays = 10
        wf = self.getWavePointSet(angle,wave,nrays,refopt)
        se = wf.fitSeidel(angle)
        return se
        
    def drawAberrationPlot(self,angle,wave = Default,colour=["r","g","b"],legend = "lower left"):
        """
        Form and draw the plots at specified angle
        
        :param angle: the ray angle
        :type angle: float or Angle or Unit3d
        :param colour: line colours is three elemnts list, Default = ["r","g","b"])
        :param legend: location of ledgend, (Default = "lower left")
        :type legend: str
        """
        if isinstance(angle,float):
            u = Unit3d(Angle(angle))                     # direction of beam
        else:
            u = Unit3d(angle)
            angle = Angle(u).theta

        self.ref = self.lens.imagePoint(u,self.design)    # Get image point at design wavelength
        

        nrays = 50
        ca = self.lens.entranceAperture()
        dr = ca.maxRadius/(nrays + 0.1)                   # Ray separation

        
        rvals = []                # Radius values
        mvals = []                # Meridional
        svalsx = []               # Sagittal x
        svalsy = []               # Sagittal y

        #              Start of loop to make rays
        for i in range(-nrays,nrays + 1):
            r = i*dr                           # Radial poition
            rvals.append(r/ca.maxRadius)       # Record normalsied position
            #
            #         Make the m and s rays at test wavelength
            mray = ray.IntensityRay(ca.point + Vector3d(0.0, r, 0.0), u, wave)
            sray = ray.IntensityRay(ca.point + Vector3d(r, 0.0, 0.0), u, wave)
        
            #       Add to pencil and propagate both back to clear lens
            pencil = ray.RayPencil(mray,sray).propagate(-ca.maxRadius)
            #         propagate through lens to image surafce
            pencil *= self.lens
            pencil *= self.ip

            #            If rays valid (so not blocked), abberations to list
            if mray:
                mpos = mray.position - self.ref
                mvals.append(mpos.y)
            else:
                mvals.append(float("nan"))
            if sray:
                spos = sray.position - self.ref
                svalsx.append(spos.x)
                svalsy.append(spos.y)
            else:
                svalsx.append(float("nan"))
                svalsy.append(float("nan"))
            

        # plots with suitable labels to the current axis.

        plt.plot(rvals,mvals, color = colour[0], label="Meridional")
        plt.plot(rvals,svalsx,color = colour[1], label="Sagittal x")
        plt.plot(rvals,svalsy,color = colour[2], label="Sagittal y")
        plt.title("{0:s}: a: {1:4.2f} w: {2:4.2f} d: {3:4.2f}".\
                  format(self.lens.title,math.degrees(angle),wave,\
                  self.design))
        plt.xlim(-1.0,1.0)
        plt.legend(loc=legend,fontsize="small")
        plt.grid()




class WaveFront(object):
    """
    Basic wavefront class, abstract at this point
    """
    def __init__(self,radius = 1.0, wavelength = Default):
        """
        Constuctor to form WaveFront object, only wavelength and maxradius is set here.
        """
        self.wavelength = float(wavelength)
        self.radius = float(radius)

    def getValue(self,x,y = None):
        """
        Get the valu at secified point

        :param x: the x value of Vector2d
        :type x: float or Vector2d
        :param y: the y value or None
        :type y: float of None
        :return: the value, initically set to "NaN"
        """
        if isinstance(x,Vector2d):
            y = x.y
            x = x.x
            
        x /= self.radius              # Normalise
        y /= self.radius

        return self._getValue(x,y)

    def _getValue(self,x,y):
        """
        Abstarct getValue, needs to be defined for extending classes
        """
        print("_getValue called by accideent")
        return float("nan")

    def getImage(self,size = 256):
        """
        Get the phase image of the expansion as an array of np.array floats of specified size.

        Note: the pixel elements outside the unit circle will be set to NaN.

        :param size: size of image (Default = 256)
        :type size: int
        :return: two dimensional np.ndarray of type float.

        """
        im = np.empty((size,size),dtype = float)      # Make Empty array
            
        xmax,ymax = im.shape
        ycentre = ymax/2.0
        xcentre = xmax/2.0

        for (i,j),v in np.ndenumerate(im):
            y = (j - ycentre)*self.radius/ycentre       # In range -1.0 to 1.0
            x = (i - xcentre)*self.radius/xcentre       # In range -1.0 to 1.0
            im[i,j] = self.getValue(x,y)                # Will be NaN for outside circle

        return im



    def getPSF(self,size = 256, log = True):
        """
        Get the PSF by fourier means.

        :param size: size of PSF array, Default = 256
        :type size: int
        :param log: logical if log of intensity taken, Default = True
        :type log: bool
        :return: two dimensional np.array 
        """
        im = self.getImage(size)
        r = np.cos(im)     # Real part
        i = np.sin(im)     # Imaginary part
        z = r + 1j*i       # Combine to form complex array

        z = np.nan_to_num(z) # Set NaN to Zero 
        psf = np.fft.fft2(z)
        psf = np.fft.fftshift(psf)  # Shift to centre
        psf = abs(psf)
        if log :                # Take the log if required.
            psf = np.log(psf + 1.0)

        return psf


    def plotPSF(self,size = 256, log = True, key = "b"):
        """
        Make plot of the PSF along the x and y axis. Plot to current defaults plt.plot()
      
        :param size: size of plot (Default = 256)
        :type size: int
        :param log: is log taken before plot, Default = True
        :type log: bool
        :param key: what is plotted, can be "h", "v" or "b" 
        :type key: str


        """
        psf = self.getPSF(size,log)
        xmax,ymax = psf.shape

        pData = np.linspace(-1.0,1.0,size)
        if key.startswith("b") or key.startswith("h"):      # Do horizontal plot 
            plt.plot(pData,psf[ymax//2],label = "Horizontal")
        if key.startswith("b") or key.startswith("v"):      # Do vertical plot
            vd = np.zeros(ymax)
            ic = xmax//2
            for j in range(0,ymax):
                vd[j] = psf[j,ic]
            plt.plot(pData,vd,label = "Vertical")
            
        plt.xlim(-1.0,1.0)
        plt.grid()
        plt.xlabel("Normalised position")
        plt.ylabel("PSF")
        plt.title("Plot of PSF")
        plt.legend(loc="upper right",fontsize="small")
        


    def getOTF(self,size = 128, key = "h"):
        """
        Get the one-dimenensioal normalsied OFT as np array in either horizontal or vertical diection.

        Note: this is calcualted in real space and can be slow for large size (128 gives sensible results).

        :param size: the number of point in the OTF (Default = 128)
        :type size: int
        :param key: horizontal or vertical "h" or "v"
        :type key: str

        """
        im = self.getImage(size)     # Get the phase image
        # Make the complex image and it complex conjugate
        r = np.cos(im)
        i = np.sin(im)
        z = r + 1j*i
        zc = np.conj(z)

        xsize,ysize = im.shape

        horizontal = key.startswith("h")    # Set logicval
        
        if horizontal:              # Sort out ditection of shift
            shiftSize = xsize
            fullRange = range(0,ysize)
        else:
            shiftSize = ysize
            fullRange = range(0,xsize)
        

        otfData = np.zeros(shiftSize)  # np array to hold the OFT

        #      Loop over the shifts
        for shift in range(0,shiftSize): 
            otf = 0.0
            shiftRange = range(shift,shiftSize)

            #
            #       Do 2d sum over overlap area
            for i in shiftRange:
                for j in fullRange:
                    ish = i - shift
                    if horizontal:
                        zr = z[i,j]
                        zl = zc[ish,j]
                    else:
                        zr = z[j,i]
                        zl = zc[j,ish]
                    if not (cmath.isnan(zr) or cmath.isnan(zl)) :
                        otf += (zr * zl).real
            otfData[shift] = otf

        #        Normalise this output
        max = np.amax(otfData)
        otfData /= max
        return otfData    
        



    def plotOTF(self,size = 128, key = "h" , ideal = True ):
        """
        Calcualte and plot OFT with sensible plot paramters. 
        It will plot hotizontal / vertical / both normalsied OTF with optional ideal plot

        :param size: number of point in shift, (Default = 128). 
        :type size: int
        :param key: Plot key, may be "h", "v" or "b" for horizontal / vertical / both
        :type key: str
        :param ideal: plot ideal OTF for reference (Default = True)
        :type ideal: bool

        """
        shiftData = np.linspace(0.0,1.0,size)

        if key.startswith("h") or key.startswith("b") :
            otfData = self.getOTF(size,key = "h")   # Get the OTF
            plt.plot(shiftData,otfData,label="Horizontal")
        if key.startswith("v") or key.startswith("b") :
            otfData = self.getOTF(size,key = "v")   # Get the OTF
            plt.plot(shiftData,otfData,label="Vertical")

        #
        #       Add ideal for reference
        if ideal:
            idealFn = lambda x: 2.0/math.pi*(np.arccos(x) - x*np.sqrt(1 - x**2))
            plt.plot(shiftData,idealFn(shiftData),"k--",label="Ideal")
            
        plt.xlim(0.0,1.0)
        plt.grid()
        plt.xlabel("Normalised spatial frequency")
        plt.ylabel("OFT")
        plt.title("Plot of OTF")
        plt.legend(loc="upper right",fontsize="small")

        

    def plotImage(self,size = 256):
        """
        Plot data is a np.array in extent +/- 1.0

        :param size: the size of the image in pixel, (Default = 256)
        :type size: int

        """
        im = self.getImage(size)
        plt.imshow(im,cmap=plt.cm.gray,extent=(-1.0,1.0,-1.0,1.0))

    

    def fromFile(self,fn = None):
        """
        Read a wavefront from a file
        """
        if fn == None:
            wfile = tio.openFile("Wavefront file","r","wf")
        else:
            fn = tio.getExpandedFilename(fn)   # Sort out logicals
            if not fn.endswith("wf"):        # Append ".wf" if not given
                fn += ".wf"
            wfile= open(fn,"r")             # open file

        #          read file and process one line at a time
        #
        coef = []                          # Local coefficients
        wave = self.wavelength
        rad = self.radius
        type ="zzz"
        theta = 0.0
        for line in wfile.readlines():
            line = line.strip()
            if not line.startswith("#") and len(line) > 0:   # Kill comments and blanks
                token = line.split()

                if token[0].startswith("type"):
                    type = str(token[1])
                elif token[0].startswith("radius"):
                    rad = float(token[1])
                elif token[0].startswith("wave"):
                    wave = float(token[1])
                elif token[0].startswith("field"):   # Only for Seidel
                    theta = float(token[1])      
                else: # Assume a coefficient
                    c = float(token[0])
                    coef.append(c)

        #          Work cout what we have
        if type.startswith("seid"):
            return SeidelWaveFront(coef,theta,rad,wave)
        elif type.startswith("zernike"):
            return ZernikeWaveFront(rad,wave,coef)
        elif type.startswith("poly"):
            return PolynomialWaveFront(rad,wave,coef)
        elif type.startswith("king"):
            return KingslakeWaveFront(rad,wave,coef)
        else:
            print("WaveFront.readFromFile: unknown type : " + str(type))
            return None


class KingslakeWaveFront(WaveFront):
    """
    Class to implement a simple Kingslake with 6 coefficeints supplied as a list. Coefficents are is wavelengths
    to match the functions in Malacara.

    :param radius: the maximum radius (Default = 1.0)
    :param wavelength: the wavelnegth (Default = w.Default)
    :param \*args: the 6 coefficents as parameters of in a list.
    """

    def __init__(self, radius = 1.0, wavelength = Default, *args):

        
        WaveFront.__init__(self,radius,wavelength)
        self.coef = []           # List of coefficients
        for z in args:
            if isinstance(z,list):
                self.coef.extend(z)
            elif isinstance(z,float):
                self.coef.append(z)

    def __str__(self):
        " Five the string"
        s = "A: {0:5.3f} B: {1:5.3f} C: {2:5.3f} D: {3:5.3f} E: {4:5.3f} F: {5:5.3f}". \
            format(self.coef[0],self.coef[1],self.coef[2],self.coef[3],self.coef[4],self.coef[5])
        return s
    

    def _getValue(self,x,y):
        """
        Internal method to get the value at x/y, assume to be normalsied
        """

        rSqr = x*x + y*y
        if rSqr > 1.0:
            return float("nan")
        
        v = self.coef[0]*rSqr*rSqr + self.coef[1]*y*rSqr + self.coef[2]*(x*x + 3.0*y*y) + self.coef[3]*rSqr + \
            self.coef[4]*x + self.coef[5]*y
        return 2.0*math.pi*v

        
seidelNames = ("Defocus","Spherical Aberration","Coma","Astigmatism","Field Curvature","Distortion")

class SeidelWaveFront(WaveFront):
    """
    Class to hold the Sidel abberations.
    
    :param radius: Radius of plane, (Default = 1.0)
    :type radius: float
    :param theta: field angle
    :type theta: float
    :param coef: 
    """

    def __init__(self,radius = 1.0, wavelength = Default, theta = 0.0, *args):
        """
        Form the siedel class with the coefficients, note the coefficients are in microns.
        """
        WaveFront.__init__(self,radius,wavelength)
        self.theta = float(theta)
        self.coef = []
        for z in args:
            if isinstance(z,list):
                self.coef.extend(z)
            elif isinstance(z,float):
                self.coef.append(z)
                

    def __str__(self):
        """
        The str to print out values on single line.
        """
        s = "S0: {0:5.3f} S1: {1:5.3f} S2: {2:5.3f} S3: {3:5.3f} S4 : {4:5.3f} S5: {5:5.3f} Theta: {6:5.3f}". \
            format(self.coef[0],self.coef[1],self.coef[2],self.coef[3],self.coef[4],self.coef[5],self.theta)
        return s

    def __repr__(self):
        """   Detailed of the Seidel in tabular form.
        """
        s = "Seidel Aberrations \n"
        for i in range(0,len(self.coef)):
            s += "{0:<24}: {1:7.5f}\n".format(seidelNames[i],self.coef[i])
        s += "{0:<24}: {1:7.5f}".format("Field Angle",self.theta)
        return s
        

    def _getValue(self,x,y):
        """
        Get the value of phase at specified poistion.

        :param x: x position of Vector2d
        :type x: float or Vector2d
        :param y: y position or None is Vector2d given
        :type y: float or None
        :return: Phase of aberration 

        """
        rSqr = x*x + y*y
        if rSqr > 1.0:
            return float("nan")

        phi = 0.5*self.coef[0]*rSqr + 0.125*self.coef[1]*rSqr*rSqr
        if self.theta != 0.0:
            phi += 0.5*self.coef[2]*y*rSqr*self.theta + \
                   0.5*self.coef[3]*y*y*self.theta*self.theta + \
                   0.25*(self.coef[3] + self.coef[4])*rSqr*self.theta*self.theta + \
                   0.5*self.coef[5]*y*self.theta**3
        return phi

class ZernikeWaveFront(WaveFront):
    """
    Class to hold a zernike wavefront , being a list of optical zernike components. There are also method to evaluate and
    display the expansion.

    :param radius: the radius (Default = 1.0)
    :type radius: float
    :param wavelength: the wavelnegth (Default = 0.55)
    :param \*args: coefficiencs as set of parameters or list, may be blank.
    """
    def __init__(self,radius = 1.0,wavelength = Default,*args):
        """
        Form the Zernike class with the coefficients, coefficeints are in microns.
        """
        WaveFront.__init__(self,radius,wavelength)

        self.coef = []
        for z in args:
            if isinstance(z,list):
                self.coef.extend(z)
            elif isinstance(z,float):
                self.coef.append(z)

    def __str__(self):
        """
        Print out list inclduing the component names.
            """ 
        s =  "r: {0:6.4f} w: {1:6.4f}".format(self.radius,self.wavelength)
        for i in range(len(self.coef)):
            s += "\n{0:s} : {1:8.4e}, ".format(opticalZernikeName(i),self.coef[i])
            
        return s

    def __repr__(self):
        """
        Return repr of class, being class name + str(self)
        """
        return "{0:s} ".format(self.__class__.__name__) + str(self)


    def _getValue(self,x,y):
        """
        Get the value of the Zernike Expansion at location x,y with are assumes to be normalsied

        :param x: x value
        :type x: float
        :param y: y vaue
        :type y: float
        :return: the float value

        Note: if x/y outside range (do outside circle specifed by self.radius) this will return "NaN".
        """
        value = 0.0
        for i,z in enumerate(self.coef):
            value += opticalZernike(z,i,x,y)

        return value
    

class PolynomialWaveFront(WaveFront):
    """
    Polynomial expansion of a wavefront in terms of x/y
    """
    def __init__(self,coeff,radius = 1.0,wavelength = Default):
        """
        Form the Polynomial class with the coefficients, coefficeints are in microns.
        """
        WaveFront.__init__(self,radius,wavelength)
        self.coef = coeff

        self.name = ("bias","xtilt","ytilt","x^2","xy","y^2","x^3","x^2 y","x y^2","y^3"\
                     "x^4","x^3 y","x^2 y^2","x y^3","y^4")

    def __repr__(self):
        """   Detailed of the Zernike in tabular form.
        """
        s = "Polynomial Expansion\n"
        for i in range(0,len(self.coef)):
            s += "{0:<24}: {1:7.5f}\n".format(self.name[i],self.coef[i])
        return s

    def polynomial(self,v,i,x,y):
        """
        Internal function to form one polynomial complonent
        """    

        # Trap trivial case
        if v == 0.0:
            return 0.0

        if x*x + y*y > 1.0:
            return float("nan")

        #         hard code the basic set
        if i == 0:
            return v
        elif i == 1:      # First order
            return x*v
        elif i == 2:
            return y*v
        elif i == 3:      # Second order
            return v*x*x
        elif i == 4:
            return v*x*y
        elif i == 5:
            return v*y*y
        elif i == 6:       # Third order
            return v*x**3
        elif i == 7:
            return v*x*x*y
        elif i == 8:
            return v*x*y*y
        elif i == 9:
            return v*y**3
        elif i == 10:     # Fourth order
            return v*x**4
        elif i == 11:
            return v*x**3*y
        elif i == 12:
            return v*x**2*y**2
        elif i == 13:
            return v*x*y**3
        elif i ==14:
            return v*y**4
        else:
            return float("nan")
        
        
    def getValue(self,x,y):
        """
        Get the value as specified poistion and field angle
        """
        x /= self.maxradius              # Normalise
        y /= self.maxradius

        rSqr = x*x + y*y
        if rSqr > 1.0:
            return float("nan")

        phi = 0.0
        for i in range(0,len(self.coef)):
            phi += self.polynomial(self.coef[i],i,x,y)
            
        return 2.0*math.pi*phi/self.wavelength




class ScalarPSF(object):
    """
    Class to form the Scalar PSF from a wavefront.
    """

    def __init__(self,size = 256):
        """
        Creat the object of specified size
        """
        self.size = int(size)
        self.image = np.zeros((self.size,self.size),dtype = complex)

    def addWaveFront(self,wf):
        """
        Add the wavefront, this will calsulate the ScarPSF and return the image,
        but will also hold in interallally so that it can be rendered
        """
        xr = range(0,self.size)

        for j in xr:
            y = 2.0*wf.maxradius*(j - self.size/2)/self.size
            for i in xr:
                 x = 2.0*wf.maxradius*(i - self.size/2)/self.size
                 if x*x + y*y <= wf.maxradius*wf.maxradius:
                     phi =wf.getValue(x,y) 
                     self.image[i,j] = complex(math.cos(phi),math.sin(phi))

        
        self.image = np.fft.fft2(self.image)
        self.image = np.fft.fftshift(self.image)
        self.image = np.absolute(self.image)
        return self.image

    def draw(self):
        """
        Render the ScalaSPF using imshow with grayscale colour map
        """
        plt.imshow(self.image,cmap=plt.cm.gray)



class Interferometer(object):
    """
    Class to implement and Interferometer to display WaveFront or Phase image held 2-numpy array
    """
    def __init__(self,wf = None ,xtilt = 3.0, ytilt = 0.0, size = 256, type = "twyman"):

        self.size = size
        self.setWaveFront(wf)
        self.setTilt(xtilt,ytilt)
        

    def setWaveFront(self,wf):
        """
        Method to set or reset the wavefront
        """
        if isinstance(wf,WaveFront):
            self.phase = wf.getImage(self.size)
        else:
            self.phase = wf
        

    def setTilt(self,xtilt = None, ytilt = None):
        """
        Sets (or rests the tilts)
        """
        if xtilt != None:
            self.xtilt = float(xtilt)
        if ytilt != None:
            self.ytilt = float(ytilt)


    def draw(self):
        """
        Method to render the fringe image any show via plt.image
        """
        if self.phase.all() == None:
            print("Phase array not defined")
            return
        
        xsize,ysize = self.phase.shape 
        im = np.empty(self.phase.shape,dtype = float) # Make empty array to hold fringe pattern

        xcentre = xsize / 2
        ycentre = ysize /2

        for (i,j),v in np.ndenumerate(self.phase):       # Loop through array
            if not math.isnan(v):                        # Valid points only
                x = (i - xcentre)/xcentre
                y = (j - ycentre)/ycentre
                im[i,j] = 1.0 + math.cos(2.0*math.pi*(x*self.xtilt + y*self.ytilt) + v)

        #        Display in greyscale and with unit entend
        plt.imshow(im,cmap=plt.cm.gray,extent=(-1.0,1.0,-1.0,1.0))

        
        
