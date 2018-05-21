"""
   Set of classes to analyse geometric PSF and produce star test plots
"""
import ray
from wavelength import WavelengthColour
from vector import Vector2d,Vector3d
import matplotlib.pyplot as plt
import math

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
    Class to implment two dimensional moments to second order
    with more efficient code than Moments above.
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
        if isinstance(plane,float):
            pl = plane
        else:
            pl = plane.getPoint().z
        mom = MomentsFixed().addRay(pencil,pl)
        c = mom.centroid()          # Get the centre
        self.set(c.x,c.y,pl)     # Cet position in 3d

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
        

        #            set self with initial condition
        psf = self.setWithRays(pencil,plane)
        area = psf.area()
        zp = psf.z
        wave = pencil[0].wavelength/1000   # Wavelengh in mm
        delta = -0.25
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
    def draw(self,colour = "k" ):
        """
        Draw the psf an ellipse to the curret MatPlotLib axis with 
        default colour "k" (black"
        """
        
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

        plt.plot(xval,yval,colour)

class SpotDiagram(object):
    """
        Create and return a SpotDiagram
    """

    def __init__(self,pencil):
        """
        Constrcutor, takes the specified RayPencil, use the .draw() to
        render the spot diagram
        """
        self.raypencil = pencil       # Record the RayPencil

        
    def draw(self,plane):
        """ 
        Draw the spot disagram to the current active MatPlotLib as circles
        with the centre given my the wavelelength of the first ray.
        """
        xData = []           # X and Y point locations
        yData = []
        
        for r in self.raypencil:
            if r:
                pt = r.pointInPlane(plane)     # Find the point in the plane
                xData.append(pt.x)
                yData.append(pt.y)

        #     Get the colour of the ray as a hex string
        col = WavelengthColour(self.raypencil[0].wavelength).hexString()
        #     Scatter plot to the current figure
        plt.axis('equal')
        plt.scatter(xData,yData,c=col)


