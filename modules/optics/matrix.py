"""
   Set of classes to implement paraxial matrix optics. These classes
   implement the matrix methods, the paraxial rays are in optics.ray.

   These classes are standalone can be used without using the complexities of the full optics package. 
"""

from vector import Vector3d,Angle,Unit3d
from matplotlib.pyplot import plot,legend
import math
import sys
import tio


class ParaxialMatrix(object):
    """
    Class to implement the a ParaxialMatrix with 4 float components and a thickness.

    :param a: A elemment or ParaxialMatrix (defult = 1.0)
    :type a: float or ParaxialMatrix
    :param b: B element (default = 0.0)
    :type b: float
    :param c: C element (default = 0.0)
    :type c: float
    :param d: D element (default = 1.0)
    :type d: float
    :param t: thickness (defaults = 0.0)
    :type t: float

    Defaults to unit matrix of zero thickness.

    """
    def __init__(self, a = 1.0, b = 0.0, c = 0.0, d = 1.0, t = 0.0):
        """
       
        """
        if isinstance(a,ParaxialMatrix):       # Make a copy
            self.A = a.A
            self.B = a.B
            self.C = a.C
            self.D = a.D
            self.thickness = a.thickness
        else:            
            self.A = float(a)
            self.B = float(b)
            self.C = float(c)
            self.D = float(d)
            self.thickness = float(t)

    def __str__(self):
        """
        Return string representation with the four components and thickess displayed in 7.5f format.
        """
        return "m: [{0:7.4f}, {1:7.4f}, {2:7.4f}, {3:7.4f}] t: {4:7.4f}".\
            format(self.A,self.B,self.C,self.D,self.thickness)

    def __repr__(self):
        """
        Return repr of class, being class name + str(self)
        """
        return "{0:s} ".format(self.__class__.__name__) + str(self)

    def copy(self):
        """
        Return copy of current ParaxialMatrix

        :return: Copy of the current Paraxialmatrix

        """
        return ParaxialMatrix(self)

    def trace(self):
        """
        Trace of the matrix. 

        :return: float, being trace of the matrix

        """
        return self.A + self.D

    def determinant(self):
        """
        Determinant of the matrix.

        :return: float, being determinant of matrix.

        """
        return self.A*self.D - self.B*self.C

    def inverse(self):
        """
        Return the inverse of the matrix as a new ParaxialMatrix, also the thickness will be negated.

        :return: Intererse, as a mew Paraxial Matrix

        """
        det = self.determinant()
        a = self.D/det
        b = -self.B/det
        c = -self.C/det
        d = self.A/det
        t = -self.thickness
        return ParaxialMatrix(a,b,c,d,t)

    def scale(self,a):
        """
        Scale the current matrix in place by linear factor. 
        Element B and thickness scales by a, C by /a, A and D unchanged.

        :param a: scale factor:
        :type a: float

        This is equivalent to scaling the optical component that the matrix represents.
        """
        self.B *= a
        self.C /= a
        self.thickness *= a
        return self

    def backPower(self):
        """
        Get the back power, so power in image space (assume the matrix is for imaging system)

        :return: back power as a float

        """
        return -self.C

    def backFocalLength(self):
        """
        Get the back focal length, so focal length in image space (assumes the matrix is for imaging system)

        :return: back focal length as a float.
        """
        return 1.0/self.backPower()

    def backFocalPlane(self):
        """
        Get position back Focal Plane relative to the output plane (assumes the matrix is for imaging system)

        :return: position of back focal plane relatative to output plane.

        """
        return -self.A/self.C

    def backPrincipalPlane(self):
        """
        Get position of back principal plane relative to the output plane. 
        (assumes the matrix is for imaging system)
        
        :return: position of back principle plane relative to output plane.

        """
        return (1.0 - self.A)/self.C

    def frontPower(self):
        """
        Get the front power, so power in object space.

        :return: front power as a float

        """
        return self.C/self.determinant()

    def frontFocalLength(self):
        """
        Get the front focal length, so focal length in object space. 

        :treturn: front focal length as a float.

        """
        return 1.0/self.frontPower()

    def frontFocalPlane(self):
        """
        Get position front Focal Plane relative to the input plane.

        :return: position of front focal plane relative to input plane as a float

        """
        return self.D/self.C

    def frontPrincipalPlane(self):
        """
        Get position of front principal plane relative to the input plane.  

        :return: position of front principal plane relative to input plane as a float.

        """
        return (self.D - self.determinant())/self.C

    def setFocalLength(self,focal):
        """
        Scale the matrix to have a specified back focal length.

        :param focal: target focal length.
        :type flocal: float
        :return: self
        
        """
        f = self.backFocalLength()
        self.scale(focal/f)
        return self

    def __imul__(self,m):
        """
        Method to pre-multiply the current matrix by a another Paraxialmatrix in place, giving
        self \*= m. If m is a float or int, it will be scaled.

        :param m: the ParaxialMatrix
        :type m: ParaxialMatrix or float\
        :return: self 

        Note: this is a pre-multiply and not the normal matrix multiply.
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

    def __mul__(self,m):
        """
        Method to pre-multiply the current matrix by a another 
        Paraxialmatrix, or if a float it is scaled. Implememnts a = self \* m
 
        :param m: ParaxialMatrix to premultiply
        :type m: ParaxialMatrix or float
        :return: a  new ParaxialMatrix

        Note: this is a pre-multiply and NOT a normal matrix multiply.
        """
        r = self.copy()
        r *= m
        return r
        
    
    def __add__(self,d):
        """
        Implement adding a distance d, so will return a new matrix after pre-multiply by propagation
        matrix of distance d

        :param: d the distance.
        :type d: float or int
        :return: ParaxialMatrix

        """
        r = self.copy()
        r += d
        return r

    def __iadd__(self, d):
        """
        Implement a propagation a distance d in place by pre-multiply of a propagation matrix of distance d.
        
        :param d: the distance
        :type d: float or int
        :return: self

        """
        m = PropagationMatrix(d)
        self *= m
        return self


class PropagationMatrix(ParaxialMatrix):
    """
    ParaxialMatrix for propagation a specifed distance

    :param d: the propagation distance.
    :type d: float

    """
    #      
    def __init__(self,d):
        """
        Constructor with single parameter.
        """
        ParaxialMatrix.__init__(self,1.0 , float(d) , 0.0 , 1.0, d)


class DielectricMatrix(ParaxialMatrix):
    """
    Matrix for flat or curved interface dilectric interface.

    :param nl:  refractive index on left of interface.
    :type nl: float
    :param nr: refractive index on right of interface
    :type nl: float
    :param c: curvature of interface. (default = 0.0 for flat surface)
    :type c: float
    
    """
    #       
    def __init__(self, nl, nr, c = 0.0):
        """
        Constructor with three parameters
        
        """
        ParaxialMatrix.__init__(self, 1.0 , 0.0, c*(nl - nr)/nr , nl/nr, 0.0)

class ThinLensMatrix(ParaxialMatrix):
    """
    ParaxialMatrix for a thin lens, will take ether one parameter (focal length) or 
    three paraeters, being curvatues and refratcive index. 

    :param f_or_cl: flocal length or left curvature.
    :type f_or_cl: float
    :param n: refractive index (defaults None), if None, then assume focal length given
    :type n: float
    :param cr: right curature (defaults to None), if n is None, this is not accessed.
    :type cr: float

    """
    #       
    def __init__(self,f_or_cl,n = None,cr = None):
        """
        Constuctor with one or three paramters
        """
        if n == None:                   # Only one parameter
            ParaxialMatrix.__init__(self,1.0 , 0.0 , -1.0/float(f_or_cl) , 1.0, 0.0)
        else:
            a = DielectricMatrix(1.0,n,f_or_cl) 
            b = DielectricMatrix(n,1.0,cr) 
            ParaxialMatrix.__init__(self,a*b)

class ThickLensMatrix(ParaxialMatrix):
    """
     ParaxialMatrix for a thick lens with 4 required parameters.

    :param cl: left curvature
    :type cl: float
    :param n: refractive index
    :type n: float
    :param t:  thickness of lens
    :type t: float
    :param cr: right curvature
    :type cr: float
    
    """
    #
    def __init__(self, cl , n , t , cr):
        """
        Constructor with 4 parameters, all reqired
        
        """
        a = DielectricMatrix(1.0,n,cl)     # Front surface
        b = PropagationMatrix(t)           # Thickness
        c = DielectricMatrix(n,1.0,cr)     # Back surface
        ParaxialMatrix.__init__(self,a*b*c)    # Set self


class DoubletMatrix(ParaxialMatrix):
    """
    Paraxial matrix for double lens with common surface.
    
    :param cl: left curvature
    :type cl: float
    :param nl: left refractive index
    :type nl: float
    :param tl: left thickness
    :type tl: float
    :param cm: middle curvature
    :type cm: float
    :param nr: right refractive index
    :type nr: float
    :param tr: right thickness
    :type tr: float
    :param cr: right curvatute
    :type cr: float
    
    """
    def __init__(self, cl, nl, tl, cm, nr, tr, cr):
        """
        Constructor for a doublet, all required
        
        """
        a = DielectricMatrix(1.0,nl,cl)
        a += tl
        a *= DielectricMatrix(nl,nr,cm)
        a += tr
        a *= DielectricMatrix(nr,1.0,cr)
        ParaxialMatrix.__init__(self,a)

class MirrorMatrix(ParaxialMatrix):
    """
    Paraxial Matrix of a mirror specified by curvature.
    
    :param c: the curvature of the mirror.
    :type c: float

    """
    def __init__(self,c):
        """
        Constructor with single parameter, the curvature of the mirror
        
        """
        ParaxialMatrix.__init__(self,1.0,0.0,-2.0*float(c),1.0,0.0)

class CavityMatrix(ParaxialMatrix):
    """
    Paraxial matrix of a laser cavity with left / right mirrors of specified curvature with
    specified separation. This assume an empty cavity with refractive index of unity.

    :param lc: left mirror curvature
    :type lc: float
    :param t: cavity length
    :type t: float
    :param rc: right mirror curvature
    :type rc: float
    
    """
    def __init__(self,lc,t,rc):
        """
        Constuctor with three paramters
        
        """
        lm = MirrorMatrix(lc)
        d = PropagationMatrix(t)
        rm = MirrorMatrix(rc)
        ParaxialMatrix.__init__(self,d*rm*d*lm)
        self.thickness = 0.0              # Special for cavity
        
        


class ParaxialGroup(ParaxialMatrix):
    """
    Class to represent a ParaxialGroup which extents ParaxialMatrix to include input/output planes
    in global coordinates.

    :param p: input plane on z-axis, (defaults to 0.0)
    :type p: float
    :param matrix: the ParaxialMatrix (defaults to unity matrix)
    :type matrix: ParaxialMatrix
    :param in_height: the height of the input plane (default float("inf"))
    :type in_height: float
    :param out_height: the height of the output plane (default float("inf)
    :type out_height: float
    :param title: title string (defaults to None)
    :type title: str

    """
    def __init__(self,p = 0.0, matrix = ParaxialMatrix(), in_height = float("inf"), out_height = None, title = None):
        """
        Constructor for an ParaxialGroup
        
        """
        ParaxialMatrix.__init__(self,matrix)    # Set matrix
        self.setInputPlane(p)             # Input plane
        self.inputPlaneHeight = float(in_height)
        if out_height == None:
            self.outputPlaneHeight = self.inputPlaneHeight
        else:
            self.outputPlaneHeight = float(out_height)
        self.title = title


    def __str__(self):
        """
        Implement the str() to print out all information including the underlying matrix.
        """
        if self.title != None:
            s = "t : {0:s} ".format(self.title)
        else:
            s = ""
        return s + "i: {0:7.2f} hi: {1:5.3f} ho: {2:5.3f} {3:s}".format(self.input_plane,\
                self.inputPlaneHeight,self.outputPlaneHeight, ParaxialMatrix.__str__(self))

     #          Method to make a deep copy of the current Paraxial Group
    def copy(self):
        return ParaxialGroup(self.input_plane,self,self.inputPlaneHeight,self.outputPlaneHeight,self.title)



    def setInputPlane(self,p):
        """
        Method to set the input plane

        :param p: new input plane
        :type p: float

        """
        self.input_plane = float(p)
        return self

    def incrementInputPlane(self,delta):
        """
        Incremment the input plane
   
        :param delta: the shift
        :type delta: float

        """
        self.input_plane += float(delta)
        return self

    def __imul__(self,m):
        """
        Redefine self \*= m to deal with m being a ParaxialGroup
        """
        if isinstance(m,ParaxialGroup):
            d = m.inputPlane() - self.outputPlane()
            ParaxialMatrix.__iadd__(self,d)
            self.outputPlaneHeight = m.outputPlaneHeight
            
        ParaxialMatrix.__imul__(self,m)
        return self

    #          
    def inputPlane(self):
        """
        Method to get input plane.
        
        :return: potiton of the input plane in global coordinates.

        """
        return self.input_plane

    #          
    def outputPlane(self):
        """
        Method to get output plane. 

        :return: position of the output plane is global coordinates.

        """
        return self.input_plane + self.thickness    # always calculate


    def maxRadius(self):
        """
        Method to get the maximum radius, typically the inputPlaneHeight

        :return: maximum radius as a float.

        Typically called by other classes.

        """
        return self.inputPlaneHeight

    def getPoint(self):
        """
        Method to get the group point, used by other optics classes.

        :return: get the group point Vector3d

        For compatibility for OpticalGroup.

        """
        return Vector3d(0.0,0.0,self.input_plane)
          
    def scale(self,a):
        """
        Method to scale the matrix and plane heights.

        :param a: scale factor
        :type a: float
        
        """
        ParaxialMatrix.scale(self,a)
        self.inputPlaneHeight *= a
        self.outputPlaneHeight *= a


    def backFocalPlane(self):
        """
        Method to get the back focal plane in global coodinates

        :return: location of front focal plane in global coordinates.

        """
        return self.outputPlane() + ParaxialMatrix.backFocalPlane(self)

    def backNodalPoint(self):
        """
        Get the back Nodal point, (normally same as back principal plane)

        :return: location of back focal point in global coordinates.

        """
        return self.backFocalPlane() + self.frontFocalLength()


    def backPrincipalPlane(self):
        """
        Method to get the back principal plane 

        :return: location of back principle plane is global coordinates.

        """
        return self.outputPlane() + ParaxialMatrix.backPrincipalPlane(self)

    #          
    def frontFocalPlane(self):
        """
        Method to get the front focal plane in global coodinates

        :return: location of front focal plane in global coordinates.

        """
        return self.inputPlane() + ParaxialMatrix.frontFocalPlane(self)

    #          
    def frontPrincipalPlane(self):
        """
        Method to get the front principal plane in global coordinates
    
        :return: location of front principal plane in global coordinates.

        """
        return self.inputPlane() + ParaxialMatrix.frontPrincipalPlane(self)

    def frontNodalPoint(self):
        """
        Get Front Modal point, (normall the same as front Prinicpal Plane)

        :return: location of front nodal point in global coordinates.

        """
        return self.frontFocalPlane() + self.backFocalLength()

    def cardinalPoints(self):
        """ 
        Method to get the 6 cardinal points as a list.

        :return: [Front Focal Point, Back Focal Point, Front principal plane, Back principal plane,, Front nodal point,  Back nodal point]


        """
        return [self.frontFocalPlane(),self.backFocalPlane(),\
                self.frontPrincipalPlane(),self.backPrincipalPlane(),\
                self.frontNodalPoint(),self.backNodalPoint()]

    def getInfo(self):
        """
        Get detailed infor of the Paraxial Group and formatted string

        :return: formatted string with details of ParaxialGroup.

        """
        pt = self.cardinalPoints()
        return repr(self) + "\nfl: {0:7.5f}".format(self.backFocalLength()) + \
            "\nffp: {0:7.5f}\nbfp: {1:7.5f}\nfpp: {2:7.5f}\nbpp: {3:7.5f}\nfnp: {4:7.5f}\nbnp: {5:7.5f}".\
            format(pt[0],pt[1],pt[2],pt[3],pt[4],pt[5])
        

    
    def imagePlane(self, op):
        """
        Method to get the image plane for specified object plane using geometric lens fomula
        and properties of the current group.
        
        :param op: location on object plane on optical axis
        :type op: float
        :return: position of image plane in global coordinates.

        """
        u = self.frontPrincipalPlane() - float(op)      # distance from front principal plane
        v = u/(self.backPower()*u - 1.0)           # distance from back principal plane
        return  self.backPrincipalPlane() + v          # where the image is
    
    
    def planePair(self,height,mag):
        """
        Calcualte the object image plane pair for specifed magnification in global coordinates.

        :param height: height of object plane
        :type height: float
        :param mag: the magnification (note most imaging system mag is -ve)
        :type mag: float
        :return: list of [obj,ima] being the Paraxial Plane of the object and image respectively.

        """
        f = self.backFocalLength()
        u = f*(1.0 - 1.0/mag)
        v = f*(1.0 - mag)
        op = self.frontPrincipalPlane() - u     # Position of object plane
        ip = self.backPrincipalPlane() + v      # Position of image plane
        obj = ParaxialPlane(op,height)          # Make the planes
        ima = ParaxialPlane(ip,abs(height*mag))
        return [obj,ima]
    

    def setWithPlanes(self,obj,ima):
        """
        Set the position of the ParaxialGroup to match the supplied object and image
        planes. Also if tghe object plane has a sopecified height the image plane height
        will also be set to match the magmification.

        :param obj: The object plane
        :type obj: ParaxialPlane
        :param ima: The image plane
        :type ima: ParaxialPlane
        :return: the magnification

        """
        f = self.backFocalLength()
        d = ima.inputPlane() - obj.inputPlane()         # Distace between planes
        pt = self.backPrincipalPlane() - self.frontPrincipalPlane() # between pp

        alpha = 1.0/(pt - d)
        a = alpha/f
        b = 1/f
        if a == 0:          # Special case of infinite object plane
            v = f
            mag = 0.0
        else:
            v = (-b + math.sqrt(b*b + 4*a))/(2*a) # Positive root of quadratic
            u = 1.0/alpha + v
            mag = v/u
            ima.inputPlaneHeight = abs(mag*obj.inputPlaneHeight)  # Set image plane height

        np = ima.inputPlane() - v - self.backPrincipalPlane()      # Amount to move lens
        self.incrementInputPlane(np)
        return mag

    def imagePoint(self,op):
        """
        Method to calcualte three-dimensional image of a point in object space in global 
        coordinates using geometric optics.
        
        :param op: Position of object point
        :type op: vector.Vector3d OR  vector.Unit3d or vector. Angle where it will assume an finite object
        :return: Position the image point as vector.Vector3d

        """

        if isinstance(op,float):
            return self.imagePoint(Unit3d(Angle(op)))
        
        elif isinstance(op,Unit3d):                         # Infinte object
            p = Vector3d(0,0,self.backNodalPoint())
            return Vector3d(p + op*(self.backFocalLength()/op.z))
        
        elif isinstance(op,Angle):                        # Also infinite object
            return self.imagePoint(Unit3d(op))

        elif isinstance(op,Vector3d):                       # Finite object
            u = self.frontPrincipalPlane() - op.z
            v = u/(self.backPower()*u - 1)
            ip = self.backPrincipalPlane() + v               # Image location
            mag = -v/u
            return Vector3d(mag*op.x , mag*op.y , ip)
       
        else:
            raise TypeError("matrix.ParaxialGroup.pointImage: called with unknown type {0:s}".format(str(op))) 

    def draw(self,showlegend = False):
        """
        Draw the input/output planes and the 4 cardinal planes using plt.plot() in the current axis.

        :param showlegend: flag to display legend on plot, (Default = False)
        :type showlegend: Bool

        """
        if math.isinf(self.inputPlaneHeight) :
            height = 10.0
        else:
            height = self.inputPlaneHeight
        
        y = [-self.inputPlaneHeight,self.inputPlaneHeight]        # The plane heights
        
        #                             Front Focal plane
        ff = self.frontFocalPlane()       
        zff = [ff,ff]
        yfp = [-0.5*height,0.5*height]
        
        #                             Input plane
        ip = self.inputPlane()
        z_ip = [ip,ip]
        yip = [-height,height]
        
        #                             Front Principal
        fp = self.frontPrincipalPlane()
        zfp = [fp,fp]
        ypp = [-1.1*height,1.1*height]
        #                            Back Principle
        bp = self.backPrincipalPlane()
        zbp = [bp,bp]

        #                            output plane
        op = self.outputPlane()
        zop = [op,op]
        #                            Back focal
        bf = self.backFocalPlane()
        zbf =  [ bf,bf]

        #          Do the actual plot by drawing vertical lines.
        plot(zff,yfp,"#0000FF",label="Front Focal")
        plot(zbf,yfp,"#00005F",label="Back Focal")
        plot(z_ip,yip,"#000000",label="Input Plane")
        plot(zop,yip,"#505050",label="Output Plane")
        plot(zfp,ypp,"r",label="Front Principal")
        plot(zbp,ypp,"g",label="Back Principal")
        if showlegend:
            legend(loc="lower right",fontsize="xx-small")



class ParaxialAperture(ParaxialGroup):
    """
    Form  a Paraxial aperture
    """ 
    def __init__(self,p ,h):
        m = ParaxialMatrix()                # Default identity matrix
        ParaxialGroup.__init__(self,m,p,h)  # set matrix, position and input height


class ParaxialThinLens(ParaxialGroup):
    """
    Paraxial group holding a thin lens.
    
    :param p: input plane position on optical axis 
    :type p: float
    :param f_or_cl: focal length OR left curvature.
    :type f_or_cl: float
    :param n: refarctive index (Default to None), if None then f_or_cl is focal lenth
    :type n: float
    :param cr: right curvature (Default = None)
    :type cr: float
    :param radius: radius of lens (Defaults to inf)
    :type radius: float
    :param title: the title (Default = None)
    :type title: str

    """
    def __init__(self,p,f_or_cl,n = None, cr = None, radius = float("inf"),title = None):
        m = ThinLensMatrix(f_or_cl,n,cr)
        ParaxialGroup.__init__(self,p,m,radius,title = title)

class ParaxialThickLens(ParaxialGroup):
    """
    Paraxial Group to hold a Thick Lens
    
    :param p: input plane position on optical axis 
    :type p: float
    :param cl: left curvature of lens
    :type cl: float
    :param n: refarctive index
    :type n: float
    :param t: thickness of lens
    :type t: float
    :param rl:  right curvature of lens
    :type rl: float
    :param radius:  radius of lens (Defaults in inf)
    :type radius: float
    :param title: the title (Defaults to None)
    :type title: str

    """
    def __init__(self,p,cl,n,t,cr,radius = float("inf"),title = None):
        m = ThickLensMatrix(cl,n,t,cr)
        ParaxialGroup.__init__(self,p,m,radius,title = title)

class ParaxialDoublet(ParaxialGroup):
    """
    Paraxial group to hold a doublet lens.
    
    :param p: input plane position on optical axis
    :type p: float
    :param cl: left curvature
    :type cl: float
    :param nl: left refractive index
    :type nl: float
    :param tl: left thickness
    :type tl: float
    :param cm: middle curvature
    :type cm: float
    :param nr: right refractive index
    :type nr: float
    :param tr: right thickness
    :type tr: floar
    :param cr: right curvatute
    :type cr: float
    :param radius: radius of lens (Default = inf)
    :type radius: float
    :param title: the title (Default = Nond)
    :type title: str
    

    """
    def __init__(self, p, cl, nl, tl, cm, nr, tr, cr,radius = float("inf"),title = None):
        """
        Constructor for a doublet, all required
       
        """
        m = DoubletMatrix(cl,nl,tl,cm,nr,tr,cr)
        ParaxialGroup.__init__(self,p,m,radius,title = title)
    


class ParaxialMirror(ParaxialGroup):
    """
    Paraxial Group consisting of a curved mirror
    
    :param p: input plane position on optical axis 
    :type p: float
    :param c: curvature of mirror
    :type c: float
    :param radius: radius of the mirror (Defaults to inf)
    :type radius: float
    :param title: the title (Defaults is None)
    :type title: str

    """
    def __init__(self,p,c,radius = float("inf"),title = None):
        m = MirrorMatrix(c)
        ParaxialGroup.__init__(self,p,m,radius,title = title)


class DataBaseMatrix(ParaxialGroup):
    """
    Class to read a ParaxialGroup from a matrix file.

    :param file: the file or filename to read ParaxialGroup from (Default = None)
    :type file: str or file

    If file = None, user will be prompted via tio.openFile

    """
    def __init__(self, file = None):
        """
        Read a paraxial group from a file, if a str, then the file is opened.
        """
        ParaxialGroup.__init__(self,0.0)

        if file == None:                  #    No file given
            file = tio.openFile("Matrix file","r","matrix")
        elif isinstance(file,str):  #    File name given as a string
            file = tio.getExpandedFilename(file)
            if not file.endswith("matrix"):
                file += ".matrix"
            try:
                file = open(file,"r")
            except:
                print("ParaxialGroup.readFile() failed to open file " + file)
                file = tio.openFile("Matrix file","r","matrix")

        
        fno = None
        #
        #     Read in line at a time
        try:
            for line in file.readlines() :
                line = line.strip()
                if not line.startswith("#") and len(line) > 0:  #     Kill comments and blank lines
                    token = line.split()                        #     Split to token 
 
                    if token[0].startswith("matrix"):           # a matrix
                        a = float(token[1])
                        b = float(token[2])
                        c = float(token[3])
                        d = float(token[4])
                        th = float(token[5])
                        self *= ParaxialMatrix(a,b,c,d,th)

                    elif token[0].startswith("prop"):          # propagation distance 
                        self +=float(token[1])

                    elif token[0].startswith("dielectric"):    # dialetric
                        nl = float(token[1])
                        nr = float(token[2])
                        if len(token) == 4:
                            c = float(token[3])
                        else:
                            c = 0.0
                        self *= DielectricMatrix(nl,nr,c)

                    elif token[0].startswith("thinlens"):      # Thin lens
                        if len(token) == 2:                    # flocal lenth
                            focal = float(token[1])
                            m = ThinLensMatrix(focal)
                        else:
                            lc = float(token[1])               # three parameter lens
                            n = float(token[2])
                            rc = float(token[3])
                            m = ThinLensMatrix(lc,n,rc)
                        self *= m

                    elif token[0].startswith("thicklens"):     # Thick lens
                        lc = float(token[1])
                        n = float(token[2])
                        th = float(token[3])
                        rc = float(token[4])
                        m = ThickLensMatrix(lc,n,th,rc)
                        self *= m

                    elif token[0].startswith("doublet"):      # Doublet
                        lc = float(token[1])
                        nl = float(token[2])
                        tl = float(token[3])
                        mc = float(token[4])
                        nr = float(token[5])
                        tr = float(token[6])
                        rc = float(token[7])
                        m = DoubletMatrix(lc,nl,tl,mc,nr,tr,rc)
                        self *= m

                    elif token[0].startswith("mirror"):       # Mirror
                        lc = float(token[1])
                        m = MirrorMatrix(lc)
                        self *= m

                    elif token[0].startswith("cavity"):        #  Cavity
                        lc = float(token[1])
                        th = float(token[2])
                        rc = float(token[3])
                        m = CavityMatrix(lc,th,rc)
                        self *= m


                    elif token[0].startswith("inputplane"):    #    Deal with input planes and height
                        self.input_plane = float(token[1])

                    elif token[0].startswith("inputheight"):
                        self.inputPlaneHeight = float(token[1])

                    elif token[0].startswith("outputheight"):
                        self.outputPlaneHeight = float(token[1])

                    elif token[0].startswith("height") or token[0].startswith("radius") :
                        self.inputPlaneHeight = float(token[1])
                        self.outputPlaneHeight = self.inputPlaneHeight

                    elif token[0].startswith("fno") :
                        fno = float(token[1])
                        

                    elif token[0].startswith("title"):
                        self.title = ""
                        for t in token[1:]:
                            self.title += str(t) + " "
                                
                        
                    else:
                        print("DataBasematrix: unknown token: " + str(token[0]))
        except:
            print("DataBaseMatrix: failed to read from file on line : " + str(line) + str(sys.exc_info()))
        
        if fno != None:        # Fno has been used
            h = abs(0.5*self.backFocalLength() / fno)
            self.inputPlaneHeight = h
            self.outputPlaneHeight = h



class ParaxialPlane(ParaxialGroup):
    """   
    Class to represent an image / object plane

    :param p: poistion on plane in global coordinates (Default = 0.0)
    :type p: float
    :param h: height of plane (Default = inf)
    :type h: float

    """
    def __init__(self, p = 0.0, h = float("inf")):
        """
        Form a Paraxial plane jist being a unit matrix at a  specfed position
        """
        ParaxialGroup.__init__(self,p,in_height = h)

    def getInfo(self):
        """
        Overload getInfo() since there are no cardinal points

        :return: in formation string

        """
        return repr(self) + "\nPlane position: {0:7.4f}\nHeight: {1:7.4f}".format(self.inputPlane(),self.inputPlaneHeight)

    def draw(self,legend = False):
        """
        Draw the plane  using plt.plot() to the current axis

        :param legend: add legend (Default = False)
        :type legend: Bool

        """
        if math.isinf(self.inputPlaneHeight) :
            height = 10.0
        else:
            height = self.inputPlaneHeight
        
        y = [-height,height]
        ip = self.inputPlane()
        z = [ip,ip]
        plot(z,y,"#000000",label="Plane")

        
        

class ParaxialSystem(list):
    """
    Class to hold a Paraxial sytem being a list of Paraxial Groups.
    """
    
    def __init__(self,*args):
        """
        Constructor to which may have a number of Paraxial Groups, each of
        which are appended.
        """
        list.__init__(self)
        #
        for pg in args:
            if isinstance(pg,list):
                self.extend(pg)
            else:
                self.append(pg)
        
    def getInputPlane(self):
        """
         Method to get input plane (input plane of first group)
        """
        return self[0].inputPlane()

    def getPoint(self):
        """
        Method to get the input point. used by orher optics classes
        """
        return Vector3d(0.0,0.0,self.getInputPlane())

    def maxRadius(self):
        """
        Get the max radius of the first electment
        """
        return self[0].maxRadius()
        
    def getOutputPlane(self):
        """
         Method to get the output plane (output plane of last element)
        """
        return self[-1].outputPlane()

    
    def getParaxialGroup(self):
        """
         Method to get the overall ParaxialGroup of the system.
        """
        gr = self[0].copy()         # copy of first element
        for g in self[1:]:          # 
            d = g.inputPlane() - gr.outputPlane() # Distance to next input place
            gr += d                 # propagate that distance 
            gr *= g                 # do mult of matrix
            gr.outputPlaneHeight = g.outputPlaneHeight  # set ouput height

        return gr                   # Return the group


    def draw(self,legend = False):
        """
        Draw the whole sytem
        """
        for pg in self:
            pg.draw(legend)           # Draw each componet in turn

