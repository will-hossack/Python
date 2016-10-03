"""
   Set of classes to implement paraxial matrix optics. These classes
   implement the matrix methods, the paraxial rays are in optics.ray

   Author: Will Hossack, The University of Edinburh
"""

from vector import Vector3d,Angle,Unit3d
import matplotlib.pyplot as plt
import math


class ParaxialMatrix(object):
    """
    Class to implement the base ParaxialMatrix with 4 float components and a thickness.
    """
    def __init__(self, a = 1.0, b = 0.0, c = 0.0, d = 1.0, t = 0.0):
        """
        Constructor with up to 5 parameters,
        param a A elemment or ParaxialMatrix (defult = 1.0)
        param b B element (default = 0.0)
        param c C element (default = 0.0)
        param d D element (default = 1.0)
        param t thickness (defaults = 0.0)

        Defaults to unit matrix of zero thickness.
        """
        if isinstance(a,ParaxialMatrix):
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
        return "[a = {0:7.5f} , b = {1:7.5f} , c = {2:7.5f} , d = {3:7.5f} ] t= {4:7.5f}".\
            format(self.A,self.B,self.C,self.D,self.thickness)

    def __repr__(self):
        """
        Return repr of class, being class name + str(self)
        """
        return "{0:s} ".format(self.__class__) + str(self)


    def copy(self):
        """
        Return copy of current ParaxialMatrix
        """
        return ParaxialMatrix(self)

    def trace(self):
        """
        Return the trace of the matrix. 
        """
        return self.A + self.D

    def determinant(self):
        """
        Return the determinant of the matrix
        """
        return self.A*self.D - self.B*self.C

    def inverse(self):
        """
        Return athe inverse of the matrix as a new ParaxialMatrix, also the thickness will be negated.
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
        Scale the current matrix by linear factor. 
        Element B and thickness scales by a, C by /a, A and D unchanged.
        """
        self.B *= a
        self.C /= a
        self.thickness *= a
        return self

    def backPower(self):
        """
        Get the back power, so power in image space (assume the matrix is for imaging system)
        """
        return -self.C

    def backFocalLength(self):
        """
        Get the back focal length, so focal length in image space (assumes the matrix is for imaging system)
        """
        return 1.0/self.backPower()

    def backFocalPlane(self):
        """
        Get position back Focal Plane relative to the output plane (assumes the matrix is for imaging system)
        """
        return -self.A/self.C

    def backPrincipalPlane(self):
        """
        Get position of back principal plane relative to the output plane. (assumes the matrix is for imaging system)
        """
        return (1.0 - self.A)/self.C

    def frontPower(self):
        """
        Get the front power, so power in object space  (assumes the matrix is for imaging system)
        """
        return self.C/self.determinant()

    def frontFocalLength(self):
        """
        Get the front focal length, so focal length in object space.  (assumes the matrix is for imaging system)
        """
        return 1.0/self.frontPower()

    def frontFocalPlane(self):
        """
        Get position front Focal Plane relative to the input plane.  (assumes the matrix is for imaging system)
        """
        return self.D/self.C

    def frontPrincipalPlane(self):
        """
        Get position of front principal plane relative to the input plane.  (assumes the matrix is for imaging system)
        """
        return (self.D - self.determinant())/self.C


    def __mul__(self,m):
        """
        Method to pre-multiply the current matrix by a another 
        Paraxialmatrix, or if a float it is scaled.
 
        return new ParaxialMatrix

        Note: this is a pre-multiply and NOT a normal matrix multiply.
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
            raise TypeError("ParaxialMatrix * call with unknown type " + str(m))
        
        return ParaxialMatrix(a,b,c,d,t)

        
    def __imul__(self,m):
        """
        Method to pre-multiply the current matrix by a another Paraxialmatrix in place,
        if parameter is float or int then matrix will be scaled.

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

    def __add__(self,d):
        """
        Implement adding a distance d, so will return a new matrix after pre-multiply by propagation
        matrix of distance d
        """
        m = PropagationMatrix(float(d))
        return self*m

    def __iadd__(self, d):
        """
        Implement a propagation a distance d in place by pre-multiply of a propagagtion matrix of distance d.
        param d the distance
        """
        m = PropagationMatrix(float(d))
        self *= m
        return self


class PropagationMatrix(ParaxialMatrix):
    """
    ParaxialMatrix for propagation. 
    """
    #      
    def __init__(self,d):
        """
        Constructor:
        param d float the propagation distance
        """
        ParaxialMatrix.__init__(self,1.0 , d , 0.0 , 1.0, d)


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
        param cr right curature (defaults to none)
        """
        if n == None:                   # Only one parameter
            ParaxialMatrix.__init__(self,1.0 , 0.0 , -1.0/float(f_or_cl) , 1.0, 0.0)
        else:
            a = DielectricMatrix(1.0,n,f_or_cl) 
            b = DielectricMatrix(n,1.0,cr) 
            ParaxialMatrix.__init__(self,a*b)

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
        ParaxialMatrix.__init__(self,a*b*c)    # Set self


class MirrorMatrix(ParaxialMatrix):
    """
    Paraxial Matrix of a mirror specified curfature
    """
    def __init__(self,c):
        """
        Constructor with single parameter, the curvature of the mirror
        param c float the curvature of the mirror.
        """
        ParaxialMatrix.__init__(self,1.0,0.0,-2.0*float(c),1.0,0.0)

class CavityMatrix(ParaxialMatrix):
    """
    Paraxial matrix of a laser cavity with left / right mirrors of specified curvature with
    specified separation. This assume an empty cavity with refractive index of unity.
    """
    def __init__(self,lc,t,rc):
        """
        Constuctor with three paramters
        param lc left mirror curvature
        param t cavity length
        param rc right mirror curvature
        """
        lm = MirrorMatrix(lc)
        d = PropagationMatrix(t)
        rm = MirrorMatrix(rc)
        ParaxialMatrix.__init__(self,d*rm*d*lm)
        self.thickness = 0.0              # Special for cavity
        
        


class ParaxialGroup(ParaxialMatrix):
    """
    Class to represent a ParaxialGroup which extents ParaxialMatrix to include input/output planes
    in global coordintes.
    """
    
    def __init__(self,p = 0.0, matrix = ParaxialMatrix(), in_height = float("inf"), out_height = None):
        """
        Constructor for an ParaxialGroup
        param p float input plane on z-axis, (defaults to 0.0)
        param matrix the ParaxialMatrix (defaults to unity matrix)
        param in_height float the height of the input plane (default float("inf"))
        param out_height float the height of the output plane (default float("inf)
        """
        ParaxialMatrix.__init__(self,matrix)   # Set matrix
        self.input_plane = float(p)             # Input plane
        self.inputPlaneHeight = float(in_height)
        if out_height == None:
            self.outputPlaneHeight = self.inputPlaneHeight
        else:
            self.outputPlaneHeight = float(out_height)


    def __str__(self):
        """
        Implement the str() to print out all information including the underlying matrix.
        """
        return "i : {0:7.5f} hi: {1:7.5f} ho: {2:7.5f} : {3:s}".format(self.input_plane,\
                self.inputPlaneHeight,self.outputPlaneHeight,\
                ParaxialMatrix.__str__(self))


    def __repr__(self):
        """
        Implment repr(), with full class name + str()
        """
        return "{0:s} ".format(self.__class__) + str(self)

     #          Method to make a deep copy of the current Paraxial Group
    def copy(self):
        return ParaxialGroup(self.inputPlane,self,self.inputPlaneHeight,self.outputPlaneHeight)


    def __mul__(self,m):
        """
        Redefine a*b to return ParaxialGroup
        """
        n = ParaxialMatrix.__mul__(self,m)
        return ParaxialGroup(self.inputPlane,n,self.inputPlaneHeight,self.outputPlaneHeight)

    def __add__(self,d):
        """
        Redefine self + d to return ParaxialGroup
        """
        n = ParaxialMatrix.__add__(self,d)
        return ParaxialGroup(self.inputPlane,n,self.inputPlaneHeight,self.outputPlaneHeight)

    #          
    def inputPlane(self):
        """
        Method to get input plane.
        """
        return self.input_plane

    #          
    def outputPlane(self):
        """
        Method to get output plane. 
        """
        return self.input_plane + self.thickness    # always calculate

    #          
    def scale(self,a):
        """
        Method to scale the matrix and plane heights 
        """
        ParaxialMatrix.scale(self,a)
        self.inputPlaneHeight *= a
        self.outputPlaneHeight *= a


    def backFocalPlane(self):
        """
        Method to get the back focal plane in global coodinates
        """
        return self.outputPlane() + ParaxialMatrix.backFocalPlane(self)

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
        return self.outputPlane() + ParaxialMatrix.backPrincipalPlane(self)

    #          
    def frontFocalPlane(self):
        """
        Method to get the front focal plane in global coodinates
        """
        return self.inputPlane() + ParaxialMatrix.frontFocalPlane(self)

    #          
    def frontPrincipalPlane(self):
        """
        Method to get the front principal plane in global coordinates
        """
        return self.inputPlane() + ParaxialMatrix.frontPrincipalPlane(self)

    def frontNodalPoint(self):
        """
        Get Front Modal point, (normall the same as front Prinicpal Plane)
        """
        return self.frontFocalPlane() + self.backFocalLength()

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
    
    
    def planePair(self,mag):
        """
        Calcualte the object image plane pair for specifed magnification in global coordinates.
        param mag float the magnification (note most imaging system mag is -ve)
        return truple of [op,ip] being the location of the object and image places respectively
        """
        f = self.backFocalLength()
        u = f*(1.0 - 1.0/mag)
        v = f*(1.0 - mag)
        op = self.frontPrincipalPlane() - u
        ip = self.backPrincipalPlane() + v
        return [op,ip]


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
            raise TypeError("matrix.ParaxialGroup.pointImage: called with unknown type {0:s}".format(str(op))) 

    def draw(self,legend = False):
        """
        Draw the input/output planes and the 4 cardinal planes.
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

        # plt.plot(zff,yfp,"b",z_ip,yip,"k",zfp,ypp,"r",zbp,ypp,"g",zop,yip,"k",zbf,yfp,"b")
        plt.plot(zff,yfp,"#0000FF",label="Front Focal")
        plt.plot(zbf,yfp,"#00005F",label="Back Focal")
        plt.plot(z_ip,yip,"#000000",label="Input Plane")
        plt.plot(zop,yip,"#505050",label="Output Plane")
        plt.plot(zfp,ypp,"r",label="Front Principal")
        plt.plot(zbp,ypp,"g",label="Back Principal")
        if legend:
            plt.legend(loc="lower right",fontsize="xx-small")
        plt.title("Paraxial Group Diagram")
        plt.xlabel("Optical Axis")

class ParaxialAperture(ParaxialGroup):
    """
    Form  a Paraxial aperture
    """ 
    def __init__(self,p ,h):
        m = ParaxialMatrix()                # Default identity matrix
        ParaxialGroup.__init__(self,m,p,h)  # set matrix, position and input height


class ParaxialThinLens(ParaxialGroup):
    """
    Paraxial group to hold a thin lens
    param p (float) input plane position on optical axis 
    param f_or_cl (float) focal length OR left curvature
    param n refarctive index, may be None
    param cr right curvature
    param radius radius of lens (defaults in inf)
    """
    def __init__(self,p,f_or_cl,n = None, cr = None, radius = float("inf")):
        m = ThinLensMatrix(f_or_cl,n,cr)
        ParaxialGroup.__init__(self,p,m,radius)

class ParaxialThickLens(ParaxialGroup):
    """
    Paraxial Group to hold a Thick Lens
    param p (float) input plane position on optical axis 
    param cl (float) left curvature
    param n refarctive index
    param t thickness
    param rl right curvature
    param radius radius of lens (defaults in inf)
    """
    def __init__(self,p,cl,n,t,cr,radius = float("inf")):
        m = ThickLensMatrix(cl,n,t,cr)
        ParaxialGroup.__init__(self,p,m,radius)


class ParaxialMirror(ParaxialGroup):
    """
    Paraxial Group consisting of a curved mirror
    param p (float) input plane position on optical axis 
    param c curvature
    param radius radius of the mirror (defaults ti inf)
    """
    def __init__(self,p,c,radius = float("inf")):
        m = MirrorMatrix(c)
        ParaxialGroup.__init__(self,p,m,radius)
    

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
        
    #     Method to get input plane (input plane of first group)
    def getInputPlane(self):
        return self[0].inputPlane()
        
    #    Method to get the output plane (output plane of last element)
    def getOutputPlane(self):
        return self[-1].uutputPlane()

    #    Method to get the overall ParaxialGroup of the system.
    def getParaxialGroup(self):
        gr = self[0].copy()         # copy of first element
        for g in self[1:]:          # 
            d = g.inputPlane() - gr.outputPlane() # Distance to input
            gr += d
            gr *= g                  # do mult of matrix
            gr.outputPlaneHeight = g.outputPlaneHeight  # set ouput height

        return gr                   # Return the group


    def draw(self):
        for pg in self:
            pg.draw()

