"""
Set of classes for optical ray tracing.

Author: Will Hossack, The Univesrity of Edinburgh
"""
import math
import ray as r
import wavelength as wl


#               Class for Paraxial Rays
class ParaxialRay(r.Ray):
    
    
    def __init__(self,height = 0.0, angle = 0.0, plane = 0.0, \
                 wavelength = wl.Default ,intensity = 1.0):
        """
        Constuctor with 5 optional arguments
        param height (defaults to 0.0) height from optical axis
        param angle (defaults to 0.0) angle in radians from optical axis
        param plane (defaults to 0.0) location of plane along optical axis
        param wavelength (defaults to Default) wavelength in microns
        paramintensity (defaults to 1.0) intensity of ray
        """
        r.Ray.__init__(self,wavelength,intensity)   # Set wavelength and intensity
        self.height = float(height)               
        self.angle = float(angle)
        self.plane = float(plane)


    #
    def __repr__(self):
        """
        Implment repr() with fulldetail
        """
        return "paraxial.ParaxialRay({0:7.5f}, {1:7.5f}, {2:7.5f}, {3:7.5f}, {4:7.5f})".\
            format(self.height,self.angle,self.plane,self.wavelength,self.intensity)

    #            
    def copy(self):
        """
        Method to make a deep copy of a current ParaxialRay
        """        
        return ParaxialRay(self.height,self.angle,self.plane,\
                           self.wavelength,self.intensity)
    #           
    def setInvalid(self) :
        """
        Method to set an ParaxialRay to inValid
        """
        self.angle = float("nan")       # set angle to be NaN 
        return self

    #           
    def isValid(self):
        """
        Method to test of a Paraxial Ray is valid
        """
        return not math.isnan(self.angle)

    #            
    def propagate(self, distance):
        """
        Method to propagate a ray a specified distance
        param distance, the distance the ray is propagated
        returns True is sucessful, False is Ray is invalid
        """
        if self.isValid() :
            self.plane += distance                # update plane
            self.height += self.angle*distance    # calculate new height
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
            h = self.height*m.A + self.angle*m.B
            a = self.height*m.C + self.angle*m.D
            p = self.plane + m.thickness
            return ParaxialRay(h,a,p,self.wavelength,self.intensity)
        else:
            return self            

    
    
    def multBy(self,m):
        """
        Method to multiply ParaxialRay by ParaxialMatrix in place.
        """
        if self.isValid():
            h = self.height*m.A + self.angle*m.B
            a = self.height*m.C + self.angle*m.D
            self.height = h
            self.angle = a
            self.plane += m.thickness
            return True
        else:
            return False

    #            Method to propagate a ray throgh and ParaxialGroup or 
    #            list or ParaxialGroups
    #            pg the ParaxialGroup or list of ParaxialGroups
    #            
    def propagateThrough(self, pg):

        if isinstance(pg,list):          #   Is it a list
            for group in pg:             #   Process each element in the list
                b = self.propagateThrough(group)
                if b == False:           #    Trap invalid propagation
                    return False
            return True
            
        #         Its a ParaxialGroup, so progess it
        if self.isValid() :            
            self.propagateToPlane(pg.inputPlane)     #   Propagate to input
            if self.height > pg.inputPlaneHeight:    #   Blocked at input
                self.setInvalid()
                return False
            self.multBy(pg.matrix)                   #   Propagate through matrix
            if self.height > pg.outputPlaneHeight:   #   Blocked at output
                self.setInvalid()
                return False

            return True                              #  Propaged OK
        return False                                 #  Was invalid

    #           Method to locate where ray crosses optical axis is global coordinates
    def crossesZero(self):
        if self.height == 0.0 :
            return self.plane            # already there
        elif self.angle == 0:
            return float('inf')          # Trap infinity
        else:
            return self.plane - self.height/self.angle

    #         Method to locate where two ParaxialRays cross in global coordinates
    #         other is other ParaxialRay
    def crosses(self,other):
        dtheta = other.angle - self.angle
        if dtheta == 0.0:
            return float('inf')        #   Trap infinity where rays are parallel
        else:
            return (self.height - other.height - self.plane*self.angle + other.plane*other.angle)/dtheta


    #         toString method to print out information
    def __str__(self):
        return "ParaxialRay: h: " + str(self.height) + " a: " + str(self.angle) + \
            " p: " + str(self.plane) +  " l: " + str(self.wavelength) + " i: " + \
            str(self.intensity)
        
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

    #      
    def copy(self):
        """
        Method to make a copy of the ParaxialMatrix
        """
        return ParaxialMatrix(self.A,self.B,self.C,self.D,self.thickness)

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
        self.scale(f/self.frontFocalLength())

    #          
    def setBackFocalLength(self, f):
        """
        Method to set the back focal length by scaling 
        """
        self.scale(f/self.backFocalLength())  


    #         
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
    def multBy(self,m):
        """
        Method to pre-multiply the currenr matrix by a another Paraxialmatrix in place
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
    def __str__(self):
        """
        Convert to a string
        """
        return "Matrix : [ " + str(self.A) + " , " + str(self.B) + " , " + str(self.C) + \
            " , " + str(self.D) + " ] t : " + str(self.thickness)

#
#            ParaxialMatrix for a thin lens 
class ThinLensMatrix(ParaxialMatrix):

    #        Constuctor with one or three paramters
    def __init__(self,f_or_cl,n=None,cr=None):
        if n == None:                   # Only one parameter
            f = float(f_or_cl)          # focal length
            ParaxialMatrix.__init__(self,1.0 , 0.0 , -1/f , 1.0, 0.0)
        else:
            cl = float(f_or_cl)         # Three paramters
            a = DielectricMatrix(1.0,n,cl) # Front surface
            b = DielectricMatrix(n,1.0,cr) # Back surface
            s = a.mult(b)                  # Total martix
            ParaxialMatrix.__init__(self,s)

#
#         ParaxialMatrix for a thick lens, 
class ThickLensMatrix(ParaxialMatrix):

    #     Constructor with 4 parameters, all reqired
    #     cl      left curvature
    #     n       refractive index
    #     t       thickness of lens
    #     cr      right curvature
    def __init__(self, cl , n , t , cr):
        a = DielectricMatrix(1.0,n,cl)     # Front surface
        b = PropagationMatrix(t)           # Thickness
        c = DielectricMatrix(n,1.0,cr)     # Back surface
        s = a.mult(b).mult(c)              # Total matrix
        ParaxialMatrix.__init__(self,s)    # Set self

            
            
#
#            ParaxialMatrix for propagation
class PropagationMatrix(ParaxialMatrix):

    #        Constuctor with one pararmeter, the focal length
    def __init__(self,d):
        ParaxialMatrix.__init__(self,1.0 , d , 0.0 , 1.0, d)

#
#             DielectricMatrix for flat or curved interface
class DielectricMatrix(ParaxialMatrix):

    #        Constructor with three parameters
    #        nLeft refractive index on left of interface
    #        nRight refractive index on right of interface
    #        curvature (defaults to 0.0) curvature of interface.
    def __init__(self, nLeft, nRight, curvature = 0.0):
        ParaxialMatrix.__init__(self, 1.0 , 0.0, curvature*(nLeft - nRight)/nRight , nLeft/nRight, 0.0)


#            ParaxialGroup class to represent paraxial group with input plane, ParaxialMatrix
#            and input and outplane heights.
class ParaxialGroup(object):
    #
    #        Conststructor 
    def __init__(self, m, p = 0.0, inHeight = None, outHeight = None):
        self.matrix = m.copy()             # Local copy of Matrix
        self.inputPlane = float(p)         # Input plane
        if inHeight == None:
            self.inputPlaneHeight = float("inf")
        else:
            self.inputPlaneHeight = inHeight
        if outHeight == None:
            self.outputPlaneHeight = float("inf")
        else:
            self.outputPlaneHeight = outHeight


    #          Method to make a deep copy of the current Paraxial Group
    def copy(self):
        return ParaxialGroup(self.matrix,self.inputPlane,self.inputPlaneHeight,self.outputPlaneHeight)

    #          method to get input plane
    def getInputPlane(self):
        return self.inputPlane

    #          Method to get output plane
    def getOutputPlane(self):
        return self.inputPlane + self.matrix.thickness    # always calculate


    #          Method to scale the matrix 
    def scale(self,a):
        self.matrix.scale(a)

    #          Method to get the back power
    def backPower(self):
        return self.matrix.backPower()

    #          Method to get the back focal length
    def backFocalLength(self):
        return self.matrix.backFocalLength()

    #          Method to get the back focal plane in global coodinates
    def backFocalPlane(self):
        return self.getOutputPlane() + self.matrix.backFocalPlane()

    #          Method to get the back principal plane (refaltive to the output plane)
    def backPrincipalPlane(self):
        return self.getOutputPlane() + self.matrix.backPrincipalPlane()


    #          Method to get the front power
    def frontPower(self):
        return self.matrix.frontPower()

    #          Method to get the frontfocal length
    def frontFocalLength(self):
        return self.matrix.frontFocalLength()

    #          Method to get the front focal plane in global coodinates
    def frontFocalPlane(self):
        return self.getInputPlane() + self.matrix.frontFocalPlane()

    #          Method to get the front principal plane in global coordinates
    def frontPrincipalPlane(self):
        return self.getInputPlane() + self.matrix.frontPrincipalPlane()

    #          Method to set the front focal length by scaling matrix
    def setFrontFocalLength(self, f):
        self.matrix.setFrontFocalLength(f)

    #          Method to set the back focal length by scaling matrix
    def setBackFocalLength(self, f):
        self.matrix.setBackFocalLength(f) 

    #          Method to get the image plane for specified object plane using geometric
    #          lens fomula
    def imagePlane(self, op):
        u = self.frontPrincipalPlane() - float(op)      # distance from front principal plane
	v = u/(self.backPower()*u - 1.0)           # distance from back principal plane
	return  self.backPrincipalPlane() + v          # where the image is


    #        Define str method
    def __str__(self):
        return "ParaxialGroup : i : " + str(self.inputPlane) + " " + str(self.matrix)

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

    


        




