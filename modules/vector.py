""" 
Set of classes to support two and three dimensional vector manipulation. 

There are four classes:

Vector2d to implement 2 dimensional vectors.
Vector3d to implement 3 dimensional vectors.
Axis2d to implement a 2 dimensional axis
Axis3d to implement a 3 dimensional axis.

Author:  Will Hossack, The University of Edunburgh.
"""

import math

class Vector2d(object):
    """  
    Class to implement two dimensional vector manipulation. 
    """
    
    def __init__(self,x_or_v = 0.0, y = 0.0):
        """ 
        Constructor.
        param x_or_v float x component (default = 0.0)
        param y float y component (default = 0.0)
        OR
        x_or_v Vector2d both componets copied
        OR
        x_or_v    list, [0] = x, [1] = y
        """
        if isinstance(x_or_v,Vector2d):    # Vector given
            self.x = x_or_v.x
            self.y = x_or_v.y
        elif isinstance(x_or_v,list) or isinstance(x_or_v,tuple):
            self.x = float(x_or_v[0])
            self.y = float(x_or_v[1])
        else:                              # two numbers or nothing given
            self.x = float(x_or_v)         # Force to floats
            self.y = float(y)

    #
    def set(self,x_or_v = 0.0, y = 0.0):
        """
        Method to set componets of the vector in various formats.
        param x_or_v float x component (default = 0.0)
        param y float y component (default = 0.0)
        OR
        x_or_v Vector2d both componets copied
        OR
        x_or_v    list, [0] = x, [1] = y
        """
        if isinstance(x_or_v,Vector2d):    # Vector given
            self.x = x_or_v.x
            self.y = x_or_v.y
        elif isinstance(x_or_v,list) or isinstance(x_or_v,tuple):
            self.x = float(x_or_v[0])
            self.y = float(x_or_v[1])
        else:                              # two numbers or nothing given
            self.x = float(x_or_v)         # Force to floats
            self.y = float(y)
            

    def __str__(self):
        """ 
        Implement str() to return a string with components in 8.4e format.
        """
        return "({0:8.4e} , {1:8.4e})".format(self.x,self.y)

    def __repr__(self):
        """ 
        Implments repr() to return a string with full call.
        """
        return "vector.Vector2d{0:s}".format(str(self))
        
    #
    #
    def __len__(self):
        """
        Implement len() to return 2
        """
        return 2

    #
    def __getitem__(self,key):
        """
        Implement indexing on read using [i] syntax.
        """
        if key == 0:
            return self.x
        elif key == 1:
            return self.y
        else:
            raise IndexError("Vector2d invalid read index of : {0}".format(key))
    #
    def __setitem__(self,key,value):
        """
        Implement indexing on write using [i] syntax
        """
        if key == 0:
            self.x = float(value)
        elif key == 1:
            self.y = float(value)
        else:
            raise IndexError("Vector2d invalid write index of : {0}".format(key))
    #
    #
    def copy(self):
        """
        Return a copy of the current Vector2d.
        """
        return Vector2d(self)

    #        
    #
    def polar(self):
        """
        Return a copy of the current vector in polar (r,theta) form
        """
        r = self.abs()
        if r != 0.0 :
            theta = math.atan2(self.y,self.x)
            return Vector2d(r,theta)
        else:
            return Vector2d()    # Default to zero vector
    #
    #
    def rect(self):
        """
        Return a copy of the current vector in rect (x,y) form assumeing it is in 
        polar (r,theta) form.

        Note: there is NO checking, so if the current Vector is NOT in polar form 
        you will get rubbish.
        """
        x = self.x*math.cos(self.y)
        y = self.x*math.sin(self.y)
        return Vector2d(x,y)
    #
    #
    def getComplex():
        """
        Return a copy of the vector as a complex(x,y)
        """
        return complex(self.x,self.y)
    #
    #
    def absSquare(self):
        """
        Return the absSquare of the vector2d as a float. Does not use pow or **2
        """
        return self.x*self.x + self.y*self.y
    # 
    #
    def absCube(self):
        """
        Return the absCube of the Vector2d as a float defined as abs(x**3) + abs(y**3).
        Does not use pow or **3
        """
        x = abs(self.x)
        y = abs(self.y)
        return x*x*x + y*y*y
    #
    #
    def abs(self):
        """
        Return abs of vectors2d as a float.
        """
        return math.sqrt(self.absSquare())
    #
    #
    def __abs__(self):
        """
        Implement abs() method for vector2d
        """
        return math.sqrt(self.absSquare())

    #
    #
    def normalise(self):
        """
        Method to normalised vector in place. Will return selfso can be used in chain.
        Note: if current vector is abs() = 0, the current vector will be set inValid()
        """       
        a = self.abs()
        if a != 0.0:      # Vector must be zero
            self /= a
        else:
            self.setInvalid()
        return self
    #      
    #
    def setLength(self,d):
        """
        Method to set the length of a vector to specified length. Will return self.
        param d float length vector is set to.
        """
        a = self.abs()
        if a != 0.0:          # is current length not zero 
            self *= (d/a)     # scale by multiply
        return self

    #      
    def absNormalised(self):
        """
        Method to return length f vector and normalsied vector as a pair
        """
        a = self.abs()
        if (a == 0.0):
            return 0.0,Vector2d()
        else:
            n = self / a
            return a,n
    #
    #
    def negate(self):
        """
        Method to negate the current vector in place, will return self to can be used in chain.
        """
        self.x = -self.x
        self.y = -self.y
        return self
    #      
    def __neg__(self):
        """
        Method to implement the __neg__ method to return a new -ve vector.
        Note current not changed.
        """
        return Vector2d(-self.x,-self.y)

    #
    #
    def round(self,figs = 0):
        """
        Method to round the current Vector2d to number of decimal figures.
        param figs number of figures to round to (default is 0)
        """
        self.x = round(self.x, figs)
        self.y = round(self.y, figs)
        return self
    # 
    def __gt__(self,b):
        """
        Implement abs(self) > abs(b)
        """
        return abs(self) > abs(b)
    #
    def __ge__(self,b):
        """
        Implement abs(self> >= abs(b)
        """
        return abs(self) >= abs(b)
    #
    def __lt__(self,b):
        """
        Implement abd(self) < abs(b)
        """
        return abs(self) < abs(b)
    #
    def __le__(self,b):
        """
        Implement abs(self) <= abs(b)
        """
        return abs(self) <= abs(b)
    #
    #      
    def setInvalid(self):
        """
         Method to set to invalid current vector2d to Invalid by setting both compoents 
        to float("nan")
        """
        self.x = float("nan")
        self.y = float("nan")
        return self
    #
    #
    def isValid(self):
        """
        Method to deterime if Vector2d is valid, so .x != Nan
        """
        return not math.isnan(self.x)
    #
    #      
    def __nonzero__(self):
        """
        Implment nonzero for logical test if Valid
        """
        return not math.isnan(self.x)
    #      
    #
    def rotate(self,gamma):
        """
        Method to implementate a rotation in place.
        param gamma rotatian angle.
        """
        cos = math.cos(gamma)
        sin = math.sin(gamma)
        x = self.x*cos + self.y*sin
        y = self.y*cos - self.x*sin
        self.set(x,y)      # Set self
        return self
    #
    #
    def __iadd__(self,v):
        """
        Method to implement the += to add Vector2d to current in place, if v is a 
        Vectord2d the components will be added, while if it is a float it will be 
        added to each component.
        """
        if isinstance(v,Vector2d):
            self.x += v.x
            self.y += v.y
        else:
            self.x += v
            self.y += v
        return self         # Return self

    #
    #
    def __isub__(self,v):
        """
        Method to implement the -= to subtract Vector2d from current in place, if v is a 
        Vectord2d the components will be subtrated, while if it is a float it will be 
        subrtacted from each component.
        """
        if isinstance(v,Vector2d):
            self.x -= v.x
            self.y -= v.y
        else:               # assume is constant and sub the each element
            self.x -= v
            self.y -= v
        return self         # Return self

    #
    #
    def __imul__(self,v):
        """
        Method to implement the *= to multiply Vector2d by current in place, if v is a 
        Vectord3d the component s will be multiplied, while if it is a float it will 
        multiply each component.
        """ 
        if isinstance(v,Vector2d):
            self.x *= v.x
            self.y *= v.y
        else:
            self.x *= v
            self.y *= v
        return self

    #
    #
    def __idiv__(self,v):
        """
        Method to implement the /= to divide current vector in place, if v is a 
        Vectord2d the components will be divided, while if it is a float it will 
        divide each component.
        """ 
        if isinstance(v,Vector2d):
            self.x /= v.x
            self.y /= v.y
        else:
            self.x /= v
            self.y /= v
        return self
    
    #
    def __add__(self, b):
        """
        Method to implments the c = self + b to add 2 Vector2d, of b is a float it
        will be added to each component.
        return new Vector2d
        """
        if isinstance(b,Vector2d):
            return Vector2d(self.x + b.x, self.y + b.y)
        else:
            return Vector2d(self.x + b , self.y + b)
    #
    #
    def __radd__(self, b):
        """
        Method to implments the c = b + self to add 2 Vector2d, of b is a float it will 
        be added to each component.
        return new Vector2d
        """
        if isinstance(b,Vector2d):
            return Vector2d(self.x + b.x, self.y + b.y)
        else:
            return Vector2d(self.x + b , self.y + b)

    #
    #
    def __sub__(self, b):
        """
        Method to implments the c = self - b to add 2 Vector2d, of b is a float it
        will be added to each component.
        return new Vector2d
        """
        if isinstance(b,Vector2d):
            return Vector2d(self.x - b.x, self.y - b.y)
        else:
            return Vector2d(self.x - b , self.y - b)

    #
    #
    def __rsub__(self, b):
        """
        Method to implments the c = b - self to add 2 Vector2d, of b is a float it 
        will be added to each component.
        return new Vector2d
        """
        if isinstance(b,Vector2d):
            return Vector2d(b.x - self.x, b.y - self.y)
        else:
            return Vector2d(b - self.x, b - self.y)
    #
    #
    def __mul__(self,b):
        """
        Method to implemnt c = self * b for b being Vector 2d or float, if float it 
        is appled to each component. 
        returns Vector2d.
        """
        if isinstance(b,Vector2d):
            return Vector2d(self.x * b.x, self.y * b.y)
        else:
            return Vector2d(self.x * b , self.y * b)
    #
    #
    def __rmul__(self,b):
        """
        Method to implemnt c = self * b for b being Vector 2d or float, if float 
        it is appled to each component. 
        returns Vector2d.
        """
        if isinstance(b,Vector2d):
            return Vector2d(self.x * b.x, self.y * b.y)
        else:
            return Vector2d(self.x * b , self.y * b)
    #
    #     
    def __div__(self,b):
        """
        Method to implemnt c = self / b for b being Vector 2d or float, if float it 
        is appled to each component. 
        returns Vector2d.
        """
        if isinstance(b,Vector2d):
            return Vector2d(self.x / b.x, self.y / b.y)
        else:
            return Vector2d(self.x / b , self.y / b)
    #
    #     
    def __rdiv__(self,b):
        """
        Method to implemnt c = b / self for b being Vector 2d or float, if float 
        it is appled to each component. 
        returns Vector2d.
        """
        if isinstance(b,Vector2d):
            return Vector2d(b.x / self.x, b.y / self.y)
        else:
            return Vector2d(b / self.x, b / self.y)
    #
    #
    def dot(self,b):
        """
        Method to get .dot product between current and specified Vector3d
        b second Vector2d
        return float the dot product
        """
        return self.x * b.x + self.y * b.y
    #     
    #
    def distanceSquare(self, b):
        """
        Method to get distanceSquare between two Vector2ds Note does NOT use **2 or pow
        param b the second Vector2d
        return float the square of the distance between vectors.
        """
        dx = b.x - self.x
        dy = b.y - self.y
        return dx*dx + dy*dy
    #
    #
    def distanceCube(self, b):
        """
        Method to get distanceCube between two Vector2ds Note does NOT use **2 or pow
        param b the second Vector2d
        return float the square of the distance between vectors.
        """
        dx = abs(b.x - self.x)
        dy = abs(b.y - self.y)
        return dx*dx*dx + dy*dy*dy
    #
    #
    def distance(self,b):
        """
        Method to get distance between two Vector2ds Note does NOT use **2 or pow
        param b the second Vector2d
        return float the square of the distance between vectors.
        """
        return math.sqrt(self.distanceSquare(b))
    #
    #
    def errorSquare(self,b):
        """
        Method to set the normalsied square error between two vectors
        """
        a = self.absSquare()
        b = b.absSquare()
        ds = self.distanceSquare(b)
        n = a*b                    # Normalisation
        if n != 0.0:
            return ds/math.sqrt(n)
        elif a != 0.0:
            return ds/math.sqrt(a)
        elif b != 0.0:
            return ds/math.sqrt(b)
        else:
            return ds             # which must be zero
    #
    #
    def angleBetween(self,b):
        """
        Method to get the angle between two Vector2d
        param b second Vector2d
        param return float, angle in range -pi/2 and pi
        Note: is abs(current) and abs(b) is zero, zero is retunned
        """
        s = abs(self)*abs(b)
        if s == 0.0:
            return 0.0
        else:
            cq = self.dot(b)/s
            return math.acos(cq)
    #
    #     Method to get the area between two Vectors3d
    #     b second Vector3d
    #     return float, area defined by triangle formed by the two vectors 
    #def areaBetween(self,b):
    #    v = self.cross(b)
    #    return 0.5*abs(v)    
    #
    def inverseSquare(self,b):
        """
        Method to get the vector from current to b scaled to inverse square of
        the distance between them, for implementation of inverse square law forces.
        """
        v = b - self
        d = v.absCube()
        v /= d
        return v
#
#
class Vector3d(object):
    """  Class to implement three  dimensional vector manipulation. """
    #
    #
    def __init__(self,x_or_v = 0.0, y = 0.0, z = 0.0):
        """
        Constructor to create and set vector.
        param x_or_v  float  x component (default = 0.0)
        param  y float y component (default = 0.0)
        param  z float z component (default = 0.0)
        OR
        param x_or_v    Vector3d, all three componets copied
        OR
        param x_or_v  list, [0] = 1, [1] = y, [2] = z
        """
        if isinstance(x_or_v,Vector3d):    # Vector given
            self.x = x_or_v.x
            self.y = x_or_v.y
            self.z = x_or_v.z
        elif isinstance(x_or_v,list) or isinstance(x_or_v,tuple):
            self.x = float(x_or_v[0])
            self.y = float(x_or_v[1])
            self.z = float(x_or_v[2])
        else:                              # three numbers or nothing given
            self.x = float(x_or_v)
            self.y = float(y)
            self.z = float(z)
    # 
    #        
    def set(self,x_or_v = 0.0, y = 0.0,z = 0.0):
        """
        Method to set vector with various augument types.
        param x_or_v  float  x component (default = 0.0)
        param  y float y component (default = 0.0)
        param  z float z component (default = 0.0)
        OR
        param x_or_v    Vector3d, all three componets copied
        OR
        param x_or_v  list, [0] = 1, [1] = y, [2] = z
        """
        if isinstance(x_or_v,Vector3d):    # Vector given
            self.x = x_or_v.x
            self.y = x_or_v.y
            self.z = x_or_v.z
        elif isinstance(x_or_v,list) or isinstance(x_or_v,tuple):
            self.x = float(x_or_v[0])
            self.y = float(x_or_v[1])
            self.z = float(x_or_v[2])
        else:                              # three numbers or nothing given
            self.x = float(x_or_v)
            self.y = float(y)
            self.z = float(z)
            
    #     
    def __str__(self):
        """
        Implement str() with components formatted with 8.4e 
        """
        return "({0:8.4e} , {1:8.4e}, {2:8.4e})".format(self.x,self.y,self.z)
    #        
    def __repr__(self):
        """
        Impment repr() with class name and components formatted with 8.4e
        """
        return "vector.Vector3d{0:s}".format(str(self))
    #
    #        
    def __len__(self):
        """
        Get len() defined at 3
        """
        return 3
    #
    #
    def __getitem__(self,key):
        """
        Implement indexing on read in [i] syntax.
        """
        if key == 0:
            return self.x
        elif key == 1:
            return self.y
        elif key == 2:
            return self.z
        else:
            raise IndexError("Vector3d invalid read index of : {0:d}".format(key))
    #
    #      
    def __setitem__(self,key,value):
        """
        Implement indexing on write in s[i] syntax.
        """
        if key == 0:
            self.x = float(value)
        elif key == 1:
            self.y = float(value)
        elif key == 2:
            self.z = float(value)
        else:
            raise IndexError("Vector3d invalid write index of : {0:d}".format(key))
    #      
    #
    def copy(self):
        """
        Return a copy of the current Vector3d.
        """
        return Vector3d(self)
    #       
    #
    def polar(self):
        """
        Return a copy of current vector in polar (r,theta,psi) form.
        """
        r = self.abs()
        if r != 0.0 :
            theta = math.acos(self.z/r)
            psi = math.atan2(self.y,self.x)
            return Vector3d(r,theta,psi)
        else:
            return Vector3d()    # Default to zero vector
    #
    #
    def rect(self):
        """
        Return a copy of the current vector in rect (x,y,z) form from current assume to 
        be in polar (r,theta,psi) form.
        Note: there is no error checking, so if the current Vector is NOT in polar 
        form you will get rubbish.
        """
        sinTheta = math.sin(self.y)
        x = self.x*sinTheta*math.cos(self.z)
        y = self.x*sinTheta*math.sin(self.z)
        z = self.x*math.cos(self.y)
        return Vector3d(x,y,z)

    #
    #
    def absSquare(self):
        """
        Return the absSquare of the vector3d as a float. Note does NOT use pow or **2
        """
        return self.x*self.x + self.y*self.y + self.z*self.z
    #
    #
    def absCube(self):
        """
        Return the absCube of the Vector3d as a float defined as 
        abs(x**3) + abs(y**3) + abs(z**3). Note does NOT use pow or **3
        """
        x = abs(self.x)
        y = abs(self.y)
        z = abs(self.z)
        return x*x*x + y*y*y + z*z*z
    #
    #
    def abs(self):
        """
        Return the abs of the Vector3d as a float.
        """
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

    def __abs__(self):
        """
        Implement abs() to return the abs length of Vector3d as a float.
        """
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
    #      
    #       
    def normalise(self):
        """
        Method to normalised vector in place.
       
        Note: if abs() = 0, the current vector will be set inValid()
        """ 
        a = self.abs()
        if a != 0.0:      # Vector must be zero
            self /= a
        else:
            self.setInvalid()
        return self
    #      
    #
    def setLength(self,d):
        """
        Method to set the length (or abs) of the current vector to specified length 
        by scaling.
        param d float length vector is set to.
        """
        a = self.abs()
        if a != 0.0:          # is current length not zero 
            self *= (d/a)     # scale by multiply
        return self
    #    
    #
    def absNormalised(self):
        """
        Method to return abs() and normalsied copy of current vector as a 
        list pair (current vector not changed)
        return [a,n] where a abs() and n is normalsied Vector3d
        """
        a = self.abs()
        if (a == 0.0):
            return 0.0,Vector3d()
        else:
            n = self / a
            return a,n
    #
    #
    def negate(self):
        """
        Method to negate the current vector in place.
        """
        self.x = -self.x
        self.y = -self.y
        self.z = -self.z
        return self
    #      
    def __neg__(self):
        """
        Implement the __neg__ method to return a new vector being the -ve of the current. 
        Note current not changed.
        """
        return Vector3d(-self.x,-self.y,-self.z)

    #
    #
    def round(self,figs = 0):
        """
        Method to round the current Vector2d to number of decimal figures.
        param figs number of figures to round to (default is 0)
        """
        self.x = round(self.x, figs)
        self.y = round(self.y, figs)
        self.z = round(self.z, figs)
        return self
    # 
    #
    def __gt__(self,b):
        """
        Implement > to compare abs() of each vectors
        """
        return abs(self) > abs(b)
    #
    #
    def __ge__(self,b):
        """
        Implement >= to compare abs() of each vector
        """
        return abs(self) >= abs(b)
    #
    #
    def __lt__(self,b):
        """
        Implement < to compare abs() of each vector.
        """
        return abs(self) < abs(b)
    #
    #
    def __le__(self,b):
        """
        Implement <= to compare abs() of each vector.
        """
        return abs(self) <= abs(b)

    #
    #      
    def setInvalid(self):
        """ 
        Method to set current vector3d to Invalid by setting all three 
        compoents to float("nan").
        """
        self.x = float("nan")
        self.y = float("nan")
        self.z = float("nan")
        return self
    #
    #
    def isValid(self):
        """
        Method to test if vector is Valid
        returns True for Valid, else False
        """
        return not math.isnan(self.x)
    #
    def __nonzero__(self):
        """
        Implement the logical no-zero test if a vector is valid. True is self.n != Nan
        """
        return not math.isnan(self.x)
    #
    #
    def rotateAboutX(self,alpha):
        """
        Method to implementate rotation about x axis in place.
        param alpha rotatian about x axis in radians
        """
        cos = math.cos(alpha)
        sin = math.sin(alpha)
        y = self.y*cos + self.z*sin
        z = self.z*cos - self.y*sin
        self.y = y                      # Overwrite y
        self.z = z                      # Overwrite z
        return self
    #
    #       
    def rotateAboutY(self,beta):
        """
        Method to implementate a rotation about y axis in place.
        param beta rotatian about y axis in radians.
        """
        cos = math.cos(beta)
        sin = math.sin(beta)
        x = self.x*cos - self.z*sin
        z = self.x*sin + self.z*cos
        self.x = x                      # Overwrite x
        self.z = z                      # Overwrite z
        return self
    #       
    #
    def rotateAboutZ(self,gamma):
        """
        Method to implementate a rotation about z axis in place.
        param gamma rotatian about z axis in radians.
        """
        cos = math.cos(gamma)
        sin = math.sin(gamma)
        x = self.x*cos + self.y*sin
        y = self.y*cos - self.x*sin
        self.x = x                      # Overwrite x
        self.y = y                      # Overwrite y
        return self
    #
    #       General Rotate about x , y, z. The rotation order is x,y then z
    #       alpha rotation ablut x axis
    #       beta rotation about y axis
    #       gamma rotation about z axis
    #
    def rotate(self,alpha,beta,gamma):
        """
        Method to implmeent general Rotate about x , y, z in place. 
        The rotation order is x,y then z.
        param alpha rotation about x axis in radians
        param beta rotation about y axis in radians
        param gamma rotation about z axis in radians
        """
        self.rotateAboutX(alpha)
        self.rotateAboutY(beta)
        self.rotateAboutZ(gamma)
        return self
    #       
    #
    def __iadd__(self,v):
        """
        Method to implement the += to add Vector3d to current in place, if v is a 
        Vectord3d the components will be added, while if it is a float it will be 
        added to each component.
        """
        if isinstance(v,Vector3d):
            self.x += v.x
            self.y += v.y
            self.z += v.z
        else:               # assume is constant and add the each element
            self.x += v
            self.y += v
            self.z += v
        return self         # Return self

    #
    #
    def __isub__(self,v):
        """
        Method to implement the -= to subtract Vector3d from current in place, if v is 
        a Vectord3d the components will be subtrated, while if it is a float it will
        be subrtacted from each component.
        """
        if isinstance(v,Vector3d):
            self.x -= v.x
            self.y -= v.y
            self.z -= v.z
        else:               # assume is constant and sub the each element
            self.x -= v
            self.y -= v
            self.z -= v
        return self         # Return self

    #     
    #
    def __imul__(self,v):
        """
        Method to implement the *= to multiply Vector3d by current in place, if v is a 
        Vectord3d the components will be multiplied, while if it is a float it will 
        multiply each component.
        """ 
        if isinstance(v,Vector3d):
            self.x *= v.x
            self.y *= v.y
            self.z *= v.z
        else:
            self.x *= v
            self.y *= v
            self.z *= v
        return self

    #   
    #
    def __idiv__(self,v):
        """
        Method to implement the /= to divide current vector in place, if v is a 
        Vectord3d the components will be devided, while if it is a float it will 
        divide each component.
        """ 
        if isinstance(v,Vector3d):
            self.x /= v.x
            self.y /= v.y
            self.z /= v.z
        else:
            self.x /= v
            self.y /= v
            self.z /= v
        return self
    #
    #
    def __add__(self, b):
        """
        Implement c = self + b for Vector3d, if b is float, then it will be added to 
        each component.
        returns new Vector3d
        """
        if isinstance(b,Vector3d):
            return Vector3d(self.x + b.x, self.y + b.y, self.z + b.z)
        else:
            return Vector3d(self.x + b , self.y + b, self.z + b)
    #
    #
    def __radd__(self, b):
        """
        Impement c = b + self for b not a Vector3d, (typically a float)
        returns new Vector3d
        """
        if isinstance(b,Vector3d):
            return Vector3d(self.x + b.x, self.y + b.y, self.z + b.z)
        else:
            return Vector3d(self.x + b , self.y + b, self.z + b)

    #
    #
    def __sub__(self, b):
        """
        Implement c = self - b for Vector3d, if b is float it will be subtracted 
        from each element.
        returns new Vector3d
        """
        if isinstance(b,Vector3d):
            return Vector3d(self.x - b.x, self.y - b.y, self.z - b.z)
        else:
            return Vector3d(self.x - b , self.y - b, self.z - b)
    #
    #
    def __rsub__(self, b):
        """
        Implement c = b - self for Vector3d, if b is float it will be subtracted 
        from each element.
        returns new Vector3d
        """
        if isinstance(b,Vector3d):
            return Vector3d(b.x - self.x, b.y - self.y, b.z - self.z)
        else:
            return Vector3d(b - self.x , b - self.y, b - self.z)
    #     
    #
    def __mul__(self,b):
        """
        Implement c = self * b for Vectors3d, if b is float it will multiply each element.
        return new Vector3d
        """
        if isinstance(b,Vector3d):
            return Vector3d(self.x * b.x, self.y * b.y, self.z * b.z)
        else:
            return Vector3d(self.x * b , self.y * b, self.z * b)
    #
    #
    def __rmul__(self,b):
        """
        Implement c = b * self for Vector3d, if b is float it will multiply each element.
        """
        if isinstance(b,Vector3d):
            return Vector3d(self.x * b.x, self.y * b.y, self.z * b.z)
        else:
            return Vector3d(self.x * b , self.y * b, self.z * b)
    #
    #    
    def __div__(self,b):
        """
        Implement c = self / b for Vector3d, if b is float it will divide each element
        """ 
        if isinstance(b,Vector3d):
            return Vector3d(self.x / b.x, self.y / b.y, self.z / b.z)
        else:
            return Vector3d(self.x / b , self.y / b, self.z / b)
    #
    #
    def __rdiv__(self,b):
        """
        Implement c = b / self for Vector3d, if b is float it will divide each element
        """ 
        if isinstance(b,Vector3d):
            return Vector3d(b.x / self.x, b.y / self.y , b.z / self.z )
        else:
            return Vector3d(b / self.x, b / self.y, b / self.z )
    #

    def propagate(self,d,u):
        """
        Return a new vectors that is self + d*u.
        Added for efficency in options calcualations.
        """
        return Vector3d(self.x + d*u.x , self.y + d*u.y, self.z + d*u.z)
    #     
    def dot(self,b):
        """
        Method to form the .dot product of self . b 
        returns float the dot product.
        """
        return self.x * b.x + self.y * b.y + self.z * b.z
    #
    #
    def cross(self,b):
        """
        Method to form the cross product c = self x b
        return Vector3d the cross product
        """
        tx = self.y*b.z - self.z*b.y
        ty = self.z*b.x - self.x*b.z
        tz = self.x*b.y - self.y*b.x
        return Vector3d(tx,ty,tz)

    #
    #
    def distanceSquare(self, b):
        """
        Method to get distanceSquare between two Vector3d, Note does NOT use **2 or pow.
        param b the second Vector3d
        return float the square of the distance between vectors
        """
        dx = b.x - self.x
        dy = b.y - self.y
        dz = b.z - self.z
        return dx*dx + dy*dy + dz*dz
    #
    #     Method to get distanceCube between two Vector3d, 
    #     Note does NOT use **2 or pow
    #
    def distanceCube(self, b):
        """
        Method to get distanceCube between two Vector3d, defined by sum |a.i - b.i|^3. 
        Note does NOT use **2 or pow.
        param b the second Vector3d
        return float the cube of the distance between vectors.
        """
        dx = abs(b.x - self.x)
        dy = abs(b.y - self.y)
        dz = abs(b.z - self.z)
        return dx*dx*dx + dy*dy*dy + dz*dz*dz
    #
    #
    def distance(self,b):
        """
        Method to det the distance between two Vector3d
        param b second Vector3d
        return float distance between two vectors.
        """
        return math.sqrt(self.distanceSquare(b))
    #
    #
    def errorSquare(self,b):
        """
        Method to get the normalsied square error between two Vector3d 
        param b Vector3d, the second vector
        return float the normalsied square error
        """
        a = self.absSquare()
        b = b.absSquare()
        ds = self.distanceSquare(b)     # square distance
        n = a*b                         # Normalisation
        if n != 0.0:
            return ds/math.sqrt(n)      
        elif a != 0.0:
            return ds/math.sqrt(a)
        elif b != 0.0:
            return ds/math.sqrt(b)
        else:
            return ds             # which must be zero
    #     
    #
    def angleBetween(self,b):
        """
        Method to get the angle between two Vector3d
        param b Vector3d second Vector3d
        return float, angle in range -pi/2 and pi
        Note: is abs(current) and abs(b) is zero, zero is returned.
        """
        s = abs(self)*abs(b)
        if s == 0.0:
            return 0.0
        else:
            cq = self.dot(b)/s
            return math.acos(cq)
    #
    #    
    def areaBetween(self,b):
        """
        Method to get the area between two Vectors3d defined by a = 0.5 * abs (self x b)
        param b Vector3d second Vector3d
        return float, area defined by triangle formed by the two vectors 
        """
        v = self.cross(b)
        return 0.5*abs(v)
    #     
    #
    def inverseSquare(self,b):
        """
        Method to get the Vector3d from current to Vector3d b scaled by inverse square of
        the distance between them, for implementation of inverse square law forces.
        Formula implements is        v = (b - self) |b - self |^3
        param b Vector3d, the second vector. 
        returned Vectors3d the scaled vector
        """
        v = b - self
        d = v.absCube()
        v /= d
        return v
#
#        
       
#
class Axis3d(object):
    """
     Axis3d class that define a coordinate axis system with origin and axis vectors
    """
    #
    #
    def __init__(self,origin = Vector3d(), u1_or_axis = Vector3d(1,0,0), \
                 u2 = Vector3d(0,1,0), u3 = Vector3d(0,0,1)):
        """
        Constructor for Axis3d
        param origin Vector3d the nex axis origin (defaults to (0,0,0))
        param u1_or_axis the u1 unit vector, (defaults to (1,0,0))
        param u2 Vector3d u2 unit vector (defaults to (0,1,))
        param u3 Vector3d u3 unit vectors (defaults to (0,0,1)
        OR
        u1_or_axis list[] of three vectors
        
        Note:      supplied u_i vectors will be automatically normalsied
        """
        self.origin = Vector3d(origin)
        self.axis = []
        if (isinstance(u1_or_axis,list)):    # Axis supplied as a list.
            for a in u1_or_axis:
                self.axis.append(Vector3d(a))
        else:                                # Supplied at 3 vectors
            self.axis.append(Vector3d(u1_or_axis))
            self.axis.append(Vector3d(u2))
            self.axis.append(Vector3d(u3))

        for u in self.axis:                 # Normalise axis components
            u.normalise()
    #
    #    
    def __repr__(self):
        """
        Implement the repr() method to show state of axis"
        """
        return "Axis: Orgin at : {0:s}\n{1:s}\n{2:s}\n{3:s}".\
        format(self.origin,self.axis[0],self.axis[1],self.axis[2])

    #
    #
    def transform(self,vec):
        """
        Method to transform a Vectors3d into this axis
        param vec Vector3d to be transformned
        return new Vector3d in the new axis
        """
        v = vec - self.origin     # Shift the origin first
        x = v.dot(self.axis[0])   # x component
        y = v.dot(self.axis[1])   # y component
        z = v.dot(self.axis[2])   # z component
        return Vector3d(x,y,z)
        
        
        

