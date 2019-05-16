"""
   Basic use of Vector3d class

Author: Will Hossack, The University of Edinburgh
"""
import tio as t
import vector as v

def main():

    a = v.Vector3d(1.0,2.0,3.0)      # Vector3d with three components specified.
    b = v.Vector3d([3.0,4.0,5.0])    # Vector3d with components as list
    c = a + b                       # Vector addition
    d = v.Vector3d(c)

    t.tprint("Scalar product of a.b is ",a.dot(b))
    t.tprint("Cross product of b x d is :",b.cross(d))

    u = v.Unit3d(d)
    t.tprint("Direction of d is : ",u)


    m,dir = d.unitPair()

    t.tprint("Modulus is :",m," and direction : ",dir)
main()

