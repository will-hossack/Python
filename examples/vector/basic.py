"""
   Basic use of vector3d 

Author: Will Hossack, The Unievsrity of Edinburgh
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

main()

