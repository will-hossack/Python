"""
   Basic use of vector3d and tio classes.

Author: Will Hossack, The Unievsrity of Edinburgh
"""
import tio as t
import vector as v

def main():

    a = t.getVector3d("Firt vector")
    b = t.getVector3d("Second Vector")
    t.tprint("vector are ",a,b)
    t.tprint("Sum is ", a+b)
    c = a.cross(b)
    t.tprint("Cross product is ",c)

main()

