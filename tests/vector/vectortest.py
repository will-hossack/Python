"""
       Test methods for vector package
"""
import vector as v
import tio

def main():
    a = tio.getVector3d("Give a vector",[1,2,3])
    tio.tprint(repr(a))

    r,u = a.unitPair()

    tio.tprint("Abs is  : ",r," Unit is : ",u)

    u = v.Unit3d(a)
    tio.tprint("Unit vector is " , repr(u))


    ang = v.Angle(u)
    tio.tprint("Angle of u is ", repr(ang))

    ang = v.Angle(a)
    tio.tprint("Angle of a is ", repr(ang))

    w = v.Unit3d(ang)
    tio.tprint("Unit 3d from angle is ",w)

    d = tio.getAngle("Give angle pair",[0.5,0.6])
    tio.tprint(repr(d))
    

main()
